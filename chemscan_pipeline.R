###############################################################################
# ChemScan Pipeline
# Retrieves and integrates chemical-protein interaction data from CTD, IEDB,
# UniProt, PhosphoSitePlus, and Human Protein Atlas to build a prioritised
# PhIP-seq peptide library of chemically-modified human protein epitopes.
#
# Dependencies: tidyverse, httr2, data.table, Biostrings, cli, cachem,
#               BiocManager, jsonlite, digest
###############################################################################

# ── Package loading ──────────────────────────────────────────────────────────

library(tidyverse)
library(httr2)
library(data.table)
library(cli)
library(cachem)
library(jsonlite)
library(digest)

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings", ask = FALSE, update = FALSE)
}
library(Biostrings)

# ── Global constants ─────────────────────────────────────────────────────────

.TILE_LENGTH    <- 62L
.TILE_OVERLAP   <- 14L
.TILE_STEP      <- .TILE_LENGTH - .TILE_OVERLAP  # 48
.CTD_URL        <- "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz"
.UNIPROT_REST   <- "https://rest.uniprot.org"
.IEDB_REST      <- "https://query-api.iedb.org/epitope_search"
.HPA_API        <- "https://www.proteinatlas.org/api/search_download.php"

# Kyte-Doolittle hydropathy values (higher = more hydrophobic)
.KD_SCALE <- c(
  A =  1.8, R = -4.5, N = -3.5, D = -3.5, C =  2.5,
  Q = -3.5, E = -3.5, G = -0.4, H = -3.2, I =  4.5,
  L =  3.8, K = -3.9, M =  1.9, F =  2.8, P = -1.6,

  S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V =  4.2
)

# ── Cache helper ─────────────────────────────────────────────────────────────

#' Initialise a disk-backed cache in the given directory.
#'
#' @param cache_dir Character path to cache directory.
#' @return A cachem disk_cache object.
init_cache <- function(cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cachem::cache_disk(dir = cache_dir, max_age = 60 * 60 * 24 * 7)
}

#' Retrieve a value from cache or compute it.
#'
#' @param cache   A cachem cache object.
#' @param key     Character cache key.
#' @param fn      A zero-argument function that produces the value.
#' @return The cached or freshly computed value.
cache_get_or_set <- function(cache, key, fn) {
  val <- cache$get(key)
  if (!cachem::is.key_missing(val)) {
    cli::cli_alert_info("Cache hit for {.val {key}}")
    return(val)
  }
  val <- fn()
  cache$set(key, val)
  val
}

###############################################################################
# MODULE 1: CTD Chemical-Gene Interactions
###############################################################################

#' Download and filter CTD chemical-gene interactions.
#'
#' @param chemicals Optional character vector of chemical names / CAS numbers
#'   to subset on. NULL returns all human direct interactions.
#' @param cache A cachem cache object.
#' @return A tibble with columns: ChemicalName, CAS_RN, GeneSymbol, UniProtID,
#'   InteractionType, PubMedIDs.
fetch_ctd_interactions <- function(chemicals = NULL, cache = NULL) {

  cli::cli_h2("Module 1: CTD Chemical-Gene Interactions")

  # -- 1a. Download / read the bulk file -------------------------------------
  read_ctd <- function() {
    cli::cli_alert_info("Downloading CTD bulk interaction file ...")
    tmp <- tempfile(fileext = ".tsv.gz")
    tryCatch(
      {
        httr2::request(.CTD_URL) |>
          httr2::req_timeout(600) |>
          httr2::req_perform(path = tmp)
      },
      error = function(e) {
        cli::cli_abort("Failed to download CTD file: {conditionMessage(e)}")
      }
    )

    cli::cli_alert_info("Parsing CTD file ...")
    # The file has 29 comment lines starting with '#'
    raw <- data.table::fread(
      tmp,
      sep       = "\t",
      header    = FALSE,
      skip      = 29,
      quote     = "",
      fill      = TRUE,
      na.strings = ""
    )

    # Assign column names from CTD documentation
    expected_cols <- c(
      "ChemicalName", "ChemicalID", "CAS_RN", "GeneSymbol", "GeneID",
      "GeneForms", "Organism", "OrganismID", "Interaction",
      "InteractionActions", "PubMedIDs"
    )
    if (ncol(raw) >= length(expected_cols)) {
      data.table::setnames(raw, seq_along(expected_cols), expected_cols)
    }

    tibble::as_tibble(raw)
  }

  ctd_raw <- if (!is.null(cache)) {
    cache_get_or_set(cache, "ctd_bulk_raw", read_ctd)
  } else {
    read_ctd()
  }

  # -- 1b. Filter for human, direct molecular interactions -------------------
  cli::cli_alert_info("Filtering for human direct interactions ...")
  ctd <- ctd_raw |>
    dplyr::filter(
      as.character(OrganismID) == "9606",
      stringr::str_detect(
        InteractionActions,
        stringr::regex("binding|modification", ignore_case = TRUE)
      )
    )

  # Remove rows that are only indirect / pathway-level annotations
  if ("Interaction" %in% names(ctd)) {
    ctd <- ctd |>
      dplyr::filter(
        !stringr::str_detect(
          Interaction,
          stringr::regex("pathway|expression|activity", ignore_case = TRUE)
        ) | is.na(Interaction)
      )
  }

  # -- 1c. Optionally subset by user-supplied chemicals ----------------------
  if (!is.null(chemicals)) {
    chemicals_upper <- toupper(chemicals)
    ctd <- ctd |>
      dplyr::filter(
        toupper(ChemicalName) %in% chemicals_upper |
          CAS_RN %in% chemicals
      )
  }

  # -- 1d. Map GeneSymbols -> UniProt IDs ------------------------------------
  genes <- unique(ctd$GeneSymbol)
  cli::cli_alert_info("Mapping {length(genes)} gene symbols to UniProt IDs ...")
  mapping <- map_genes_to_uniprot(genes, cache = cache)

  ctd <- ctd |>
    dplyr::left_join(mapping, by = "GeneSymbol") |>
    dplyr::select(
      ChemicalName, CAS_RN, GeneSymbol, UniProtID,
      InteractionType = InteractionActions, PubMedIDs
    ) |>
    dplyr::distinct()

  cli::cli_alert_success("CTD: {nrow(ctd)} interaction rows retained")
  ctd
}

#' Map a vector of human gene symbols to reviewed UniProt accessions using the
#' UniProt ID-mapping REST API.
#'
#' @param genes Character vector of gene symbols.
#' @param cache Optional cachem cache.
#' @return A tibble with columns GeneSymbol, UniProtID.
map_genes_to_uniprot <- function(genes, cache = NULL) {

  do_mapping <- function() {
    # Process in batches to stay under API limits
    batch_size <- 500L
    batches <- split(genes, ceiling(seq_along(genes) / batch_size))
    results <- vector("list", length(batches))

    for (i in seq_along(batches)) {
      batch <- batches[[i]]
      cli::cli_alert_info(
        "  UniProt ID-mapping batch {i}/{length(batches)} ({length(batch)} genes)"
      )

      # Submit the mapping job
      submit_resp <- tryCatch(
        {
          httr2::request(paste0(.UNIPROT_REST, "/idmapping/run")) |>
            httr2::req_body_form(
              from = "Gene_Name",
              to   = "UniProtKB",
              ids  = paste(batch, collapse = ","),
              taxId = "9606"
            ) |>
            httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
            httr2::req_perform()
        },
        error = function(e) {
          cli::cli_warn("UniProt ID-mapping submit failed: {conditionMessage(e)}")
          return(NULL)
        }
      )

      if (is.null(submit_resp)) next

      job_id <- httr2::resp_body_json(submit_resp)$jobId
      if (is.null(job_id)) next

      # Poll for completion
      result_url <- NULL
      for (attempt in 1:30) {
        Sys.sleep(2)
        status_resp <- tryCatch(
          httr2::request(
            paste0(.UNIPROT_REST, "/idmapping/status/", job_id)
          ) |>
            httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
            httr2::req_perform(),
          error = function(e) NULL
        )
        if (is.null(status_resp)) next
        status_body <- httr2::resp_body_json(status_resp)
        if (!is.null(status_body$results) ||
            !is.null(status_body$failedIds) ||
            isTRUE(status_body$jobStatus == "FINISHED")) {
          result_url <- paste0(
            .UNIPROT_REST, "/idmapping/uniprotkb/results/", job_id,
            "?format=tsv&fields=accession,gene_names,reviewed&size=500"
          )
          break
        }
      }

      if (is.null(result_url)) {
        cli::cli_warn("UniProt ID-mapping timed out for batch {i}")
        next
      }

      # Fetch results
      res_resp <- tryCatch(
        httr2::request(result_url) |>
          httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
          httr2::req_perform(),
        error = function(e) NULL
      )

      if (is.null(res_resp)) next

      res_text <- httr2::resp_body_string(res_resp)
      if (nchar(res_text) < 10) next

      res_tbl <- tryCatch(
        readr::read_tsv(I(res_text), show_col_types = FALSE),
        error = function(e) NULL
      )
      if (is.null(res_tbl) || nrow(res_tbl) == 0) next

      # Prefer reviewed (Swiss-Prot) entries
      res_tbl <- res_tbl |>
        dplyr::rename_with(~ c("From", "Entry", "Gene_Names", "Reviewed")[seq_along(.)]) |>
        dplyr::arrange(dplyr::desc(Reviewed == "reviewed")) |>
        dplyr::distinct(From, .keep_all = TRUE) |>
        dplyr::transmute(
          GeneSymbol = From,
          UniProtID  = Entry
        )

      results[[i]] <- res_tbl
    }

    dplyr::bind_rows(results) |> dplyr::distinct()
  }

  if (!is.null(cache)) {
    key <- paste0("uniprot_gene_map_", digest::digest(sort(genes)))
    cache_get_or_set(cache, key, do_mapping)
  } else {
    do_mapping()
  }
}

###############################################################################
# MODULE 2: IEDB Modified Epitope Retrieval
###############################################################################

#' Query IEDB for linear B-cell epitopes from human proteins bearing chemical
#' modifications.
#'
#' @param cache Optional cachem cache.
#' @return A tibble with columns: EpitopeID, Sequence, ModifiedPositions,
#'   UniProtID, AssayType, AntibodyIsotype, PubMedIDs, HighConfidenceIgG.
fetch_iedb_modified_epitopes <- function(cache = NULL) {

  cli::cli_h2("Module 2: IEDB Modified Epitope Retrieval")

  mod_terms <- c(
    "citrullin", "carbamyl", "acetyl", "oxidis",
    "phospho", "malondialdehyde", "adduct", "hapten"
  )

  fetch_for_term <- function(term) {
    cli::cli_alert_info("  Querying IEDB for term: {.val {term}}")

    resp <- tryCatch(
      {
        httr2::request(.IEDB_REST) |>
          httr2::req_url_query(
            epitope_structure_type = "Linear peptide",
            source_organism_id    = "9606",         # Homo sapiens
            host_organism_id      = "9606",
            epitope_modification  = term,
            response_type         = "B cell",
            output_format         = "json"
          ) |>
          httr2::req_timeout(120) |>
          httr2::req_retry(max_tries = 3, backoff = ~ 3) |>
          httr2::req_perform()
      },
      error = function(e) {
        cli::cli_warn("IEDB query failed for {.val {term}}: {conditionMessage(e)}")
        return(NULL)
      }
    )

    if (is.null(resp)) return(tibble::tibble())

    body <- tryCatch(
      httr2::resp_body_json(resp, simplifyVector = TRUE),
      error = function(e) {
        cli::cli_warn("IEDB JSON parse error for {.val {term}}: {conditionMessage(e)}")
        return(NULL)
      }
    )

    if (is.null(body) || length(body) == 0) return(tibble::tibble())

    # Normalise to tibble
    tbl <- tryCatch(tibble::as_tibble(body), error = function(e) tibble::tibble())
    if (nrow(tbl) == 0) return(tibble::tibble())

    # Extract relevant fields defensively
    safe_col <- function(df, col) {
      if (col %in% names(df)) df[[col]] else NA_character_
    }

    tibble::tibble(
      EpitopeID         = safe_col(tbl, "epitope_id"),
      Sequence          = safe_col(tbl, "linear_sequence"),
      Description       = safe_col(tbl, "description"),
      UniProtID         = safe_col(tbl, "source_antigen_accession"),
      AssayType         = safe_col(tbl, "assay_type"),
      AntibodyIsotype   = safe_col(tbl, "isotype"),
      PubMedIDs         = safe_col(tbl, "pubmed_id"),
      Qualitative       = safe_col(tbl, "qualitative_measure"),
      SearchTerm        = term
    )
  }

  all_epitopes <- if (!is.null(cache)) {
    cache_get_or_set(cache, "iedb_modified_epitopes", function() {
      purrr::map(mod_terms, fetch_for_term) |> dplyr::bind_rows()
    })
  } else {
    purrr::map(mod_terms, fetch_for_term) |> dplyr::bind_rows()
  }

  if (nrow(all_epitopes) == 0) {
    cli::cli_warn("No IEDB epitopes returned; returning empty tibble.")
    return(tibble::tibble(
      EpitopeID = character(), Sequence = character(),
      ModifiedPositions = character(), UniProtID = character(),
      AssayType = character(), AntibodyIsotype = character(),
      PubMedIDs = character(), HighConfidenceIgG = logical()
    ))
  }

  # Deduplicate
  all_epitopes <- all_epitopes |> dplyr::distinct(EpitopeID, .keep_all = TRUE)

  # Flag high-confidence IgG positives
  all_epitopes <- all_epitopes |>
    dplyr::mutate(
      HighConfidenceIgG = stringr::str_detect(
        AntibodyIsotype, stringr::regex("IgG", ignore_case = TRUE)
      ) & stringr::str_detect(
        Qualitative, stringr::regex("Positive", ignore_case = TRUE)
      ),
      HighConfidenceIgG = dplyr::coalesce(HighConfidenceIgG, FALSE),
      # Attempt to parse modification positions from description
      ModifiedPositions = stringr::str_extract_all(
        Description,
        "\\d+"
      ) |> purrr::map_chr(~ paste(.x, collapse = ","))
    )

  result <- all_epitopes |>
    dplyr::select(
      EpitopeID, Sequence, ModifiedPositions, UniProtID,
      AssayType, AntibodyIsotype, PubMedIDs, HighConfidenceIgG
    )

  cli::cli_alert_success("IEDB: {nrow(result)} modified epitopes retrieved")
  result
}

###############################################################################
# MODULE 3: PTM Site Mapping
###############################################################################

#' Retrieve PTM annotations for a set of UniProt accessions from UniProt
#' features API and, optionally, PhosphoSitePlus bulk files.
#'
#' @param uniprot_ids Character vector of UniProt accession IDs.
#' @param ptm_file    Optional path to PhosphoSitePlus bulk download TSV.
#' @param cache       Optional cachem cache.
#' @return A tibble with columns: UniProtID, ResiduePosition, ResidueAA,
#'   ModificationType, EvidenceCount, IsDiseaseSite.
fetch_ptm_sites <- function(uniprot_ids, ptm_file = NULL, cache = NULL) {

  cli::cli_h2("Module 3: PTM Site Mapping")
  uniprot_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])

  if (length(uniprot_ids) == 0) {
    cli::cli_warn("No UniProt IDs supplied to Module 3.")
    return(tibble::tibble(
      UniProtID = character(), ResiduePosition = integer(),
      ResidueAA = character(), ModificationType = character(),
      EvidenceCount = integer(), IsDiseaseSite = logical()
    ))
  }

  # -- 3a. UniProt Features API -----------------------------------------------
  fetch_uniprot_ptms <- function(acc) {
    resp <- tryCatch(
      {
        httr2::request(
          paste0(.UNIPROT_REST, "/uniprotkb/", acc, ".json")
        ) |>
          httr2::req_timeout(30) |>
          httr2::req_retry(max_tries = 3, backoff = ~ 2) |>
          httr2::req_perform()
      },
      error = function(e) {
        cli::cli_warn("UniProt features request failed for {acc}: {conditionMessage(e)}")
        return(NULL)
      }
    )

    if (is.null(resp)) return(tibble::tibble())

    body <- tryCatch(
      httr2::resp_body_json(resp, simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (is.null(body)) return(tibble::tibble())

    features <- body$features
    if (is.null(features) || length(features) == 0) return(tibble::tibble())

    ptm_types <- c("Modified residue", "Cross-link", "Lipidation")
    ptm_features <- purrr::keep(features, function(f) {
      !is.null(f$type) && f$type %in% ptm_types
    })

    if (length(ptm_features) == 0) return(tibble::tibble())

    purrr::map_dfr(ptm_features, function(f) {
      pos <- as.integer(f$location$start$value %||% NA)
      desc <- f$description %||% NA_character_
      # Check for disease annotation in evidences
      evidence_tags <- purrr::map_chr(
        f$evidences %||% list(),
        ~ .x$source$name %||% ""
      )
      is_disease <- any(stringr::str_detect(
        evidence_tags, stringr::regex("disease|patholog", ignore_case = TRUE)
      ))
      evidence_count <- length(f$evidences %||% list())

      tibble::tibble(
        UniProtID        = acc,
        ResiduePosition  = pos,
        ResidueAA        = NA_character_,
        ModificationType = desc,
        EvidenceCount    = as.integer(max(evidence_count, 1L)),
        IsDiseaseSite    = is_disease
      )
    })
  }

  cli::cli_alert_info("Fetching PTM features from UniProt for {length(uniprot_ids)} proteins ...")

  # Process with progress bar
  ptm_uniprot <- if (!is.null(cache)) {
    cache_get_or_set(
      cache,
      paste0("uniprot_ptms_", digest::digest(sort(uniprot_ids))),
      function() {
        cli::cli_progress_bar("UniProt PTMs", total = length(uniprot_ids))
        results <- purrr::map(uniprot_ids, function(acc) {
          res <- fetch_uniprot_ptms(acc)
          cli::cli_progress_update()
          Sys.sleep(0.15)
          res
        })
        cli::cli_progress_done()
        dplyr::bind_rows(results)
      }
    )
  } else {
    cli::cli_progress_bar("UniProt PTMs", total = length(uniprot_ids))
    results <- purrr::map(uniprot_ids, function(acc) {
      res <- fetch_uniprot_ptms(acc)
      cli::cli_progress_update()
      Sys.sleep(0.15)
      res
    })
    cli::cli_progress_done()
    dplyr::bind_rows(results)
  }

  # -- 3b. PhosphoSitePlus bulk file (optional) --------------------------------
  ptm_psp <- tibble::tibble()
  if (!is.null(ptm_file) && file.exists(ptm_file)) {
    cli::cli_alert_info("Parsing PhosphoSitePlus file: {.file {ptm_file}}")
    psp <- tryCatch(
      {
        data.table::fread(ptm_file, sep = "\t", header = TRUE,
                          skip = 3, fill = TRUE, quote = "") |>
          tibble::as_tibble()
      },
      error = function(e) {
        cli::cli_warn("Failed to parse PSP file: {conditionMessage(e)}")
        tibble::tibble()
      }
    )

    if (nrow(psp) > 0) {
      # PhosphoSitePlus column names vary; attempt common schema
      # Expected: ACC_ID, MOD_RSD, ORGANISM, HTP_SITE, LTP_SITE
      acc_col <- intersect(names(psp), c("ACC_ID", "ACC#", "ACCESSION"))[1]
      mod_col <- intersect(names(psp), c("MOD_RSD", "MOD_TYPE"))[1]
      org_col <- intersect(names(psp), c("ORGANISM", "ORG"))[1]
      htp_col <- intersect(names(psp), c("MS_LIT", "HTP_LIT", "HTP#"))[1]
      ltp_col <- intersect(names(psp), c("LT_LIT", "LTP_LIT", "LTP#"))[1]

      if (!is.na(acc_col) && !is.na(mod_col)) {
        psp_filtered <- psp

        if (!is.na(org_col)) {
          psp_filtered <- psp_filtered |>
            dplyr::filter(
              stringr::str_detect(.data[[org_col]], stringr::regex("human", ignore_case = TRUE))
            )
        }

        # Parse residue position and amino acid from MOD_RSD (e.g., "S473-p")
        psp_filtered <- psp_filtered |>
          dplyr::mutate(
            UniProtID       = .data[[acc_col]],
            .rsd_raw        = .data[[mod_col]],
            ResidueAA       = stringr::str_extract(.rsd_raw, "^[A-Z]"),
            ResiduePosition = as.integer(stringr::str_extract(.rsd_raw, "\\d+")),
            ModificationType = stringr::str_extract(.rsd_raw, "(?<=-)\\w+$")
          )

        # Compute evidence count
        if (!is.na(htp_col) && !is.na(ltp_col)) {
          psp_filtered <- psp_filtered |>
            dplyr::mutate(
              EvidenceCount = as.integer(
                dplyr::coalesce(as.integer(.data[[htp_col]]), 0L) +
                  dplyr::coalesce(as.integer(.data[[ltp_col]]), 0L)
              )
            )
        } else {
          psp_filtered <- psp_filtered |>
            dplyr::mutate(EvidenceCount = 1L)
        }

        # Keep high-confidence sites (>5 references)
        psp_filtered <- psp_filtered |>
          dplyr::filter(EvidenceCount > 5L) |>
          dplyr::filter(UniProtID %in% uniprot_ids) |>
          dplyr::transmute(
            UniProtID, ResiduePosition, ResidueAA,
            ModificationType, EvidenceCount,
            IsDiseaseSite = FALSE
          )

        ptm_psp <- psp_filtered
        cli::cli_alert_success(
          "PSP: {nrow(ptm_psp)} high-confidence sites (>5 refs) loaded"
        )
      }
    }
  }

  # -- 3c. Merge UniProt + PSP -------------------------------------------------
  ptm_all <- dplyr::bind_rows(ptm_uniprot, ptm_psp) |>
    dplyr::distinct(UniProtID, ResiduePosition, ModificationType, .keep_all = TRUE)

  cli::cli_alert_success("Module 3: {nrow(ptm_all)} PTM sites total")
  ptm_all
}

###############################################################################
# MODULE 4: Protein Sequence Retrieval and 62-mer Tiling
###############################################################################

#' Fetch canonical protein sequences and generate 62-mer tiles centred on
#' modification sites.
#'
#' @param uniprot_ids    Character vector of UniProt accessions.
#' @param ptm_sites      Tibble from Module 3 with PTM site coordinates.
#' @param high_priority  Character vector of UniProt IDs for full-length tiling.
#' @param mod_code_map   Named character vector mapping ModificationType
#'   patterns to single lowercase letters for encoding. Default covers common
#'   modifications.
#' @param cache          Optional cachem cache.
#' @return A tibble with columns: UniProtID, ProteinName, PeptideSequence,
#'   ModType, ModPosition, WindowStart, WindowEnd.
tile_peptides <- function(uniprot_ids,
                          ptm_sites,
                          high_priority = NULL,
                          mod_code_map = NULL,
                          cache = NULL) {

  cli::cli_h2("Module 4: Protein Sequence Retrieval & 62-mer Tiling")

  if (is.null(mod_code_map)) {
    mod_code_map <- c(
      "Citrullin"        = "r",
      "Carbamyl"         = "k",
      "Acetyl"           = "k",
      "Phospho"          = "s",
      "Phosphoserine"    = "s",
      "Phosphothreonine" = "t",
      "Phosphotyrosine"  = "y",
      "Oxidat"           = "m",
      "Malondialdehyde"  = "k",
      "Hydroxyl"         = "p",
      "Methyl"           = "k",
      "Ubiquitin"        = "k",
      "SUMO"             = "k",
      "ADP-ribos"        = "e"
    )
  }

  uniprot_ids <- unique(uniprot_ids[!is.na(uniprot_ids)])

  # -- 4a. Fetch FASTA sequences -----------------------------------------------
  fetch_fastas <- function() {
    cli::cli_alert_info("Fetching FASTA sequences for {length(uniprot_ids)} proteins ...")
    batch_size <- 100L
    batches <- split(uniprot_ids, ceiling(seq_along(uniprot_ids) / batch_size))
    seqs <- list()

    for (i in seq_along(batches)) {
      batch <- batches[[i]]
      query <- paste0(
        "(accession:", paste(batch, collapse = " OR accession:"), ")"
      )

      resp <- tryCatch(
        {
          httr2::request(paste0(.UNIPROT_REST, "/uniprotkb/stream")) |>
            httr2::req_url_query(
              query  = query,
              format = "fasta"
            ) |>
            httr2::req_timeout(120) |>
            httr2::req_retry(max_tries = 3, backoff = ~ 3) |>
            httr2::req_perform()
        },
        error = function(e) {
          cli::cli_warn("FASTA batch {i} failed: {conditionMessage(e)}")
          return(NULL)
        }
      )

      if (is.null(resp)) next

      fasta_text <- httr2::resp_body_string(resp)
      if (nchar(fasta_text) < 5) next

      # Parse FASTA
      tmp_fa <- tempfile(fileext = ".fasta")
      writeLines(fasta_text, tmp_fa)
      aa <- tryCatch(
        Biostrings::readAAStringSet(tmp_fa),
        error = function(e) NULL
      )
      unlink(tmp_fa)

      if (!is.null(aa) && length(aa) > 0) {
        seqs[[i]] <- tibble::tibble(
          Header   = names(aa),
          Sequence = as.character(aa)
        )
      }
      Sys.sleep(0.25)
    }

    dplyr::bind_rows(seqs)
  }

  fasta_tbl <- if (!is.null(cache)) {
    cache_get_or_set(
      cache,
      paste0("fasta_seqs_", digest::digest(sort(uniprot_ids))),
      fetch_fastas
    )
  } else {
    fetch_fastas()
  }

  if (nrow(fasta_tbl) == 0) {
    cli::cli_warn("No FASTA sequences retrieved.")
    return(tibble::tibble(
      UniProtID = character(), ProteinName = character(),
      PeptideSequence = character(), ModType = character(),
      ModPosition = integer(), WindowStart = integer(),
      WindowEnd = integer()
    ))
  }

  # Parse UniProt accession and protein name from header
  fasta_tbl <- fasta_tbl |>
    dplyr::mutate(
      UniProtID   = stringr::str_extract(Header, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
      ProteinName = stringr::str_extract(Header, "(?<=\\s).+?(?=\\sOS=)")
    )

  # -- 4b. Site-centred 62-mers ------------------------------------------------
  cli::cli_alert_info("Generating site-centred 62-mer tiles ...")

  # Helper: extract a 62-mer centred on a position, pad with 'X' at termini
  extract_window <- function(seq_str, centre_pos, tile_len = .TILE_LENGTH) {
    seq_len <- nchar(seq_str)
    half <- tile_len %/% 2L
    start <- centre_pos - half
    end   <- start + tile_len - 1L

    pad_left  <- 0L
    pad_right <- 0L

    if (start < 1L) {
      pad_left <- 1L - start
      start <- 1L
    }
    if (end > seq_len) {
      pad_right <- end - seq_len
      end <- seq_len
    }

    substr_seq <- substring(seq_str, start, end)
    paste0(
      strrep("X", pad_left),
      substr_seq,
      strrep("X", pad_right)
    )
  }

  # Determine mod code for a modification type string
  get_mod_code <- function(mod_type_str, residue_aa) {
    if (is.na(mod_type_str)) return(tolower(residue_aa))
    for (pat in names(mod_code_map)) {
      if (stringr::str_detect(mod_type_str, stringr::regex(pat, ignore_case = TRUE))) {
        return(mod_code_map[[pat]])
      }
    }
    # Default: lowercase of the residue
    if (!is.na(residue_aa) && nchar(residue_aa) == 1L) {
      return(tolower(residue_aa))
    }
    "x"
  }

  # Join PTMs with sequences
  site_tiles <- ptm_sites |>
    dplyr::inner_join(
      fasta_tbl |> dplyr::select(UniProtID, Sequence, ProteinName),
      by = "UniProtID"
    ) |>
    dplyr::filter(!is.na(ResiduePosition), ResiduePosition > 0L)

  if (nrow(site_tiles) > 0) {
    site_tiles <- site_tiles |>
      dplyr::rowwise() |>
      dplyr::mutate(
        WindowStart     = max(1L, ResiduePosition - .TILE_LENGTH %/% 2L),
        WindowEnd       = min(nchar(Sequence),
                              ResiduePosition + .TILE_LENGTH %/% 2L - 1L),
        RawPeptide      = extract_window(Sequence, ResiduePosition),
        # Encode modification
        .mod_code       = get_mod_code(ModificationType, ResidueAA),
        # Position of the modified residue within the 62-mer
        .local_pos      = .TILE_LENGTH %/% 2L + 1L,
        PeptideSequence = paste0(
          substring(RawPeptide, 1, .local_pos - 1L),
          .mod_code,
          substring(RawPeptide, .local_pos + 1L)
        ),
        ModType         = ModificationType
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        UniProtID, ProteinName, PeptideSequence, ModType,
        ModPosition = ResiduePosition, WindowStart, WindowEnd
      )
  }

  # -- 4c. Full-length overlapping tiles for high-priority proteins ------------
  full_tiles <- tibble::tibble()
  if (!is.null(high_priority) && length(high_priority) > 0) {
    hp_seqs <- fasta_tbl |>
      dplyr::filter(UniProtID %in% high_priority)

    if (nrow(hp_seqs) > 0) {
      cli::cli_alert_info(
        "Generating full-length overlapping tiles for {nrow(hp_seqs)} high-priority proteins ..."
      )

      full_tiles <- purrr::pmap_dfr(hp_seqs, function(Header, Sequence,
                                                        UniProtID, ProteinName) {
        seq_len <- nchar(Sequence)
        starts  <- seq(1L, max(1L, seq_len - .TILE_LENGTH + 1L), by = .TILE_STEP)
        # Ensure last tile covers the C-terminus
        if (utils::tail(starts, 1) + .TILE_LENGTH - 1L < seq_len) {
          starts <- c(starts, seq_len - .TILE_LENGTH + 1L)
        }

        purrr::map_dfr(starts, function(s) {
          e <- min(s + .TILE_LENGTH - 1L, seq_len)
          pep <- substring(Sequence, s, e)
          # Pad if at C-terminus and shorter than tile length
          if (nchar(pep) < .TILE_LENGTH) {
            pep <- paste0(pep, strrep("X", .TILE_LENGTH - nchar(pep)))
          }
          tibble::tibble(
            UniProtID       = UniProtID,
            ProteinName     = ProteinName,
            PeptideSequence = pep,
            ModType         = NA_character_,
            ModPosition     = NA_integer_,
            WindowStart     = s,
            WindowEnd       = e
          )
        })
      })
    }
  }

  # -- 4d. Combine and deduplicate ---------------------------------------------
  all_tiles <- dplyr::bind_rows(site_tiles, full_tiles) |>
    dplyr::distinct(UniProtID, PeptideSequence, .keep_all = TRUE)

  cli::cli_alert_success(
    "Module 4: {nrow(all_tiles)} peptide tiles generated"
  )
  all_tiles
}

###############################################################################
# MODULE 5: Prioritisation Scoring
###############################################################################

#' Score and rank candidate 62-mer peptides using a composite score.
#'
#' @param peptides   Tibble from Module 4 with peptide tiles.
#' @param iedb_data  Tibble from Module 2 with IEDB epitope data.
#' @param ctd_data   Tibble from Module 1 with CTD interaction data.
#' @param cache      Optional cachem cache.
#' @param weights    Named numeric vector with component weights. Defaults:
#'   c(iedb = 0.4, ctd = 0.3, hydrophilicity = 0.2, abundance = 0.1).
#' @return The peptide tibble augmented with score columns, sorted descending
#'   by CompositeScore.
score_peptides <- function(peptides,
                           iedb_data,
                           ctd_data,
                           cache = NULL,
                           weights = c(iedb = 0.4, ctd = 0.3,
                                       hydrophilicity = 0.2,
                                       abundance = 0.1)) {

  cli::cli_h2("Module 5: Prioritisation Scoring")

  if (nrow(peptides) == 0) {
    cli::cli_warn("No peptides to score.")
    return(peptides |>
      dplyr::mutate(
        IEDB_Score = numeric(), CTD_Score = numeric(),
        Hydrophilicity_Score = numeric(), Abundance_Score = numeric(),
        CompositeScore = numeric()
      ))
  }

  # -- 5a. IEDB evidence score ------------------------------------------------
  cli::cli_alert_info("Computing IEDB evidence scores ...")

  # Build set of confirmed immunogenic sequences
  iedb_sequences <- iedb_data |>
    dplyr::filter(!is.na(Sequence), nchar(Sequence) >= 5L) |>
    dplyr::pull(Sequence) |>
    unique() |>
    toupper()

  iedb_hc_sequences <- iedb_data |>
    dplyr::filter(HighConfidenceIgG == TRUE, !is.na(Sequence)) |>
    dplyr::pull(Sequence) |>
    unique() |>
    toupper()

  peptides <- peptides |>
    dplyr::mutate(
      .pep_upper = toupper(PeptideSequence),
      IEDB_Score = purrr::map_dbl(.pep_upper, function(pep) {
        # Check if any IEDB epitope is contained within this peptide or vice versa
        exact_match <- pep %in% iedb_sequences
        overlap     <- any(purrr::map_lgl(iedb_sequences, function(ep) {
          stringr::str_detect(pep, stringr::fixed(ep)) ||
            stringr::str_detect(ep, stringr::fixed(
              stringr::str_sub(pep, 10, 53)  # core region
            ))
        }))
        hc_match <- pep %in% iedb_hc_sequences || any(purrr::map_lgl(
          iedb_hc_sequences, function(ep) {
            stringr::str_detect(pep, stringr::fixed(ep))
          }
        ))
        dplyr::case_when(
          hc_match    ~ 1.0,
          exact_match ~ 0.8,
          overlap     ~ 0.5,
          TRUE        ~ 0.0
        )
      })
    )

  # -- 5b. CTD interaction evidence score -------------------------------------
  cli::cli_alert_info("Computing CTD evidence scores ...")

  ctd_ref_counts <- ctd_data |>
    dplyr::filter(!is.na(UniProtID), !is.na(PubMedIDs)) |>
    dplyr::group_by(UniProtID) |>
    dplyr::summarise(
      .n_refs = sum(
        stringr::str_count(PubMedIDs, "\\|") + 1L, na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(.ctd_raw = log10(pmax(.n_refs, 1)))

  # Normalise to 0-1
  max_ctd <- max(ctd_ref_counts$.ctd_raw, na.rm = TRUE)
  if (is.finite(max_ctd) && max_ctd > 0) {
    ctd_ref_counts <- ctd_ref_counts |>
      dplyr::mutate(.ctd_norm = .ctd_raw / max_ctd)
  } else {
    ctd_ref_counts <- ctd_ref_counts |>
      dplyr::mutate(.ctd_norm = 0)
  }

  peptides <- peptides |>
    dplyr::left_join(
      ctd_ref_counts |> dplyr::select(UniProtID, CTD_Score = .ctd_norm),
      by = "UniProtID"
    ) |>
    dplyr::mutate(CTD_Score = dplyr::coalesce(CTD_Score, 0))

  # -- 5c. Surface accessibility (hydrophilicity) proxy -----------------------
  cli::cli_alert_info("Computing hydrophilicity scores ...")

  compute_hydrophilicity <- function(pep_seq, window = 9L) {
    aa_vec <- strsplit(toupper(pep_seq), "")[[1]]
    values <- .KD_SCALE[aa_vec]
    values[is.na(values)] <- 0
    # Invert so that hydrophilic residues score higher
    values <- -values
    if (length(values) < window) {
      return(mean(values, na.rm = TRUE))
    }
    # Sliding window mean, return the max window (most hydrophilic stretch)
    windows <- zoo::rollmean(values, k = window, fill = NA, align = "center")
    max(windows, na.rm = TRUE)
  }

  # Check for zoo; if unavailable use simple mean
  has_zoo <- requireNamespace("zoo", quietly = TRUE)

  peptides <- peptides |>
    dplyr::mutate(
      .hydro_raw = if (has_zoo) {
        purrr::map_dbl(.pep_upper, compute_hydrophilicity)
      } else {
        purrr::map_dbl(.pep_upper, function(pep) {
          aa_vec <- strsplit(pep, "")[[1]]
          vals <- .KD_SCALE[aa_vec]
          vals[is.na(vals)] <- 0
          mean(-vals, na.rm = TRUE)
        })
      }
    )

  # Normalise to 0-1 via min-max
  hydro_range <- range(peptides$.hydro_raw, na.rm = TRUE, finite = TRUE)
  if (diff(hydro_range) > 0) {
    peptides <- peptides |>
      dplyr::mutate(
        Hydrophilicity_Score = (`.hydro_raw` - hydro_range[1]) /
          diff(hydro_range)
      )
  } else {
    peptides <- peptides |>
      dplyr::mutate(Hydrophilicity_Score = 0.5)
  }

  # -- 5d. Protein abundance from Human Protein Atlas -------------------------
  cli::cli_alert_info("Fetching HPA blood protein abundance data ...")

  fetch_hpa_abundance <- function() {
    resp <- tryCatch(
      {
        httr2::request(.HPA_API) |>
          httr2::req_url_query(
            search   = "",
            format   = "json",
            columns  = "g,up,blood_concentration_nM",
            compress = "no"
          ) |>
          httr2::req_timeout(120) |>
          httr2::req_retry(max_tries = 3, backoff = ~ 3) |>
          httr2::req_perform()
      },
      error = function(e) {
        cli::cli_warn("HPA API request failed: {conditionMessage(e)}")
        return(NULL)
      }
    )

    if (is.null(resp)) return(tibble::tibble(UniProtID = character(),
                                              .hpa_nM = numeric()))

    body <- tryCatch(
      httr2::resp_body_json(resp, simplifyVector = TRUE),
      error = function(e) NULL
    )

    if (is.null(body) || length(body) == 0) {
      return(tibble::tibble(UniProtID = character(), .hpa_nM = numeric()))
    }

    tbl <- tryCatch(tibble::as_tibble(body), error = function(e) tibble::tibble())
    if (nrow(tbl) == 0) return(tibble::tibble(UniProtID = character(),
                                               .hpa_nM = numeric()))

    # Normalise column names
    up_col <- intersect(names(tbl), c("Uniprot", "up", "UP", "uniprot"))[1]
    conc_col <- intersect(names(tbl), c("Blood concentration [nM]",
                                         "blood_concentration_nM",
                                         "Blood.concentration..nM."))[1]
    if (is.na(up_col) || is.na(conc_col)) {
      # Try positional
      if (ncol(tbl) >= 3) {
        up_col <- names(tbl)[2]
        conc_col <- names(tbl)[3]
      } else {
        return(tibble::tibble(UniProtID = character(), .hpa_nM = numeric()))
      }
    }

    tbl |>
      dplyr::transmute(
        UniProtID = as.character(.data[[up_col]]),
        .hpa_nM   = as.numeric(.data[[conc_col]])
      ) |>
      dplyr::filter(!is.na(.hpa_nM), .hpa_nM > 0)
  }

  hpa_data <- if (!is.null(cache)) {
    cache_get_or_set(cache, "hpa_blood_abundance", fetch_hpa_abundance)
  } else {
    fetch_hpa_abundance()
  }

  if (nrow(hpa_data) > 0) {
    hpa_data <- hpa_data |>
      dplyr::mutate(.log_nM = log10(.hpa_nM))
    hpa_range <- range(hpa_data$.log_nM, na.rm = TRUE, finite = TRUE)
    if (diff(hpa_range) > 0) {
      hpa_data <- hpa_data |>
        dplyr::mutate(
          Abundance_Score = (.log_nM - hpa_range[1]) / diff(hpa_range)
        ) |>
        dplyr::select(UniProtID, Abundance_Score) |>
        dplyr::distinct(UniProtID, .keep_all = TRUE)
    } else {
      hpa_data <- hpa_data |>
        dplyr::mutate(Abundance_Score = 0.5) |>
        dplyr::select(UniProtID, Abundance_Score) |>
        dplyr::distinct(UniProtID, .keep_all = TRUE)
    }

    peptides <- peptides |>
      dplyr::left_join(hpa_data, by = "UniProtID") |>
      dplyr::mutate(Abundance_Score = dplyr::coalesce(Abundance_Score, 0))
  } else {
    peptides <- peptides |>
      dplyr::mutate(Abundance_Score = 0)
  }

  # -- 5e. Composite score ----------------------------------------------------
  cli::cli_alert_info("Computing composite scores ...")

  peptides <- peptides |>
    dplyr::mutate(
      CompositeScore = weights["iedb"]           * IEDB_Score +
                       weights["ctd"]            * CTD_Score +
                       weights["hydrophilicity"] * Hydrophilicity_Score +
                       weights["abundance"]      * Abundance_Score
    ) |>
    dplyr::select(
      -dplyr::starts_with("."),  # drop internal helper columns
    ) |>
    dplyr::arrange(dplyr::desc(CompositeScore))

  cli::cli_alert_success(
    "Module 5: scoring complete. Top score = {round(max(peptides$CompositeScore, na.rm = TRUE), 4)}"
  )

  peptides
}

###############################################################################
# MODULE 6: Master Pipeline Function
###############################################################################

#' Run the full ChemScan pipeline end-to-end.
#'
#' @param chemicals  Character vector of chemical names or CAS numbers to focus
#'   on. NULL uses all CTD interactions.
#' @param proteins   Character vector of UniProt IDs to focus on. NULL uses all
#'   proteins identified by CTD.
#' @param ptm_file   Optional path to a PhosphoSitePlus bulk download TSV file.
#' @param n_peptides Maximum number of peptides to include in the final library.
#'   Default 50 000.
#' @param output_dir Path to write intermediate and final outputs. Will be
#'   created if it does not exist.
#' @return A tibble: the final prioritised peptide library (also written to
#'   output_dir/chemscan_library.tsv).
chemscan_pipeline <- function(chemicals  = NULL,
                              proteins   = NULL,
                              ptm_file   = NULL,
                              n_peptides = 50000L,
                              output_dir = "chemscan_output") {

  cli::cli_h1("ChemScan Pipeline")

  # ── Setup ------------------------------------------------------------------
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir <- file.path(output_dir, "cache")
  cache <- init_cache(cache_dir)

  # ── Module 1: CTD ----------------------------------------------------------
  cli::cli_rule()
  ctd_data <- fetch_ctd_interactions(chemicals = chemicals, cache = cache)
  readr::write_tsv(
    ctd_data, file.path(output_dir, "01_ctd_interactions.tsv")
  )

  # ── Module 2: IEDB ---------------------------------------------------------
  cli::cli_rule()
  iedb_data <- fetch_iedb_modified_epitopes(cache = cache)
  readr::write_tsv(
    iedb_data, file.path(output_dir, "02_iedb_epitopes.tsv")
  )

  # Collect all UniProt IDs of interest
  all_uniprot <- unique(c(
    ctd_data$UniProtID,
    iedb_data$UniProtID,
    proteins
  ))
  all_uniprot <- all_uniprot[!is.na(all_uniprot) & nchar(all_uniprot) > 0]

  if (!is.null(proteins)) {
    all_uniprot <- unique(c(all_uniprot, proteins))
  }

  cli::cli_alert_info("Total unique UniProt IDs: {length(all_uniprot)}")

  # ── Module 3: PTM sites ----------------------------------------------------
  cli::cli_rule()
  ptm_sites <- fetch_ptm_sites(
    uniprot_ids = all_uniprot,
    ptm_file    = ptm_file,
    cache       = cache
  )
  readr::write_tsv(
    ptm_sites, file.path(output_dir, "03_ptm_sites.tsv")
  )

  # ── Module 4: 62-mer tiling ------------------------------------------------
  cli::cli_rule()

  # High-priority: proteins with IEDB evidence or high CTD support
  high_priority <- unique(c(
    iedb_data |>
      dplyr::filter(HighConfidenceIgG == TRUE) |>
      dplyr::pull(UniProtID),
    ctd_data |>
      dplyr::filter(
        !is.na(PubMedIDs),
        stringr::str_count(PubMedIDs, "\\|") >= 4L
      ) |>
      dplyr::pull(UniProtID)
  ))
  high_priority <- high_priority[!is.na(high_priority)]

  peptides <- tile_peptides(
    uniprot_ids   = all_uniprot,
    ptm_sites     = ptm_sites,
    high_priority = high_priority,
    cache         = cache
  )
  readr::write_tsv(
    peptides, file.path(output_dir, "04_peptide_tiles.tsv")
  )

  # ── Module 5: Scoring and ranking ------------------------------------------
  cli::cli_rule()
  scored <- score_peptides(
    peptides  = peptides,
    iedb_data = iedb_data,
    ctd_data  = ctd_data,
    cache     = cache
  )

  # ── Trim to library size ---------------------------------------------------
  if (nrow(scored) > n_peptides) {
    cli::cli_alert_info(
      "Trimming library from {nrow(scored)} to {n_peptides} peptides"
    )
    scored <- scored |> dplyr::slice_head(n = n_peptides)
  }

  # ── Write final output -----------------------------------------------------
  output_path <- file.path(output_dir, "chemscan_library.tsv")
  readr::write_tsv(scored, output_path)

  cli::cli_rule()
  cli::cli_alert_success(
    "ChemScan pipeline complete. Final library: {nrow(scored)} peptides"
  )
  cli::cli_alert_info("Output written to {.file {output_path}}")

  scored
}

###############################################################################
# WORKED EXAMPLE (commented)
###############################################################################

# library(chemscan)  # or source("chemscan_pipeline.R")
#
# # Run the ChemScan pipeline for three test chemicals
# results <- chemscan_pipeline(
#   chemicals  = c("acrylamide", "benzene", "cisplatin"),
#   proteins   = NULL,
#   ptm_file   = NULL,             # or "path/to/Phosphorylation_site_dataset.gz"
#   n_peptides = 50000,
#   output_dir = "chemscan_output"
# )
#
# # Inspect the top-ranked candidates
# results |>
#   dplyr::select(UniProtID, ProteinName, ModType, ModPosition,
#                 CompositeScore, IEDB_Score, CTD_Score) |>
#   head(20) |>
#   print()
#
# # Filter for high-confidence epitopes only
# high_confidence <- results |>
#   dplyr::filter(IEDB_Score >= 0.8, CTD_Score >= 0.5)
#
# cat(sprintf("\nHigh-confidence peptides: %d\n", nrow(high_confidence)))
