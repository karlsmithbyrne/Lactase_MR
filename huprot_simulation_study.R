#!/usr/bin/env Rscript
## ============================================================================
##
##  HuProt Autoantibody Array Simulation Study
##  ---------------------------------------------------------------------------
##  Evaluates Type I error and statistical power of single-protein tests
##  (Fisher Exact, Wilcoxon) versus set-based aggregation tests (Burden,
##  SKAT, SKAT-O, ACAT-V, ACAT-O) on simulated HuProt-like data.
##
##  Author:  Claude (Anthropic) for Karl Smith Byrne
##  Date:    2026-02-28
##  License: MIT
##
## ============================================================================


## ---------- 0. DEPENDENCIES & GLOBAL PARAMETERS ----------

suppressPackageStartupMessages({
  library(tidyverse)
  library(SKAT)          # SKAT / Burden / SKAT-O set-based tests
  library(Matrix)        # sparse matrix utilities
})

## Reproducibility
set.seed(42)

## -- Simulation geometry --
N_CASES       <- 20L       # number of cases per iteration
N_CONTROLS    <- 20L       # number of controls per iteration
N_PROTEINS    <- 5000L     # total proteins on the array
N_SETS        <- 50L       # number of protein sets ("genes" / "pathways")
SET_SIZE_RANGE <- c(2L, 20L) # min/max proteins per set
N_ITER        <- 100L      # Monte-Carlo iterations

## -- Signal injection --
N_SIGNAL_SETS   <- 5L      # how many of the 50 sets receive a true signal
BURDEN_SHIFT    <- 1.2     # mean log2-MFI shift for burden-style signal
SKAT_SHIFT_SD   <- 1.0     # SD of per-protein shifts for SKAT-style signal
HIT_THRESHOLD_Z <- 2.0     # Z-score above which a protein is called a "hit"
BACKGROUND_MEAN <- 6.0     # background log2 MFI mean (typical for HuProt)
BACKGROUND_SD   <- 1.5     # background log2 MFI SD

## -- Testing --
ALPHA <- 0.05              # significance threshold


## ============================================================================
## ---------- 1. DATA GENERATION ----------
## ============================================================================
##
##  simulate_huprot()
##  -----------------
##  Returns a list containing:
##    $mfi       : N × P matrix of continuous log2-MFI values
##    $hits      : N × P binary hit matrix (0/1)
##    $group     : length-N factor ("case" / "control")
##    $sets      : list of length N_SETS; each element is an integer vector
##                 of column indices into mfi / hits
##    $signal_ix : integer vector of set indices that received injected signal
##    $signal_type : character vector ("burden" or "skat") for each signal set
##

simulate_huprot <- function(n_cases       = N_CASES,
                            n_controls    = N_CONTROLS,
                            n_proteins    = N_PROTEINS,
                            n_sets        = N_SETS,
                            set_size_range = SET_SIZE_RANGE,
                            n_signal_sets = N_SIGNAL_SETS,
                            burden_shift  = BURDEN_SHIFT,
                            skat_shift_sd = SKAT_SHIFT_SD,
                            hit_z         = HIT_THRESHOLD_Z,
                            bg_mean       = BACKGROUND_MEAN,
                            bg_sd         = BACKGROUND_SD) {

  n_total <- n_cases + n_controls
  group   <- factor(c(rep("case", n_cases), rep("control", n_controls)),
                    levels = c("control", "case"))

  ## ------------------------------------------------------------------
  ## 1a. Background continuous MFI data (samples × proteins)
  ## ------------------------------------------------------------------
  ## Each protein gets its own mean drawn from the global background
  ## distribution, introducing realistic inter-protein variance.
  protein_means <- rnorm(n_proteins, mean = bg_mean, sd = bg_sd * 0.3)
  protein_sds   <- runif(n_proteins, min = bg_sd * 0.6, max = bg_sd * 1.4)

  mfi <- matrix(NA_real_, nrow = n_total, ncol = n_proteins)
  for (j in seq_len(n_proteins)) {
    mfi[, j] <- rnorm(n_total, mean = protein_means[j], sd = protein_sds[j])
  }

  ## ------------------------------------------------------------------
  ## 1b. Assign proteins to sets (non-overlapping for simplicity)
  ## ------------------------------------------------------------------
  set_sizes <- sample(set_size_range[1]:set_size_range[2],
                      size = n_sets, replace = TRUE)

  ## We only need sum(set_sizes) proteins assigned; the rest are "orphans"
  ## not belonging to any tested set (they still contribute to the array).
  total_assigned <- sum(set_sizes)
  if (total_assigned > n_proteins) {
    stop("Total assigned proteins (", total_assigned,
         ") exceeds n_proteins (", n_proteins, ").")
  }

  assigned_cols <- sample.int(n_proteins, size = total_assigned)
  sets <- vector("list", n_sets)
  idx  <- 1L
  for (s in seq_len(n_sets)) {
    sets[[s]] <- assigned_cols[idx:(idx + set_sizes[s] - 1L)]
    idx <- idx + set_sizes[s]
  }
  names(sets) <- paste0("set_", seq_len(n_sets))

  ## ------------------------------------------------------------------
  ## 1c. Inject true signal into a subset of pathways (cases only)
  ## ------------------------------------------------------------------
  signal_ix   <- sort(sample.int(n_sets, size = n_signal_sets))
  signal_type <- character(n_signal_sets)

  ## Half burden-style, half SKAT-style (rounded)
  n_burden <- ceiling(n_signal_sets / 2)
  signal_type[seq_len(n_burden)] <- "burden"
  if (n_signal_sets > n_burden) {
    signal_type[(n_burden + 1):n_signal_sets] <- "skat"
  }

  case_rows <- which(group == "case")

  for (k in seq_along(signal_ix)) {
    s_idx <- signal_ix[k]
    cols  <- sets[[s_idx]]

    if (signal_type[k] == "burden") {
      ## All proteins in the set shift upward by a fixed amount in cases
      for (j in cols) {
        mfi[case_rows, j] <- mfi[case_rows, j] + burden_shift
      }
    } else {
      ## Heterogeneous: each protein gets a random shift (pos or neg)
      per_protein_shift <- rnorm(length(cols), mean = 0, sd = skat_shift_sd)
      for (idx_j in seq_along(cols)) {
        mfi[case_rows, cols[idx_j]] <- mfi[case_rows, cols[idx_j]] +
          per_protein_shift[idx_j]
      }
    }
  }

  ## ------------------------------------------------------------------
  ## 1d. Derive binary hit matrix from Z-scored MFI
  ## ------------------------------------------------------------------
  ## Z-score each protein across all samples, then threshold.
  col_means <- colMeans(mfi)
  col_sds   <- apply(mfi, 2, sd)
  ## Guard against zero-variance columns (extremely unlikely but defensive)
  col_sds[col_sds < 1e-12] <- 1e-12

  z_mat <- sweep(mfi, 2, col_means, `-`)
  z_mat <- sweep(z_mat, 2, col_sds, `/`)

  hits <- ifelse(z_mat >= hit_z, 1L, 0L)

  ## ------------------------------------------------------------------
  ## Return the full simulation object
  ## ------------------------------------------------------------------
  list(
    mfi         = mfi,
    hits        = hits,
    group       = group,
    sets        = sets,
    signal_ix   = signal_ix,
    signal_type = signal_type
  )
}


## ============================================================================
## ---------- 2. ANALYTICAL WRAPPERS ----------
## ============================================================================


## ---------- 2a. Single-Protein Fisher Exact Test (binary hits) ----------
##
## Returns a numeric vector of length P with one p-value per protein.
## For each protein, builds a 2×2 table: (hit / no-hit) × (case / control).

fisher_per_protein <- function(hits, group) {
  n_prot   <- ncol(hits)
  p_values <- numeric(n_prot)

  for (j in seq_len(n_prot)) {
    tbl <- table(factor(hits[, j], levels = c(0, 1)),
                 factor(group,     levels = c("control", "case")))
    ## fisher.test can fail on degenerate tables; wrap defensively
    p_values[j] <- tryCatch(
      fisher.test(tbl)$p.value,
      error = function(e) NA_real_
    )
  }
  p_values
}


## ---------- 2b. Single-Protein Wilcoxon Rank-Sum Test (continuous MFI) ---
##
## Returns a numeric vector of length P with one p-value per protein.

wilcox_per_protein <- function(mfi, group) {
  n_prot   <- ncol(mfi)
  p_values <- numeric(n_prot)

  for (j in seq_len(n_prot)) {
    p_values[j] <- tryCatch(
      wilcox.test(mfi[, j] ~ group)$p.value,
      error = function(e) NA_real_
    )
  }
  p_values
}


## ---------- 2c. ACAT (Aggregated Cauchy Association Test) ----------------
##
## Liu & Xie (2020, JASA). Combines arbitrary p-values into a single
## test statistic using the Cauchy distribution.
##
##   T_ACAT = sum(w_i * tan((0.5 - p_i) * pi)) / sum(w_i)
##
## The p-value of T_ACAT is:  0.5 - arctan(T_ACAT) / pi
##
## This implementation uses equal weights by default.

acat_combine <- function(p_values, weights = NULL) {

  ## Remove NAs
  valid <- !is.na(p_values)
  p     <- p_values[valid]

  if (length(p) == 0L) return(NA_real_)

  ## Clamp p-values away from exact 0 or 1 to avoid ±Inf from tan()
  p <- pmax(p, .Machine$double.eps)
  p <- pmin(p, 1 - .Machine$double.eps)

  if (is.null(weights)) {
    weights <- rep(1, length(p))
  } else {
    weights <- weights[valid]
  }

  ## Cauchy transformation
  t_stat <- sum(weights * tan((0.5 - p) * pi)) / sum(weights)

  ## P-value from standard Cauchy survival function
  acat_p <- 0.5 - atan(t_stat) / pi

  ## Clamp to valid range
  acat_p <- max(acat_p, .Machine$double.eps)
  acat_p <- min(acat_p, 1)

  acat_p
}


## ---------- 2d. Set-Level SKAT / Burden / SKAT-O (binary hit matrix) -----
##
## We repurpose the SKAT package by treating the binary hit matrix as
## a "genotype" matrix (rare-variant-like coding: 0 = no hit, 1 = hit).
## The phenotype is case/control status.
##
## Returns a tibble with one row per set and columns:
##   set_name, skat_p, burden_p, skato_p

skat_burden_per_set <- function(hits, group, sets) {

  ## Phenotype vector: 0 = control, 1 = case (numeric for SKAT)
  y <- as.integer(group == "case")

  ## Fit the null model ONCE (no covariates, dichotomous outcome)
  null_model <- tryCatch(
    SKAT_Null_Model(y ~ 1, out_type = "D"),
    error = function(e) {
      warning("SKAT_Null_Model failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(null_model)) {
    ## Return NA for everything if null model fitting fails
    return(tibble(
      set_name  = names(sets),
      skat_p    = NA_real_,
      burden_p  = NA_real_,
      skato_p   = NA_real_
    ))
  }

  n_sets  <- length(sets)
  skat_p  <- numeric(n_sets)
  burden_p <- numeric(n_sets)
  skato_p <- numeric(n_sets)

  for (s in seq_len(n_sets)) {
    cols <- sets[[s]]
    Z    <- hits[, cols, drop = FALSE]

    ## Defensive: SKAT requires at least some variation in the genotype
    ## matrix. If every column is constant, skip.
    col_vars <- apply(Z, 2, var)
    Z <- Z[, col_vars > 0, drop = FALSE]

    if (ncol(Z) < 1L) {
      skat_p[s]   <- NA_real_
      burden_p[s] <- NA_real_
      skato_p[s]  <- NA_real_
      next
    }

    ## -- SKAT --
    skat_p[s] <- tryCatch({
      res <- SKAT(Z, null_model, kernel = "linear.weighted",
                  method = "davies")
      res$p.value
    }, error = function(e) NA_real_)

    ## -- Burden --
    burden_p[s] <- tryCatch({
      res <- SKAT(Z, null_model, kernel = "linear.weighted",
                  method = "davies", r.corr = 1)
      res$p.value
    }, error = function(e) NA_real_)

    ## -- SKAT-O (optimal rho) --
    skato_p[s] <- tryCatch({
      res <- SKAT(Z, null_model, kernel = "linear.weighted",
                  method = "optimal.adj")
      res$p.value
    }, error = function(e) NA_real_)
  }

  tibble(
    set_name = names(sets),
    skat_p   = skat_p,
    burden_p = burden_p,
    skato_p  = skato_p
  )
}


## ---------- 2e. ACAT-V: Combine per-protein Fisher p-values per set ------
##
## For each set, gather the Fisher Exact p-values from member proteins and
## aggregate them via the Cauchy combination (ACAT).

acatv_per_set <- function(fisher_pvals, sets) {
  n_sets <- length(sets)
  acatv_p <- numeric(n_sets)

  for (s in seq_len(n_sets)) {
    cols <- sets[[s]]
    acatv_p[s] <- acat_combine(fisher_pvals[cols])
  }

  tibble(
    set_name = names(sets),
    acatv_p  = acatv_p
  )
}


## ---------- 2f. ACAT-O: Omnibus combination of set-level tests -----------
##
## For each set, combine the SKAT, Burden, and ACAT-V p-values into a
## single omnibus p-value using ACAT.

acato_per_set <- function(skat_results, acatv_results) {

  ## Join on set_name
  combined <- skat_results %>%
    left_join(acatv_results, by = "set_name")

  acato_p <- numeric(nrow(combined))

  for (s in seq_len(nrow(combined))) {
    p_vec <- c(
      combined$skat_p[s],
      combined$burden_p[s],
      combined$acatv_p[s]
    )
    acato_p[s] <- acat_combine(p_vec)
  }

  combined %>%
    mutate(acato_p = acato_p)
}


## ============================================================================
## ---------- 3. SINGLE-ITERATION RUNNER ----------
## ============================================================================
##
##  run_one_iteration()
##  -------------------
##  Simulates data, runs all tests, and returns a tidy tibble of set-level
##  results with columns indicating whether each set was a true signal set.

run_one_iteration <- function(iter_id = 1L) {

  ## -- 3a. Simulate data --
  sim <- simulate_huprot()

  ## -- 3b. Single-protein tests (full array) --
  fisher_pvals  <- fisher_per_protein(sim$hits, sim$group)
  wilcox_pvals  <- wilcox_per_protein(sim$mfi,  sim$group)

  ## -- 3c. Set-level tests --
  skat_res  <- skat_burden_per_set(sim$hits, sim$group, sim$sets)
  acatv_res <- acatv_per_set(fisher_pvals, sim$sets)
  full_res  <- acato_per_set(skat_res, acatv_res)

  ## -- 3d. Aggregate single-protein tests to set level --
  ##        For Fisher and Wilcoxon, we compute the minimum p-value within
  ##        each set (Bonferroni-style) AND the ACAT combination.
  n_sets <- length(sim$sets)
  fisher_min_p  <- numeric(n_sets)
  wilcox_min_p  <- numeric(n_sets)
  fisher_acat_p <- numeric(n_sets)
  wilcox_acat_p <- numeric(n_sets)

  for (s in seq_len(n_sets)) {
    cols <- sim$sets[[s]]
    set_size <- length(cols)

    ## Minimum p-value with Bonferroni correction by set size
    fp <- fisher_pvals[cols]
    wp <- wilcox_pvals[cols]

    fisher_min_p[s]  <- min(fp, na.rm = TRUE) * set_size
    wilcox_min_p[s]  <- min(wp, na.rm = TRUE) * set_size

    ## Clamp Bonferroni-corrected p-values to [0, 1]
    fisher_min_p[s]  <- min(fisher_min_p[s], 1)
    wilcox_min_p[s]  <- min(wilcox_min_p[s], 1)

    ## ACAT combination of single-protein p-values
    fisher_acat_p[s] <- acat_combine(fp)
    wilcox_acat_p[s] <- acat_combine(wp)
  }

  ## -- 3e. Assemble the results tibble --
  full_res %>%
    mutate(
      fisher_bonf_p  = fisher_min_p,
      wilcox_bonf_p  = wilcox_min_p,
      fisher_acat_p  = fisher_acat_p,
      wilcox_acat_p  = wilcox_acat_p,
      is_signal      = seq_len(n_sets) %in% sim$signal_ix,
      signal_type    = case_when(
        seq_len(n_sets) %in% sim$signal_ix[sim$signal_type == "burden"] ~ "burden",
        seq_len(n_sets) %in% sim$signal_ix[sim$signal_type == "skat"]   ~ "skat",
        TRUE ~ "null"
      ),
      iteration = iter_id
    )
}


## ============================================================================
## ---------- 4. EXECUTION LOOP ----------
## ============================================================================

message("=== HuProt Simulation Study ===")
message("Iterations: ", N_ITER,
        " | Samples: ", N_CASES, " cases + ", N_CONTROLS, " controls",
        " | Proteins: ", N_PROTEINS,
        " | Sets: ", N_SETS,
        " (", N_SIGNAL_SETS, " with signal)")
message("Starting simulation at ", Sys.time(), "\n")

all_results <- vector("list", N_ITER)

for (i in seq_len(N_ITER)) {
  if (i %% 10 == 0 || i == 1) {
    message("  Iteration ", i, " / ", N_ITER, " [", Sys.time(), "]")
  }
  all_results[[i]] <- tryCatch(
    run_one_iteration(iter_id = i),
    error = function(e) {
      warning("Iteration ", i, " failed: ", conditionMessage(e))
      NULL
    }
  )
}

## Bind all iterations into a single tibble, dropping any that failed
results <- bind_rows(all_results)
message("\nSimulation complete at ", Sys.time())
message("Total set-level results: ", nrow(results), " rows\n")


## ============================================================================
## ---------- 5. EVALUATION: TYPE I ERROR & POWER ----------
## ============================================================================

## The test columns we want to evaluate
test_cols <- c("skat_p", "burden_p", "skato_p",
               "acatv_p", "acato_p",
               "fisher_bonf_p", "wilcox_bonf_p",
               "fisher_acat_p", "wilcox_acat_p")

## -- 5a. Type I Error Rate --
## Proportion of null (non-signal) sets that are called significant
## at alpha = 0.05 across all iterations.

null_results <- results %>% filter(!is_signal)

type1_error <- null_results %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "test", values_to = "type1_error") %>%
  arrange(test)

message("== TYPE I ERROR RATE (alpha = ", ALPHA, ") ==")
message("   Expected under null: ", ALPHA)
message(strrep("-", 50))
for (r in seq_len(nrow(type1_error))) {
  message(sprintf("   %-20s : %.4f", type1_error$test[r], type1_error$type1_error[r]))
}
message("")

## -- 5b. Overall Power --
## Proportion of true signal sets that are called significant.

signal_results <- results %>% filter(is_signal)

power_overall <- signal_results %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "test", values_to = "power") %>%
  arrange(test)

message("== STATISTICAL POWER (alpha = ", ALPHA, ") ==")
message(strrep("-", 50))
for (r in seq_len(nrow(power_overall))) {
  message(sprintf("   %-20s : %.4f", power_overall$test[r], power_overall$power[r]))
}
message("")

## -- 5c. Power stratified by signal type (burden vs. SKAT) --

power_by_type <- signal_results %>%
  group_by(signal_type) %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(-signal_type, names_to = "test", values_to = "power") %>%
  pivot_wider(names_from = signal_type, values_from = power) %>%
  arrange(test)

message("== POWER BY SIGNAL TYPE ==")
message(sprintf("   %-20s   %-10s   %-10s", "Test", "Burden", "SKAT"))
message(strrep("-", 50))
for (r in seq_len(nrow(power_by_type))) {
  message(sprintf("   %-20s : %.4f     %.4f",
                  power_by_type$test[r],
                  power_by_type$burden[r],
                  power_by_type$skat[r]))
}
message("")


## ============================================================================
## ---------- 6. DIAGNOSTIC PLOTS (Optional, saved to disk) ----------
## ============================================================================

## Only attempt plotting if we are in an interactive session or if the user
## wants output files. Wrap in tryCatch so a missing display doesn't crash.

tryCatch({

  ## 6a. Type I Error bar chart
  p1 <- ggplot(type1_error, aes(x = reorder(test, type1_error), y = type1_error)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_hline(yintercept = ALPHA, linetype = "dashed", colour = "red") +
    coord_flip() +
    labs(title = "Type I Error Rate by Test Method",
         subtitle = paste0(N_ITER, " iterations | ", N_CASES, " vs ", N_CONTROLS,
                           " | ", N_PROTEINS, " proteins | alpha = ", ALPHA),
         x = NULL, y = "Type I Error Rate") +
    theme_minimal(base_size = 12)

  ggsave("type1_error_barplot.pdf", p1, width = 8, height = 5)
  message("Saved: type1_error_barplot.pdf")

  ## 6b. Power comparison (grouped by signal type)
  p2 <- power_by_type %>%
    pivot_longer(c(burden, skat), names_to = "signal_type", values_to = "power") %>%
    ggplot(aes(x = reorder(test, power), y = power, fill = signal_type)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c(burden = "darkorange", skat = "purple"),
                      labels = c("Burden (directional)", "SKAT (heterogeneous)")) +
    labs(title = "Statistical Power by Test Method and Signal Architecture",
         subtitle = paste0(N_ITER, " iterations | alpha = ", ALPHA),
         x = NULL, y = "Power", fill = "Signal Type") +
    theme_minimal(base_size = 12)

  ggsave("power_comparison_barplot.pdf", p2, width = 9, height = 5)
  message("Saved: power_comparison_barplot.pdf")

  ## 6c. QQ-plot of null p-values for each test (check calibration)
  qq_data <- null_results %>%
    select(set_name, iteration, all_of(test_cols)) %>%
    pivot_longer(all_of(test_cols), names_to = "test", values_to = "pvalue") %>%
    filter(!is.na(pvalue)) %>%
    group_by(test) %>%
    arrange(pvalue) %>%
    mutate(expected = ppoints(n())) %>%
    ungroup()

  p3 <- ggplot(qq_data, aes(x = -log10(expected), y = -log10(pvalue))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.3, size = 0.6, colour = "steelblue") +
    facet_wrap(~ test, scales = "free") +
    labs(title = "QQ Plots of Null P-values (calibration check)",
         x = expression(-log[10](expected)),
         y = expression(-log[10](observed))) +
    theme_minimal(base_size = 10)

  ggsave("null_qq_plots.pdf", p3, width = 12, height = 8)
  message("Saved: null_qq_plots.pdf")

}, error = function(e) {
  message("Note: Plotting skipped (", conditionMessage(e), ")")
})


## ============================================================================
## ---------- 7. SAVE RESULTS ----------
## ============================================================================

write_csv(results, "simulation_results_full.csv")
message("Saved: simulation_results_full.csv")

write_csv(type1_error, "type1_error_summary.csv")
message("Saved: type1_error_summary.csv")

write_csv(power_overall, "power_overall_summary.csv")
message("Saved: power_overall_summary.csv")

write_csv(power_by_type, "power_by_type_summary.csv")
message("Saved: power_by_type_summary.csv")

message("\n=== Done. ===")
