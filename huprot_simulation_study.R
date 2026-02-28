#!/usr/bin/env Rscript
## ============================================================================
##
##  HuProt Autoantibody Array Simulation Study (v2)
##  ---------------------------------------------------------------------------
##  Evaluates Type I error and statistical power of single-protein tests
##  (limma, Firth logistic) and set-based aggregation tests (Burden, SKAT,
##  SKAT-O, Firth-Burden, ACAT-V) across a grid of sample sizes on simulated
##  HuProt-like data.
##
##  Methods:
##    Single-Protein:
##      - limma (Empirical Bayes moderated t-test) on continuous log2-MFI
##      - Firth penalized logistic regression (logistf) on binary hits
##
##    Set-Based:
##      - Continuous SKAT / Burden (SKAT package on continuous Z-scores)
##      - Binary SKAT / Burden (SKAT package on binary hit matrix)
##      - Firth-Burden (sum-score of binary hits per set → Firth logistic)
##
##    P-value Combination:
##      - ACAT-V (Cauchy combination) of per-protein limma p-values
##      - ACAT-V (Cauchy combination) of per-protein logistf p-values
##
##  Signal Architectures:
##      - Burden: Directional, cumulative shifts across the set
##      - SKAT: Heterogeneous shifts (some proteins up, some down)
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
  library(limma)         # Empirical Bayes moderated t-tests
  library(logistf)       # Firth penalized logistic regression
})

## Reproducibility — use L'Ecuyer-CMRG for parallel-safe streams
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

## -- Simulation geometry --
N_GRID        <- c(20L, 50L, 100L, 500L)  # cases = controls at each grid point
N_PROTEINS    <- 5000L     # total proteins on the array
N_SETS        <- 50L       # number of protein sets ("pathways")
SET_SIZE_RANGE <- c(2L, 20L) # min/max proteins per set
N_ITER        <- 500L      # Monte-Carlo iterations per grid point

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
## ---------- 1. DATA GENERATING MECHANISM (DGM) ----------
## ============================================================================
##
##  simulate_huprot()
##  -----------------
##  Returns a list containing:
##    $mfi         : N × P matrix of continuous log2-MFI values
##    $z_mat       : N × P matrix of Z-scored MFI (protein-wise standardised)
##    $hits        : N × P sparse binary hit matrix (0/1)
##    $group       : length-N factor ("case" / "control")
##    $y           : length-N integer (0 = control, 1 = case)
##    $sets        : list of length N_SETS; each element is an integer vector
##                   of column indices into mfi / hits
##    $signal_ix   : integer vector of set indices that received injected signal
##    $signal_type : character vector ("burden" or "skat") for each signal set
##

simulate_huprot <- function(n_cases       = 20L,
                            n_controls    = 20L,
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
  y       <- as.integer(group == "case")

  ## ------------------------------------------------------------------
  ## 1a. Background continuous MFI data (samples × proteins)
  ## ------------------------------------------------------------------
  ## Each protein gets its own mean drawn from the global background
  ## distribution, introducing realistic inter-protein variance.
  protein_means <- rnorm(n_proteins, mean = bg_mean, sd = bg_sd * 0.3)
  protein_sds   <- runif(n_proteins, min = bg_sd * 0.6, max = bg_sd * 1.4)

  ## Vectorised generation: draw from N(protein_mean_j, protein_sd_j) for each column
  mfi <- vapply(seq_len(n_proteins), function(j) {
    rnorm(n_total, mean = protein_means[j], sd = protein_sds[j])
  }, numeric(n_total))

  ## ------------------------------------------------------------------
  ## 1b. Assign proteins to sets (non-overlapping for simplicity)
  ## ------------------------------------------------------------------
  set_sizes <- sample(set_size_range[1]:set_size_range[2],
                      size = n_sets, replace = TRUE)

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
      mfi[case_rows, cols] <- mfi[case_rows, cols] + burden_shift
    } else {
      ## Heterogeneous: each protein gets a random shift (pos or neg)
      per_protein_shift <- rnorm(length(cols), mean = 0, sd = skat_shift_sd)
      ## Broadcast shift across case rows
      mfi[case_rows, cols] <- sweep(
        mfi[case_rows, cols, drop = FALSE], 2,
        per_protein_shift, `+`
      )
    }
  }

  ## ------------------------------------------------------------------
  ## 1d. Derive Z-scored matrix and binary hit matrix
  ## ------------------------------------------------------------------
  col_means <- colMeans(mfi)
  col_sds   <- apply(mfi, 2, sd)
  ## Guard against zero-variance columns
  col_sds[col_sds < 1e-12] <- 1e-12

  z_mat <- sweep(mfi, 2, col_means, `-`)
  z_mat <- sweep(z_mat, 2, col_sds, `/`)

  hits <- ifelse(z_mat >= hit_z, 1L, 0L)

  ## ------------------------------------------------------------------
  ## Return the full simulation object
  ## ------------------------------------------------------------------
  list(
    mfi         = mfi,
    z_mat       = z_mat,
    hits        = hits,
    group       = group,
    y           = y,
    sets        = sets,
    signal_ix   = signal_ix,
    signal_type = signal_type
  )
}


## ============================================================================
## ---------- 2. ANALYTICAL WRAPPERS ----------
## ============================================================================


## ---------- 2a. limma: Empirical Bayes moderated t-test (continuous) --------
##
## Fits a simple two-group (case vs control) design using limma's eBayes.
## Returns a numeric vector of length P with one p-value per protein.

limma_per_protein <- function(mfi, group) {
  design <- model.matrix(~ group)  # intercept + groupcase
  fit    <- lmFit(t(mfi), design)  # limma expects features × samples

  fit    <- eBayes(fit)
  ## Extract p-value for the "groupcase" coefficient (column 2)
  topTable(fit, coef = 2, number = ncol(mfi), sort.by = "none")$P.Value
}


## ---------- 2b. Firth logistic regression per protein (binary hits) ---------
##
## Uses logistf to fit a penalised logistic model:
##   hit_j ~ case/control
## for each protein j. This handles complete/quasi-complete separation that
## plagues standard glm() with sparse binary data.
##
## Returns a numeric vector of length P with one p-value per protein.

firth_per_protein <- function(hits, y) {
  n_prot   <- ncol(hits)
  p_values <- numeric(n_prot)

  for (j in seq_len(n_prot)) {
    h <- hits[, j]

    ## Defensive: if the column is constant (all 0 or all 1), skip
    if (var(h) < 1e-12) {
      p_values[j] <- NA_real_
      next
    }

    p_values[j] <- tryCatch({
      fit <- logistf(h ~ y, firth = TRUE, plconf = NULL)
      ## The coefficient for y is the second element; extract its p-value
      fit$prob[2]
    }, error = function(e) NA_real_)
  }
  p_values
}


## ---------- 2c. ACAT (Aggregated Cauchy Association Test) -------------------
##
## Liu & Xie (2020, JASA). Combines arbitrary p-values into a single
## test statistic using the Cauchy distribution.
##
##   T_ACAT = sum(w_i * tan((0.5 - p_i) * pi)) / sum(w_i)
##
## The p-value of T_ACAT is:  0.5 - arctan(T_ACAT) / pi
##
## Uses equal weights by default.

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
  max(min(acat_p, 1), .Machine$double.eps)
}


## ---------- 2d. Set-Level SKAT / Burden on CONTINUOUS Z-scores --------------
##
## Applies the SKAT package to the continuous Z-score matrix. The outcome
## is dichotomous (case/control). We treat each column of the Z-matrix as
## a continuous "variant" — analogous to rare-variant burden/SKAT tests
## but operating on quantitative protein measurements.
##
## Returns a tibble with one row per set and columns:
##   set_name, cont_skat_p, cont_burden_p

skat_burden_continuous <- function(z_mat, y, sets) {

  ## Fit the null model ONCE (no covariates, dichotomous outcome)
  null_model <- tryCatch(
    SKAT_Null_Model(y ~ 1, out_type = "D"),
    error = function(e) {
      warning("SKAT_Null_Model (continuous) failed: ", conditionMessage(e))
      NULL
    }
  )

  n_sets       <- length(sets)
  cont_skat_p  <- rep(NA_real_, n_sets)
  cont_burden_p <- rep(NA_real_, n_sets)

  if (!is.null(null_model)) {
    for (s in seq_len(n_sets)) {
      cols <- sets[[s]]
      Z    <- z_mat[, cols, drop = FALSE]

      ## Remove zero-variance columns
      col_vars <- apply(Z, 2, var)
      Z <- Z[, col_vars > 1e-12, drop = FALSE]

      if (ncol(Z) < 1L) next

      ## SKAT (variance component test — good for heterogeneous effects)
      cont_skat_p[s] <- tryCatch({
        res <- SKAT(Z, null_model, kernel = "linear.weighted", method = "davies")
        res$p.value
      }, error = function(e) NA_real_)

      ## Burden (collapsing test — good for directional effects)
      cont_burden_p[s] <- tryCatch({
        res <- SKAT(Z, null_model, kernel = "linear.weighted",
                    method = "davies", r.corr = 1)
        res$p.value
      }, error = function(e) NA_real_)
    }
  }

  tibble(
    set_name      = names(sets),
    cont_skat_p   = cont_skat_p,
    cont_burden_p = cont_burden_p
  )
}


## ---------- 2e. Set-Level SKAT / Burden on BINARY hit matrix ----------------
##
## Repurposes the SKAT package by treating the binary hit matrix as a
## "genotype" matrix (0 = no hit, 1 = hit). The phenotype is case/control.
##
## Returns a tibble with one row per set and columns:
##   set_name, bin_skat_p, bin_burden_p

skat_burden_binary <- function(hits, y, sets) {

  null_model <- tryCatch(
    SKAT_Null_Model(y ~ 1, out_type = "D"),
    error = function(e) {
      warning("SKAT_Null_Model (binary) failed: ", conditionMessage(e))
      NULL
    }
  )

  n_sets      <- length(sets)
  bin_skat_p  <- rep(NA_real_, n_sets)
  bin_burden_p <- rep(NA_real_, n_sets)

  if (!is.null(null_model)) {
    for (s in seq_len(n_sets)) {
      cols <- sets[[s]]
      Z    <- hits[, cols, drop = FALSE]

      ## Remove zero-variance columns (proteins never hit or always hit)
      col_vars <- apply(Z, 2, var)
      Z <- Z[, col_vars > 0, drop = FALSE]

      if (ncol(Z) < 1L) next

      ## SKAT
      bin_skat_p[s] <- tryCatch({
        res <- SKAT(Z, null_model, kernel = "linear.weighted", method = "davies")
        res$p.value
      }, error = function(e) NA_real_)

      ## Burden
      bin_burden_p[s] <- tryCatch({
        res <- SKAT(Z, null_model, kernel = "linear.weighted",
                    method = "davies", r.corr = 1)
        res$p.value
      }, error = function(e) NA_real_)
    }
  }

  tibble(
    set_name     = names(sets),
    bin_skat_p   = bin_skat_p,
    bin_burden_p = bin_burden_p
  )
}


## ---------- 2f. Firth-Burden: sum-score per set → Firth logistic -----------
##
## For each set, collapse the binary hits into a single sum-score
## (count of hits per sample), then test the association between
## this sum-score and case/control status via Firth logistic regression.
##
## This is conceptually a "burden" test that does not assume a specific
## kernel structure and is robust to sparse binary data.
##
## Returns a tibble with one row per set and columns:
##   set_name, firth_burden_p

firth_burden_per_set <- function(hits, y, sets) {
  n_sets        <- length(sets)
  firth_burden_p <- rep(NA_real_, n_sets)

  for (s in seq_len(n_sets)) {
    cols      <- sets[[s]]
    sum_score <- rowSums(hits[, cols, drop = FALSE])

    ## Defensive: skip if sum_score is constant (no variation)
    if (var(sum_score) < 1e-12) next

    firth_burden_p[s] <- tryCatch({
      fit <- logistf(y ~ sum_score, firth = TRUE, plconf = NULL)
      fit$prob[2]  # p-value for sum_score coefficient
    }, error = function(e) NA_real_)
  }

  tibble(
    set_name       = names(sets),
    firth_burden_p = firth_burden_p
  )
}


## ---------- 2g. ACAT-V: Cauchy combination of per-protein p-values ---------
##
## For each set, aggregate the single-protein p-values (either limma or
## logistf) via the ACAT Cauchy combination method.
##
## Returns a tibble with one row per set and columns:
##   set_name, <output_col>

acatv_per_set <- function(protein_pvals, sets, output_col = "acatv_p") {
  n_sets  <- length(sets)
  acat_p  <- numeric(n_sets)

  for (s in seq_len(n_sets)) {
    cols     <- sets[[s]]
    acat_p[s] <- acat_combine(protein_pvals[cols])
  }

  tibble(
    set_name = names(sets),
    !!output_col := acat_p
  )
}


## ============================================================================
## ---------- 3. SINGLE-ITERATION RUNNER ----------
## ============================================================================
##
##  run_one_iteration()
##  -------------------
##  Simulates data at a given sample size, runs all tests, and returns a tidy
##  tibble of set-level results with signal annotations.

run_one_iteration <- function(iter_id = 1L, n_cases = 20L, n_controls = 20L) {

  ## -- 3a. Simulate data --
  sim <- simulate_huprot(n_cases = n_cases, n_controls = n_controls)

  ## -- 3b. Single-protein tests (full array) --
  limma_pvals <- limma_per_protein(sim$mfi, sim$group)
  firth_pvals <- firth_per_protein(sim$hits, sim$y)

  ## -- 3c. Set-level: Continuous SKAT & Burden (on Z-scores) --
  cont_res <- skat_burden_continuous(sim$z_mat, sim$y, sim$sets)

  ## -- 3d. Set-level: Binary SKAT & Burden (on hit matrix) --
  bin_res <- skat_burden_binary(sim$hits, sim$y, sim$sets)

  ## -- 3e. Set-level: Firth-Burden (sum-score of binary hits) --
  firth_burd_res <- firth_burden_per_set(sim$hits, sim$y, sim$sets)

  ## -- 3f. ACAT-V: Cauchy combination of per-protein p-values per set --
  acatv_limma_res <- acatv_per_set(limma_pvals, sim$sets,
                                   output_col = "acatv_limma_p")
  acatv_firth_res <- acatv_per_set(firth_pvals, sim$sets,
                                   output_col = "acatv_firth_p")

  ## -- 3g. Merge all set-level results --
  full_res <- cont_res %>%
    left_join(bin_res,          by = "set_name") %>%
    left_join(firth_burd_res,   by = "set_name") %>%
    left_join(acatv_limma_res,  by = "set_name") %>%
    left_join(acatv_firth_res,  by = "set_name")

  ## -- 3h. Annotate with signal metadata --
  n_sets <- length(sim$sets)

  full_res %>%
    mutate(
      is_signal = seq_len(n_sets) %in% sim$signal_ix,
      signal_type = case_when(
        seq_len(n_sets) %in% sim$signal_ix[sim$signal_type == "burden"] ~ "burden",
        seq_len(n_sets) %in% sim$signal_ix[sim$signal_type == "skat"]   ~ "skat",
        TRUE ~ "null"
      ),
      iteration = iter_id,
      n_cases   = n_cases,
      n_controls = n_controls
    )
}


## ============================================================================
## ---------- 4. GRID-BASED EXECUTION LOOP ----------
## ============================================================================

message("=======================================================================")
message("  HuProt Simulation Study v2")
message("  Grid: N_cases = N_controls in {",
        paste(N_GRID, collapse = ", "), "}")
message("  Iterations per grid point: ", N_ITER)
message("  Proteins: ", N_PROTEINS, " | Sets: ", N_SETS,
        " (", N_SIGNAL_SETS, " with signal)")
message("  Signal: burden_shift = ", BURDEN_SHIFT,
        ", skat_shift_sd = ", SKAT_SHIFT_SD)
message("  Started at ", Sys.time())
message("=======================================================================\n")

all_results <- vector("list", length(N_GRID) * N_ITER)
result_idx  <- 0L

for (n_val in N_GRID) {
  message("--- Grid point: N_cases = N_controls = ", n_val, " ---")
  t_start <- Sys.time()

  for (i in seq_len(N_ITER)) {
    if (i %% 50 == 0 || i == 1) {
      message("    Iteration ", i, " / ", N_ITER, " [", Sys.time(), "]")
    }

    result_idx <- result_idx + 1L
    all_results[[result_idx]] <- tryCatch(
      run_one_iteration(iter_id = i, n_cases = n_val, n_controls = n_val),
      error = function(e) {
        warning("N=", n_val, " iter=", i, " failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  message("    Completed N=", n_val, " in ",
          round(as.numeric(t_elapsed), 1), " minutes\n")
}

## Bind all iterations into a single tibble, dropping any that failed
results <- bind_rows(all_results)
message("\nSimulation complete at ", Sys.time())
message("Total set-level results: ", nrow(results), " rows\n")


## ============================================================================
## ---------- 5. EVALUATION: TYPE I ERROR & POWER ----------
## ============================================================================

## The test columns we want to evaluate
test_cols <- c(
  "cont_skat_p", "cont_burden_p",     # continuous SKAT/Burden
  "bin_skat_p", "bin_burden_p",        # binary SKAT/Burden
  "firth_burden_p",                    # Firth-Burden (sum-score)
  "acatv_limma_p", "acatv_firth_p"    # ACAT-V combinations
)

## Readable labels for plots
test_labels <- c(
  cont_skat_p     = "SKAT (continuous)",
  cont_burden_p   = "Burden (continuous)",
  bin_skat_p      = "SKAT (binary)",
  bin_burden_p    = "Burden (binary)",
  firth_burden_p  = "Firth-Burden",
  acatv_limma_p   = "ACAT-V (limma)",
  acatv_firth_p   = "ACAT-V (Firth)"
)

## -- 5a. Type I Error Rate --
## Proportion of null (non-signal) sets called significant at alpha,
## stratified by sample size.

null_results <- results %>% filter(!is_signal)

type1_error <- null_results %>%
  group_by(n_cases) %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(-n_cases, names_to = "test", values_to = "type1_error") %>%
  mutate(test_label = test_labels[test]) %>%
  arrange(n_cases, test)

message("== TYPE I ERROR RATE (alpha = ", ALPHA, ") ==")
message("   Expected under null: ", ALPHA)
message(strrep("-", 70))
for (n_val in N_GRID) {
  message("\n  N_cases = N_controls = ", n_val)
  sub <- type1_error %>% filter(n_cases == n_val)
  for (r in seq_len(nrow(sub))) {
    message(sprintf("    %-25s : %.4f", sub$test_label[r], sub$type1_error[r]))
  }
}
message("")

## -- 5b. Overall Power (across both signal types) --

signal_results <- results %>% filter(is_signal)

power_overall <- signal_results %>%
  group_by(n_cases) %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(-n_cases, names_to = "test", values_to = "power") %>%
  mutate(test_label = test_labels[test]) %>%
  arrange(n_cases, test)

message("== OVERALL STATISTICAL POWER (alpha = ", ALPHA, ") ==")
message(strrep("-", 70))
for (n_val in N_GRID) {
  message("\n  N_cases = N_controls = ", n_val)
  sub <- power_overall %>% filter(n_cases == n_val)
  for (r in seq_len(nrow(sub))) {
    message(sprintf("    %-25s : %.4f", sub$test_label[r], sub$power[r]))
  }
}
message("")

## -- 5c. Power stratified by signal type (Burden vs. SKAT) --

power_by_type <- signal_results %>%
  group_by(n_cases, signal_type) %>%
  summarise(across(all_of(test_cols),
                   ~ mean(.x < ALPHA, na.rm = TRUE),
                   .names = "{.col}"),
            .groups = "drop") %>%
  pivot_longer(all_of(test_cols), names_to = "test", values_to = "power") %>%
  mutate(test_label = test_labels[test]) %>%
  arrange(n_cases, signal_type, test)

message("== POWER BY SIGNAL TYPE ==")
message(strrep("-", 70))
for (n_val in N_GRID) {
  message("\n  N_cases = N_controls = ", n_val)
  sub_wide <- power_by_type %>%
    filter(n_cases == n_val) %>%
    select(test_label, signal_type, power) %>%
    pivot_wider(names_from = signal_type, values_from = power)
  message(sprintf("    %-25s   %-10s   %-10s", "Test", "Burden", "SKAT"))
  for (r in seq_len(nrow(sub_wide))) {
    message(sprintf("    %-25s : %.4f     %.4f",
                    sub_wide$test_label[r],
                    sub_wide$burden[r],
                    sub_wide$skat[r]))
  }
}
message("")


## ============================================================================
## ---------- 6. DIAGNOSTIC PLOTS ----------
## ============================================================================

tryCatch({

  ## 6a. Power Curves: Power vs. N for each test method (Burden signal)
  p_burden <- power_by_type %>%
    filter(signal_type == "burden") %>%
    ggplot(aes(x = n_cases, y = power, colour = test_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = N_GRID, trans = "log10") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Power vs. Sample Size — Burden Signal Architecture",
      subtitle = paste0(N_ITER, " iterations per grid point | ",
                        N_PROTEINS, " proteins | ", N_SETS, " sets | ",
                        "shift = ", BURDEN_SHIFT, " | alpha = ", ALPHA),
      x = expression(N[cases] == N[controls]),
      y = "Statistical Power",
      colour = "Test Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          legend.box = "horizontal")

  ggsave("power_curve_burden.pdf", p_burden, width = 10, height = 6)
  message("Saved: power_curve_burden.pdf")

  ## 6b. Power Curves: Power vs. N for each test method (SKAT signal)
  p_skat <- power_by_type %>%
    filter(signal_type == "skat") %>%
    ggplot(aes(x = n_cases, y = power, colour = test_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = N_GRID, trans = "log10") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Power vs. Sample Size — SKAT (Heterogeneous) Signal Architecture",
      subtitle = paste0(N_ITER, " iterations per grid point | ",
                        N_PROTEINS, " proteins | ", N_SETS, " sets | ",
                        "shift SD = ", SKAT_SHIFT_SD, " | alpha = ", ALPHA),
      x = expression(N[cases] == N[controls]),
      y = "Statistical Power",
      colour = "Test Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          legend.box = "horizontal")

  ggsave("power_curve_skat.pdf", p_skat, width = 10, height = 6)
  message("Saved: power_curve_skat.pdf")

  ## 6c. Type I Error vs. N
  p_t1e <- type1_error %>%
    ggplot(aes(x = n_cases, y = type1_error, colour = test_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = ALPHA, linetype = "dashed", colour = "grey40") +
    scale_x_continuous(breaks = N_GRID, trans = "log10") +
    labs(
      title = "Type I Error Rate vs. Sample Size",
      subtitle = paste0(N_ITER, " iterations per grid point | alpha = ", ALPHA),
      x = expression(N[cases] == N[controls]),
      y = "Empirical Type I Error Rate",
      colour = "Test Method"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          legend.box = "horizontal")

  ggsave("type1_error_vs_N.pdf", p_t1e, width = 10, height = 6)
  message("Saved: type1_error_vs_N.pdf")

  ## 6d. QQ-plot of null p-values per test, faceted by N
  qq_data <- null_results %>%
    select(set_name, iteration, n_cases, all_of(test_cols)) %>%
    pivot_longer(all_of(test_cols), names_to = "test", values_to = "pvalue") %>%
    filter(!is.na(pvalue)) %>%
    mutate(test_label = test_labels[test]) %>%
    group_by(test_label, n_cases) %>%
    arrange(pvalue) %>%
    mutate(expected = ppoints(n())) %>%
    ungroup()

  p_qq <- ggplot(qq_data, aes(x = -log10(expected), y = -log10(pvalue))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.15, size = 0.4, colour = "steelblue") +
    facet_grid(n_cases ~ test_label, scales = "free",
               labeller = labeller(n_cases = function(x) paste0("N=", x))) +
    labs(
      title = "QQ Plots of Null P-values (calibration check)",
      x = expression(-log[10](expected)),
      y = expression(-log[10](observed))
    ) +
    theme_minimal(base_size = 9)

  ggsave("null_qq_plots.pdf", p_qq, width = 16, height = 10)
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
