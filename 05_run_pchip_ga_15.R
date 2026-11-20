# Title: 15-Point PCHIP-Integrated Genetic Algorithm Schedule Generator

# Purpose: PCHIP is integrated in-loop as a reconstruction engine with this
# script generating multiple schedules per cohort then creating a consensus schedule
# based on the modal sampling times. All sampling frequencies
# included (sparse: 4, semi-sparse: 7, intensive: 15). The objective function 
# itself is the same one used in Pure GA. There is an optional conditional modal
# terminal anchor that can be enabled (enforce_terminal_anchor = FALSE for 
# manuscript benchmarks; configurable to TRUE for prospective protocol deployment).

# This script uses a limited 15 point sampling grid for PCHIP optimization as
# would be available in early clinical trials, and generated data used in
# Sensitivity Analysis 2 (comparison between dense 49-point PCHIP integrated GA
# and limited 15-point PCHIP integrated GA)

# Author: Amirtha Maria

# Input: simulated_PK_cohort_XX.csv (XX = 01 to 10)
# Output: all_best_schedules_pchip_15_binary.csv, sampling_frequency_pchip_15_summary.csv,
# selected_schedules_pchip_15_nca_metrics.csv

# Load packages ====

required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "GA", "doParallel", "foreach", "ncappc", "pracma", "gridExtra")
missing_packages  <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(GA); library(doParallel); library(foreach); library(ncappc); library(pracma); library(gridExtra)
})


# Configuration and multi-core cluster ====

input_directory <- "./data"
results_dir     <- "./schedules"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cohort_ids_to_run <- 1L:10L
ga_seeds          <- c(101L, 102L, 103L, 104L, 105L, 106L, 107L, 108L, 109L, 110L)
regime_configs    <- list("sparse" = 4L, "semi-sparse" = 7L, "intensive" = 15L)

# Optional Modal Terminal Anchor:
# Set to FALSE to reproduce manuscript benchmarks (unconstrained consensus);
# Set to TRUE to enforce a conditional modal terminal anchor within terminal_window.
enforce_terminal_anchor <- FALSE
terminal_window         <- c(72, 96)

# 15 canonical clinical observation times (Tsuchiwata and Tsuji 2023 Clinical Schedule)
clinical_15_target_times <- c(0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 72, 96)

total_cores <- parallel::detectCores()
n_workers   <- max(1L, total_cores - 1L)
cl <- makeCluster(n_workers)
registerDoParallel(cl)


# Data ingestion and grid extraction ====

load_cohort_csv <- function(c_id) {
  file_name <- sprintf("simulated_PK_cohort_%02d.csv", c_id)
  file_path <- file.path(input_directory, file_name)
  if (!file.exists(file_path)) stop("Missing file: ", file_path)
  
  raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  if ("NEW_ID" %in% names(raw)) {
    raw$ID <- as.integer(raw$NEW_ID)
  } else {
    raw$ID <- as.integer(raw$ID)
  }
  
  raw$TIME <- as.numeric(raw$TIME)
  raw$DV   <- as.numeric(raw$DV)
  
  raw %>% 
    dplyr::filter(TIME >= 0, TIME <= 96) %>%
    dplyr::select(ID, TIME, DV) %>% 
    dplyr::arrange(ID, TIME)
}

sample_data        <- load_cohort_csv(cohort_ids_to_run[1])
all_possible_times <- sort(unique(round(sample_data$TIME, 6))) # 49 candidate points
n_candidate_times  <- length(all_possible_times)

diff_limit <- 5; diff_lim_Tmax <- 50; PTWT_threshold <- 10; invalid_penalty <- 1e6

# PCHIP-15 fast evaluation engine ====

fast_pchip15_nca <- function(x_obs, y_obs, target_grid) {
  if (!any(round(x_obs, 6) == 0)) {
    x_obs <- c(0, x_obs)
    y_obs <- c(0, y_obs)
  }
  ord   <- order(x_obs)
  x_obs <- x_obs[ord]
  y_obs <- y_obs[ord]
  n     <- length(x_obs)
  if (n < 3L) return(NULL)
  
  max_x <- max(x_obs)
  in_t  <- target_grid[target_grid <= max_x]
  out_t <- target_grid[target_grid > max_x]
  
  in_y <- tryCatch(pracma::pchip(x_obs, y_obs, in_t), error = function(e) NULL)
  if (is.null(in_y) || anyNA(in_y) || any(!is.finite(in_y))) return(NULL)
  
  out_y <- numeric(0)
  if (length(out_t) > 0L) {
    slope <- (y_obs[n] - y_obs[n - 1L]) / (x_obs[n] - x_obs[n - 1L])
    out_y <- y_obs[n] + slope * (out_t - x_obs[n])
  }
  
  recon_y <- pmax(c(in_y, out_y), 0)
  recon_t <- c(in_t, out_t)
  
  # Prepend baseline (0,0) if target_grid doesn't include 0
  if (!any(round(recon_t, 6) == 0)) {
    recon_t <- c(0, recon_t)
    recon_y <- c(0, recon_y)
  }
  
  cmax_idx <- which.max(recon_y)
  n_g      <- length(recon_t)
  auc      <- sum((recon_t[-1] - recon_t[-n_g]) * (recon_y[-1] + recon_y[-n_g]) / 2)
  
  c(cmax = recon_y[cmax_idx], tmax = recon_t[cmax_idx], auc = auc)
}

run_ncappc <- function(data) {
  subj_ids <- sort(unique(data$ID))
  has_t0   <- data %>% dplyr::filter(round(TIME, 6) == 0)
  if (length(unique(has_t0$ID)) < length(subj_ids)) {
    missing_t0 <- setdiff(subj_ids, has_t0$ID)
    baseline   <- data.frame(ID = missing_t0, TIME = 0, DV = 0)
    data       <- dplyr::bind_rows(data, baseline) %>% dplyr::arrange(ID, TIME)
  }
  as.data.frame(ncappc::ncappc(
    obsFile = data, concNmObs = "DV", method = "linearup-logdown",
    onlyNCA = TRUE, extrapolate = FALSE, printOut = FALSE, evid = FALSE, noPlot = TRUE
  ))
}

abs_mpe <- function(pred, true) abs(mean((pred - true) / true) * 100)
rmspe   <- function(pred, true) sqrt(mean(((pred - true) / true)^2)) * 100

calculate_standard_metrics <- function(cmax, auc, tmax, true_cmax, true_auc, true_tmax) {
  c(
    abs_MPE_Cmax = abs_mpe(cmax, true_cmax), RMSPE_Cmax = rmspe(cmax, true_cmax),
    abs_MPE_AUC  = abs_mpe(auc, true_auc),   RMSPE_AUC  = rmspe(auc, true_auc),
    abs_MPE_Tmax = abs_mpe(tmax, true_tmax), RMSPE_Tmax = rmspe(tmax, true_tmax)
  )
}

make_eval_pchip15 <- function(out, true_cmax, true_auc, true_tmax, n_subjects, target_n) {
  split_data <- split(out, out$ID)
  
  function(chromosome) {
    chromosome <- as.integer(chromosome)
    n_sel <- sum(chromosome == 1L)
    if (n_sel < 2L) return(invalid_penalty)
    
    sel_times <- all_possible_times[chromosome == 1L]
    cand_tmax <- numeric(n_subjects); cand_cmax <- numeric(n_subjects); cand_auc <- numeric(n_subjects)
    
    for (k in seq_len(n_subjects)) {
      sub_df <- split_data[[k]]
      m_idx  <- round(sub_df$TIME, 6) %in% round(sel_times, 6)
      if (sum(m_idx) < 2L) return(invalid_penalty)
      
      # Reconstruct onto the 15-point clinical grid
      res <- fast_pchip15_nca(sub_df$TIME[m_idx], sub_df$DV[m_idx], clinical_15_target_times)
      if (is.null(res)) return(invalid_penalty)
      
      cand_cmax[k] <- res["cmax"]; cand_tmax[k] <- res["tmax"]; cand_auc[k] <- res["auc"]
    }
    
    mpe_cmax <- abs_mpe(cand_cmax, true_cmax); rmspe_cmax <- rmspe(cand_cmax, true_cmax)
    mpe_auc  <- abs_mpe(cand_auc, true_auc);   rmspe_auc  <- rmspe(cand_auc, true_auc)
    mpe_tmax <- abs_mpe(cand_tmax, true_tmax); rmspe_tmax <- rmspe(cand_tmax, true_tmax)
    
    if (anyNA(c(mpe_cmax, rmspe_cmax, mpe_auc, rmspe_auc, mpe_tmax, rmspe_tmax))) return(invalid_penalty)
    
    # Objective Function Design Rationale:
    # 25% MPE / 75% RMSPE weighting suppresses inter-individual outlier dispersion.
    # wt_Tmax = 0.025 prevents peak-chasing along flat derivatives near Cmax.
    # Continuous quadratic barrier replaces step-cliffs to maintain gradient selection pressure.
    # Quadratic barrier strictly penalizes deviation from target_n
    # MPE/RMSPE of Cmax, AUC, and Tmax normalized by clinical tolerance threshold
    # Note: PCHIP is not a part of the objective function, but is integrated in the overall
    # optimization loop
    
    obj <- 100000 * ((n_sel - target_n) / PTWT_threshold)^2 +
      0.25 * 0.60 * (mpe_cmax / diff_limit) + 0.25 * 0.80 * (mpe_auc / diff_limit) + 0.25 * 0.025 * (mpe_tmax / diff_lim_Tmax) +
      0.75 * 0.60 * (rmspe_cmax / diff_limit) + 0.75 * 0.80 * (rmspe_auc / diff_limit) + 0.75 * 0.025 * (rmspe_tmax / diff_lim_Tmax)
    
    ratios <- c(mpe_cmax/diff_limit, mpe_auc/diff_limit, mpe_tmax/diff_lim_Tmax,
                rmspe_cmax/diff_limit, rmspe_auc/diff_limit, rmspe_tmax/diff_lim_Tmax)
    as.numeric(obj + 100 * sum(pmax(ratios - 1, 0)^2))
  }
}


# Reference NCA cache ====

cohort_cache <- list()
for (c_id in cohort_ids_to_run) {
  c_data   <- load_cohort_csv(c_id)
  nca_full <- run_ncappc(c_data)
  cohort_cache[[as.character(c_id)]] <- list(
    data       = c_data,
    true_cmax  = nca_full[, 6],
    true_auc   = nca_full[, 11],
    true_tmax  = nca_full[, 4],
    n_subjects = length(unique(c_data$ID))
  )
}

# Master PCHIP-15 and GA Loop ==== 

all_regimes_binary <- list(); all_regimes_freq <- list(); all_regimes_metrics <- list()

for (regime_name in names(regime_configs)) {
  target_n <- regime_configs[[regime_name]]
  message(sprintf("\n=== PCHIP-15 GA Optimization: %s (Total Samples = %d) ===", toupper(regime_name), target_n))
  
  tasks_df <- expand.grid(cohort_id = cohort_ids_to_run, seed = ga_seeds)
  
  ga_results <- foreach(
    row_i = seq_len(nrow(tasks_df)),
    .packages = c("GA", "dplyr", "pracma"),
    .export = c("all_possible_times", "n_candidate_times", "clinical_15_target_times",
                "diff_limit", "diff_lim_Tmax", "PTWT_threshold", "invalid_penalty",
                "fast_pchip15_nca", "cohort_cache", "make_eval_pchip15", "abs_mpe", "rmspe")
  ) %dopar% {
    cur_c  <- tasks_df$cohort_id[row_i]; cur_s <- tasks_df$seed[row_i]
    c_info <- cohort_cache[[as.character(cur_c)]]
    fit_fn <- make_eval_pchip15(c_info$data, c_info$true_cmax, c_info$true_auc, c_info$true_tmax, c_info$n_subjects, target_n)
    
    set.seed(cur_s)
    ga_mod <- GA::ga(
      type = "binary", nBits = n_candidate_times,
      fitness = function(chr) -fit_fn(chr),
      popSize = 200, maxiter = 100, pmutation = 0.03, elitism = 40, parallel = FALSE
    )
    list(cohort_id = cur_c, seed = cur_s, chromosome = as.integer(ga_mod@solution[1, ]))
  }
  
  regime_binary_list <- list()
  for (c_id in cohort_ids_to_run) {
    c_runs    <- ga_results[sapply(ga_results, function(x) x$cohort_id == c_id)]
    chrom_mat <- do.call(rbind, lapply(c_runs, function(x) x$chromosome))
    
    counts <- colSums(chrom_mat)
    freqs  <- counts / length(ga_seeds)
    all_regimes_freq[[length(all_regimes_freq) + 1]] <- data.frame(
      DENSITY_REGIME = regime_name, COHORT = c_id, TIME = all_possible_times,
      SELECTION_COUNT = counts, SELECTION_FREQ = freqs
    )
    
    # Rank all candidate timepoints by cross-seed selection frequency (tie-breaker: earlier time)
    ranked_idx   <- order(counts, all_possible_times, decreasing = c(TRUE, FALSE))
    top_cand_idx <- ranked_idx[1:target_n]
    
    # Optional terminal anchoring: enforce modal terminal anchor only when enabled
    if (enforce_terminal_anchor &&
        !any(all_possible_times[top_cand_idx] >= terminal_window[1] &
             all_possible_times[top_cand_idx] <= terminal_window[2])) {
      term_idx               <- which(all_possible_times >= terminal_window[1] &
                                        all_possible_times <= terminal_window[2])
      modal_term_idx         <- term_idx[order(counts[term_idx], all_possible_times[term_idx],
                                               decreasing = c(TRUE, FALSE))[1]]
      top_cand_idx[target_n] <- modal_term_idx
    }
    final_design <- sort(all_possible_times[top_cand_idx])
    
    binary_row <- as.data.frame(t(as.integer(all_possible_times %in% final_design)))
    colnames(binary_row) <- paste0("T_", all_possible_times, "h")
    binary_row <- cbind(COHORT = c_id, DENSITY_REGIME = regime_name, N_SAMPLES = length(final_design), binary_row)
    
    regime_binary_list[[as.character(c_id)]] <- binary_row
    all_regimes_binary[[length(all_regimes_binary) + 1]] <- binary_row
    
    # Validation against ncappc
    c_entry    <- cohort_cache[[as.character(c_id)]]
    split_subj <- split(c_entry$data, c_entry$data$ID)
    cand_cmax  <- numeric(length(split_subj)); cand_tmax <- numeric(length(split_subj)); cand_auc <- numeric(length(split_subj))
    for (s_i in seq_along(split_subj)) {
      sub_d <- split_subj[[s_i]]
      m_i   <- round(sub_d$TIME, 6) %in% round(final_design, 6)
      res   <- fast_pchip15_nca(sub_d$TIME[m_i], sub_d$DV[m_i], clinical_15_target_times)
      cand_cmax[s_i] <- res["cmax"]; cand_tmax[s_i] <- res["tmax"]; cand_auc[s_i] <- res["auc"]
    }
    val_met <- calculate_standard_metrics(cand_cmax, cand_auc, cand_tmax, c_entry$true_cmax, c_entry$true_auc, c_entry$true_tmax)
    all_regimes_metrics[[length(all_regimes_metrics) + 1]] <- data.frame(DENSITY_REGIME = regime_name, COHORT = c_id, t(val_met))
  }
  
  write.csv(dplyr::bind_rows(regime_binary_list), file.path(results_dir, sprintf("best_schedules_pchip15_%s_binary.csv", regime_name)), row.names = FALSE)
}

stopCluster(cl)

# Data exports ====

write.csv(dplyr::bind_rows(all_regimes_binary),  file.path(results_dir, "all_best_schedules_pchip_15_binary.csv"),   row.names = FALSE)
write.csv(dplyr::bind_rows(all_regimes_freq),    file.path(results_dir, "sampling_frequency_pchip_15_summary.csv"), row.names = FALSE)
write.csv(dplyr::bind_rows(all_regimes_metrics), file.path(results_dir, "selected_schedules_pchip_15_nca_metrics.csv"), row.names = FALSE)
message(">> PCHIP-15 GA pipeline completed successfully.")