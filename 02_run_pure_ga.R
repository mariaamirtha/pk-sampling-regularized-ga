# Title: Pure Genetic Algorithm Schedule Generator

# Purpose: This script uses a pure genetic algorithm approach with the new 
# objective function integrating the continuous quadratic barrier, new weighting
# for MPE and RMSPE, new weighting for Cmax and AUC, the inclusion of Tmax, 
# and the normalization of PK metrics by clinical tolerance threshold.
# Multiple schedules per cohort were generated, then a consensus schedule
# based on the modal sampling times was created.All sampling frequencies included 
# (sparse: 4, semi-sparse: 7, intensive: 15).There is an optional conditional modal
# terminal anchor that can be enabled (enforce_terminal_anchor = FALSE for 
# manuscript benchmarks; configurable to TRUE for prospective protocol deployment).

# Author: Amirtha Maria

# Input: simulated_PK_cohort_XX.csv (XX = 01 to 10)
# Output: all_best_schedules_pure_ga_binary.csv, sampling_frequency_pure_ga_summary.csv,
# selected_schedules_pure_ga_nca_metrics.csv 

# Load packages ====

required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "GA", "doParallel", "foreach", "ncappc", "gridExtra")
missing_packages  <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(GA); library(doParallel); library(foreach); library(ncappc); library(gridExtra)
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

total_cores <- parallel::detectCores()
n_workers   <- max(1L, total_cores - 1L)
message(sprintf(">> Starting Pure GA cluster on %d worker threads...", n_workers))

cl <- makeCluster(n_workers)
registerDoParallel(cl)

# Data ingestion and full grid extraction ====

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
all_possible_times <- sort(unique(round(sample_data$TIME, 6))) # All 49 candidate points
n_candidate_times  <- length(all_possible_times)

# Candidate times for GA chromosome (all 49 points are eligible)
diff_limit <- 5; diff_lim_Tmax <- 50; PTWT_threshold <- 10; invalid_penalty <- 1e6

# Fast in-memory NCA engine with prepended C(0)=0 for interpolation ====
fast_nca <- function(time, conc) {
  # Prepend physical baseline (t=0, C=0) if not explicitly drawn
  if (!any(round(time, 6) == 0)) {
    time <- c(0, time)
    conc <- c(0, conc)
  }
  ord <- order(time); t <- time[ord]; c <- conc[ord]; n <- length(t)
  cmax_idx <- which.max(c)
  auc      <- sum((t[-1] - t[-n]) * (c[-1] + c[-n]) / 2)
  c(cmax = c[cmax_idx], tmax = t[cmax_idx], auc = auc)
}

run_ncappc <- function(data) {
  # Ensure baseline t=0 exists per subject for reference calculation
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

make_eval_function <- function(out, true_cmax, true_auc, true_tmax, n_subjects, target_n) {
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
      
      res <- fast_nca(sub_df$TIME[m_idx], sub_df$DV[m_idx])
      cand_cmax[k] <- res["cmax"]; cand_tmax[k] <- res["tmax"]; cand_auc[k] <- res["auc"]
    }
    
    mpe_cmax <- abs_mpe(cand_cmax, true_cmax); rmspe_cmax <- rmspe(cand_cmax, true_cmax)
    mpe_auc  <- abs_mpe(cand_auc, true_auc);   rmspe_auc  <- rmspe(cand_auc, true_auc)
    mpe_tmax <- abs_mpe(cand_tmax, true_tmax); rmspe_tmax <- rmspe(cand_tmax, true_tmax)
    
    if (anyNA(c(mpe_cmax, rmspe_cmax, mpe_auc, rmspe_auc, mpe_tmax, rmspe_tmax))) return(invalid_penalty)
    
    # Objective Function Design Rationale:
    # 25% MPE / 75% RMSPE weighting suppresses inter-individual outlier dispersion.
    # weight of Tmax = 0.025 prevents peak-chasing along flat derivatives near Cmax.
    # Continuous quadratic barrier replaces step-cliffs to maintain gradient selection pressure.
    # Quadratic barrier strictly penalizes deviation from target_n
    # MPE/RMSPE of Cmax, AUC, and Tmax normalized by clinical tolerance threshold
    
    obj <- 100000 * ((n_sel - target_n) / PTWT_threshold)^2 +
      0.25 * 0.60 * (mpe_cmax / diff_limit) + 0.25 * 0.80 * (mpe_auc / diff_limit) + 0.25 * 0.025 * (mpe_tmax / diff_lim_Tmax) +
      0.75 * 0.60 * (rmspe_cmax / diff_limit) + 0.75 * 0.80 * (rmspe_auc / diff_limit) + 0.75 * 0.025 * (rmspe_tmax / diff_lim_Tmax)
    
    ratios <- c(mpe_cmax/diff_limit, mpe_auc/diff_limit, mpe_tmax/diff_lim_Tmax,
                rmspe_cmax/diff_limit, rmspe_auc/diff_limit, rmspe_tmax/diff_lim_Tmax)
    as.numeric(obj + 100 * sum(pmax(ratios - 1, 0)^2))
  }
}

# Pre-compute ground-truth reference NCA cache ====

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

# Master Pure GA optimization loop ====

all_regimes_binary <- list(); all_regimes_freq <- list(); all_regimes_metrics <- list()

for (regime_name in names(regime_configs)) {
  target_n <- regime_configs[[regime_name]]
  message(sprintf("\n=== Pure GA Optimization: %s (Total Samples = %d) ===", toupper(regime_name), target_n))
  
  tasks_df <- expand.grid(cohort_id = cohort_ids_to_run, seed = ga_seeds)
  
  ga_results <- foreach(
    row_i = seq_len(nrow(tasks_df)),
    .packages = c("GA", "dplyr"),
    .export = c("all_possible_times", "n_candidate_times", "diff_limit", "diff_lim_Tmax", 
                "PTWT_threshold", "invalid_penalty", "fast_nca", "cohort_cache", 
                "make_eval_function", "abs_mpe", "rmspe")
  ) %dopar% {
    cur_c  <- tasks_df$cohort_id[row_i]; cur_s <- tasks_df$seed[row_i]
    c_info <- cohort_cache[[as.character(cur_c)]]
    fit_fn <- make_eval_function(c_info$data, c_info$true_cmax, c_info$true_auc, c_info$true_tmax, c_info$n_subjects, target_n)
    
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
    sparse_obs <- c_entry$data %>% dplyr::filter(round(TIME, 6) %in% round(final_design, 6))
    nca_cand   <- run_ncappc(sparse_obs)
    val_met    <- calculate_standard_metrics(nca_cand[, 6], nca_cand[, 11], nca_cand[, 4],
                                             c_entry$true_cmax, c_entry$true_auc, c_entry$true_tmax)
    all_regimes_metrics[[length(all_regimes_metrics) + 1]] <- data.frame(DENSITY_REGIME = regime_name, COHORT = c_id, t(val_met))
  }
  
  write.csv(dplyr::bind_rows(regime_binary_list), file.path(results_dir, sprintf("best_schedules_%s_binary.csv", regime_name)), row.names = FALSE)
}

stopCluster(cl)

# Data exports and visualizations ====

master_binary_df <- dplyr::bind_rows(all_regimes_binary)
master_freq_df   <- dplyr::bind_rows(all_regimes_freq)
master_metric_df <- dplyr::bind_rows(all_regimes_metrics)

write.csv(master_binary_df, file.path(results_dir, "all_best_schedules_pure_ga_binary.csv"), row.names = FALSE)
write.csv(master_freq_df,   file.path(results_dir, "sampling_frequency_pure_ga_summary.csv"), row.names = FALSE)
write.csv(master_metric_df, file.path(results_dir, "selected_schedules_pure_ga_nca_metrics.csv"), row.names = FALSE)

theme_set(theme_minimal(base_size = 12))

long_master <- master_binary_df %>%
  pivot_longer(cols = starts_with("T_"), names_to = "T_Label", values_to = "Selected") %>%
  mutate(Time = as.numeric(sub("^T_(.+)h$", "\\1", T_Label)),
         Status = factor(Selected, levels = c(0, 1), labels = c("Not Selected", "Selected")))

for (r_name in names(regime_configs)) {
  sub_binary <- long_master %>% dplyr::filter(DENSITY_REGIME == r_name)
  sub_freq   <- master_freq_df %>% dplyr::filter(DENSITY_REGIME == r_name)
  
  p_final <- ggplot(sub_binary, aes(x = factor(Time), y = factor(COHORT), fill = Status)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Not Selected" = "#F0F2F5", "Selected" = "#2E7D32")) +
    labs(title = sprintf("Optimal Pure GA Schedules (%s Regime: Strict N=%d)", toupper(r_name), regime_configs[[r_name]]),
         x = "Sampling Time (h)", y = "Cohort ID") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  
  ggsave(file.path(results_dir, sprintf("plot_final_schedule_only_%s.png", r_name)), p_final, width = 12, height = 5.5, dpi = 300)
}

p_density <- ggplot(master_freq_df, aes(x = TIME, y = SELECTION_FREQ, color = DENSITY_REGIME, group = DENSITY_REGIME)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_color_manual(values = c("sparse" = "#E65100", "semi-sparse" = "#1976D2", "intensive" = "#388E3C")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = all_possible_times[seq(1, length(all_possible_times), by = 2)]) +
  labs(title = "Mean Pure GA Selection Pressure Across 96h Grid (Strict N=4, 7, 15)",
       x = "Sampling Time (h)", y = "Mean Selection Frequency", color = "Regime") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

ggsave(file.path(results_dir, "plot_sampling_density_profiles_pure_ga.png"), p_density, width = 12, height = 5.5, dpi = 300)
message(">> Pure GA pipeline completed successfully.")