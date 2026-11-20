# Title: Tsuchiwata Objective Function Schedule Generator

# Purpose: This script uses Tsuchiwata and Tsuji (2023)'s genetic algorithm
# approach and objective function to generate schedules from the same 10 standard
# cohorts used to generate other schedules for this project.
# These outputs were used for benchmarking of the new objective function (used in
# Pure GA) against Tsuchiwata and Tsuji's historical formulation in a paired 
# cohort-matched analyses against schedules generated on the identical 10 training cohorts.

# As Tsuchiwata and Tsuji's formulation was designed to target about 7 sampling
# times, only semi-sparse (N=7) schedules were generated for comparison to 
# semi-sparse Pure GA schedules. However, the modal consensus protocol to select
# the top N sampling times as the final schedule for each cohort was kept in order
# to have a fair comparison as the original Tsuchiwata and Tsuji approach meant 
# that final schedules could have more or less sampling times than the target.

# Author: Amirtha Maria

# Input: simulated_PK_cohort_XX.csv (XX = 01 to 10)
# Output: all_best_schedules_tsuchiwata_binary.csv, sampling_frequency_tsuchiwata_summary.csv,
# selected_schedules_tsuchiwata_nca_metrics.csv 

# Load packages ====

required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "scales", "GA", 
  "doParallel", "foreach", "ncappc", "gridExtra"
)
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

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

# Primary focus is Semi-Sparse (N = 7), matching Tsuchiwata's PTWT = 7
regime_configs    <- list("semi-sparse" = 7L)

total_cores <- parallel::detectCores()
n_workers   <- max(1L, total_cores - 1L)
message(sprintf(">> Starting Tsuchiwata Ablation Cluster with %d workers...", n_workers))

cl <- makeCluster(n_workers)
registerDoParallel(cl)

# Data ingestion and grid extraction ====

load_cohort_csv <- function(c_id) {
  file_name <- sprintf("simulated_PK_cohort_%02d.csv", c_id)
  file_path <- file.path(input_directory, file_name)
  if (!file.exists(file_path)) stop("Missing file: ", file_path)
  
  raw <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Prioritize NEW_ID (1-24) to avoid collisions from sampling with replacement
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

# Fast in-memory NCA engine ====

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
  subj_ids <- sort(unique(data$ID))
  has_t0   <- data %>% dplyr::filter(round(TIME, 6) == 0)
  if (length(unique(has_t0$ID)) < length(subj_ids)) {
    missing_t0 <- setdiff(subj_ids, has_t0$ID)
    baseline   <- data.frame(ID = missing_t0, TIME = 0, DV = 0)
    data       <- dplyr::bind_rows(data, baseline) %>% dplyr::arrange(ID, TIME)
  }
  as.data.frame(ncappc::ncappc(
    obsFile     = data, 
    concNmObs   = "DV", 
    method      = "linearup-logdown",
    onlyNCA     = TRUE, 
    extrapolate = FALSE, 
    printOut    = FALSE, 
    evid        = FALSE, 
    noPlot      = TRUE
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

# Implementation of the Tsuchiwata and Tsuji (2023) methodology ====

make_eval_tsuchiwata <- function(out, true_cmax, true_auc, n_subjects, target_n = 7) {
  split_data <- split(out, out$ID)
  
  function(chromosome) {
    chromosome <- as.integer(chromosome)
    n_sel <- sum(chromosome == 1L)
    
    # Mathematical minimum to calculate trapezoidal area
    if (n_sel < 2L) return(1000)
    
    sel_times <- all_possible_times[chromosome == 1L]
    cand_cmax <- numeric(n_subjects)
    cand_auc  <- numeric(n_subjects)
    
    for (k in seq_len(n_subjects)) {
      sub_df <- split_data[[k]]
      m_idx  <- round(sub_df$TIME, 6) %in% round(sel_times, 6)
      if (sum(m_idx) < 2L) return(1000)
      
      res <- fast_nca(sub_df$TIME[m_idx], sub_df$DV[m_idx])
      cand_cmax[k] <- res["cmax"]
      cand_auc[k]  <- res["auc"]
    }
    
    # Exact Tsuchiwata |MPE| formulation (labeled MAPE in paper)
    MAPE_Cmax  <- abs(mean((cand_cmax - true_cmax) / true_cmax)) * 100
    MAPE_AUC   <- abs(mean((cand_auc  - true_auc)  / true_auc))  * 100
    
    # Exact Tsuchiwata RMSPE formula (* 100 inside square root)
    RMSPE_Cmax <- (mean(((cand_cmax - true_cmax) / true_cmax)^2) * 100)^0.5
    RMSPE_AUC  <- (mean(((cand_auc  - true_auc)  / true_auc)^2)  * 100)^0.5
    
    if (anyNA(c(MAPE_Cmax, MAPE_AUC, RMSPE_Cmax, RMSPE_AUC))) return(1000)
    
    # Exact Tsuchiwata soft sample penalty + raw unweighted sum
    value <- (n_sel / target_n)^2 + MAPE_Cmax + MAPE_AUC + RMSPE_Cmax + RMSPE_AUC
    
    # Exact Tsuchiwata 5% step-cliff constraint (Diff_limit = 5.0)
    if (MAPE_Cmax <= 5.0 && MAPE_AUC <= 5.0 && RMSPE_Cmax <= 5.0 && RMSPE_AUC <= 5.0) {
      return(as.numeric(value))
    } else {
      return(1000)
    }
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

# Master Tsuchiwata optimization loop ====

all_regimes_binary  <- list()
all_regimes_freq    <- list()
all_regimes_metrics <- list()

for (regime_name in names(regime_configs)) {
  target_n <- regime_configs[[regime_name]]
  message(sprintf("\n=== Tsuchiwata GA Ablation: %s (PTWT = %d) ===", toupper(regime_name), target_n))
  
  tasks_df <- expand.grid(cohort_id = cohort_ids_to_run, seed = ga_seeds)
  
  ga_results <- foreach(
    row_i = seq_len(nrow(tasks_df)),
    .packages = c("GA", "dplyr"),
    .export = c("all_possible_times", "n_candidate_times", "fast_nca", 
                "cohort_cache", "make_eval_tsuchiwata")
  ) %dopar% {
    cur_c  <- tasks_df$cohort_id[row_i]
    cur_s  <- tasks_df$seed[row_i]
    c_info <- cohort_cache[[as.character(cur_c)]]
    
    # Initialize Tsuchiwata objective function
    fit_fn <- make_eval_tsuchiwata(
      out         = c_info$data, 
      true_cmax   = c_info$true_cmax, 
      true_auc    = c_info$true_auc, 
      n_subjects  = c_info$n_subjects, 
      target_n    = target_n
    )
    
    set.seed(cur_s)
    ga_mod <- GA::ga(
      type           = "binary", 
      nBits          = n_candidate_times,
      fitness        = function(chr) -fit_fn(chr), # GA maximizes, so negate cost
      popSize        = 200, 
      maxiter        = 100, 
      pmutation      = 0.03, 
      elitism        = 40, 
      parallel       = FALSE
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
      DENSITY_REGIME  = regime_name, 
      COHORT          = c_id, 
      TIME            = all_possible_times,
      SELECTION_COUNT = counts, 
      SELECTION_FREQ  = freqs
    )
    
    # Select the consensus top-N points across seeds
    top_cand_idx <- order(counts, all_possible_times, decreasing = c(TRUE, FALSE))[1:target_n]
    final_design <- sort(all_possible_times[top_cand_idx])
    
    binary_row <- as.data.frame(t(as.integer(all_possible_times %in% final_design)))
    colnames(binary_row) <- paste0("T_", all_possible_times, "h")
    binary_row <- cbind(
      COHORT         = c_id, 
      DENSITY_REGIME = regime_name, 
      N_SAMPLES      = length(final_design), 
      binary_row
    )
    
    regime_binary_list[[as.character(c_id)]] <- binary_row
    all_regimes_binary[[length(all_regimes_binary) + 1]] <- binary_row
    
    # Validation against ncappc
    c_entry    <- cohort_cache[[as.character(c_id)]]
    sparse_obs <- c_entry$data %>% dplyr::filter(round(TIME, 6) %in% round(final_design, 6))
    nca_cand   <- run_ncappc(sparse_obs)
    val_met    <- calculate_standard_metrics(
      cmax      = nca_cand[, 6], 
      auc       = nca_cand[, 11], 
      tmax      = nca_cand[, 4],
      true_cmax = c_entry$true_cmax, 
      true_auc  = c_entry$true_auc, 
      true_tmax = c_entry$true_tmax
    )
    all_regimes_metrics[[length(all_regimes_metrics) + 1]] <- data.frame(
      DENSITY_REGIME = regime_name, 
      COHORT         = c_id, 
      t(val_met)
    )
  }
  
  write.csv(
    dplyr::bind_rows(regime_binary_list), 
    file.path(results_dir, sprintf("best_schedules_tsuchiwata_%s_binary.csv", regime_name)), 
    row.names = FALSE
  )
}

stopCluster(cl)

# Data exports and visualizations ====

master_binary_df <- dplyr::bind_rows(all_regimes_binary)
master_freq_df   <- dplyr::bind_rows(all_regimes_freq)
master_metric_df <- dplyr::bind_rows(all_regimes_metrics)

write.csv(master_binary_df, file.path(results_dir, "all_best_schedules_tsuchiwata_binary.csv"), row.names = FALSE)
write.csv(master_freq_df,   file.path(results_dir, "sampling_frequency_tsuchiwata_summary.csv"), row.names = FALSE)
write.csv(master_metric_df, file.path(results_dir, "selected_schedules_tsuchiwata_nca_metrics.csv"), row.names = FALSE)

theme_set(theme_minimal(base_size = 12))

long_master <- master_binary_df %>%
  pivot_longer(cols = starts_with("T_"), names_to = "T_Label", values_to = "Selected") %>%
  mutate(
    Time   = as.numeric(sub("^T_(.+)h$", "\\1", T_Label)),
    Status = factor(Selected, levels = c(0, 1), labels = c("Not Selected", "Selected"))
  )

for (r_name in names(regime_configs)) {
  sub_binary <- long_master %>% dplyr::filter(DENSITY_REGIME == r_name)
  sub_freq   <- master_freq_df %>% dplyr::filter(DENSITY_REGIME == r_name)
  
  p_final <- ggplot(sub_binary, aes(x = factor(Time), y = factor(COHORT), fill = Status)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Not Selected" = "#F0F2F5", "Selected" = "#D55E00")) +
    labs(
      title    = sprintf("Tsuchiwata Ablation Schedules (%s Regime: Target N=%d)", toupper(r_name), regime_configs[[r_name]]),
      subtitle = "Optimization using Tsuchiwata et al. (2023) Objective Formulation",
      x        = "Sampling Time (h)", 
      y        = "Cohort ID"
    ) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  
  ggsave(file.path(results_dir, sprintf("plot_final_schedule_only_tsuchiwata_%s.png", r_name)), 
         p_final, width = 12, height = 5.5, dpi = 300)
}

p_density <- ggplot(master_freq_df, aes(x = TIME, y = SELECTION_FREQ, color = DENSITY_REGIME, group = DENSITY_REGIME)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  scale_color_manual(values = c("semi-sparse" = "#D55E00")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = all_possible_times[seq(1, length(all_possible_times), by = 2)]) +
  labs(
    title    = "Mean Selection Pressure across 96h Grid (Tsuchiwata Formulation)",
    subtitle = "Averaged across 10 cohorts using Tsuchiwata et al. (2023) objective function",
    x        = "Sampling Time (h)", 
    y        = "Mean Selection Frequency", 
    color    = "Regime"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

ggsave(file.path(results_dir, "plot_sampling_density_profiles_tsuchiwata.png"), p_density, width = 12, height = 5.5, dpi = 300)

message(">> Tsuchiwata ablation GA pipeline completed successfully.")