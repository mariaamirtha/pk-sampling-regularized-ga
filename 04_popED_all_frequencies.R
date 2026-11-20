# Title: PopED Schedule Generation Across All Sampling Densities
# Purpose: Generate a PopED schedule for the 10 cohorts at N=4, 7, 15 sampling
# frequencies and output results in format acceptable for data analysis

# Author: Amirtha Maria

# Inputs: simulated_PK_cohort_XX.csv (XX = 01 to 10)

# Outputs: all_best_schedules_poped_binary.csv, selected_schedules_poped_nca_metrics.csv

# Load packages ====

required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "PopED", "doParallel", "foreach", "ncappc", "gridExtra")
missing_packages  <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) install.packages(missing_packages, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales)
  library(PopED); library(doParallel); library(foreach); library(ncappc); library(gridExtra)
})

# Configuration and multi-core setup ====
input_directory <- "./data"
results_dir     <- "./schedules"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

cohort_ids_to_run <- 1L:10L
regime_configs    <- list("sparse" = 4L, "semi-sparse" = 7L, "intensive" = 15L)

total_cores <- parallel::detectCores()
n_workers   <- min(length(cohort_ids_to_run), max(1L, total_cores - 1L))
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
    dplyr::select(ID, TIME, DV, any_of(c("TVCL", "TVV2", "TVV3", "TVQ", "TVKA"))) %>% 
    dplyr::arrange(ID, TIME)
}

sample_data        <- load_cohort_csv(cohort_ids_to_run[1])
all_possible_times <- sort(unique(round(sample_data$TIME, 6)))

# Candidate discrete pool for PopED: times after lag window (> 0.125h)
discrete_space <- all_possible_times[all_possible_times > 0.125]

# PopED 2-compartment oral model (with lag time) ====

ff_2comp_oral <- function(model_switch, xt, p, poped.db) {
  CL <- p[1]
  V2 <- p[2]
  V3 <- p[3]
  Q <- p[4]
  KA <- p[5]
  DOSE  <- 10000
  ALAG1 <- 0.125
  
  # Micro and macro rate constants
  k10 <- CL / V2
  k12 <- Q / V2
  k21 <- Q / V3
  ka <- KA
  k_sum <- k10 + k12 + k21
  k_prod <- k10 * k21
  disc  <- sqrt(pmax(0, k_sum^2 - 4 * k_prod))
  alpha <- (k_sum + disc) / 2; beta <- (k_sum - disc) / 2
  
  # prevent division by 0 when ka approaches alpha or beta
  d_a  <- (ka - alpha); if (abs(d_a) < 1e-5) d_a <- 1e-5
  d_b  <- (ka - beta);  if (abs(d_b) < 1e-5) d_b <- 1e-5
  d_ab <- (beta - alpha); if (abs(d_ab) < 1e-5) d_ab <- -1e-5
  
  A <- (k21 - alpha) / (d_a * d_ab)
  B <- (k21 - beta)  / (d_b * (-d_ab))
  C <- (k21 - ka)    / ((-d_a) * (-d_b))
  
  t_eff <- pmax(0, xt - ALAG1)
  y <- (DOSE * ka / V2) * (A * exp(-alpha * t_eff) + B * exp(-beta * t_eff) + C * exp(-ka * t_eff))
  y[xt < ALAG1] <- 0
  
  list(y = matrix(pmax(as.numeric(y), 0), ncol = 1), poped.db = poped.db)
}

# Parameter model (exponential IIV)
sfg_model <- function(x, a, bpop, b, bocc) {
  c(CL = bpop[1] * exp(b[1]), 
    V2 = bpop[2] * exp(b[2]), 
    V3 = bpop[3] * exp(b[3]),
    Q  = bpop[4] * exp(b[4]), 
    KA = bpop[5] * exp(b[5]))
}

# Combined proportional and additive residual error model
feps_model <- function(model_switch, xt, parameters, epsi, poped.db) {
  returnArgs <- do.call(poped.db$model$ff_pointer, list(model_switch, xt, parameters, poped.db))
  f <- returnArgs[[1]]
  list(y = f + (f * epsi[, 1] + epsi[, 2]), poped.db = returnArgs[[2]])
}

# Reference NCA cache ====

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

cohort_cache <- list()
for (c_id in cohort_ids_to_run) {
  c_data   <- load_cohort_csv(c_id)
  nca_full <- run_ncappc(c_data %>% dplyr::select(ID, TIME, DV))
  
  # Fixed nominal population priors are shown as the reference configuration
  
  bpop <- c(CL = 0.190, V2 = 5.10, V3 = 32.5, Q = 0.336, KA = 1.201)
  d    <- c(0.04, 0.04, 0.04, 0.04, 0.09)
  
  cohort_cache[[as.character(c_id)]] <- list(
    data       = c_data,
    bpop       = bpop,
    d          = d,
    n_subjects = length(unique(c_data$ID)),
    true_cmax  = nca_full[, 6],
    true_auc   = nca_full[, 11],
    true_tmax  = nca_full[, 4]
  )
}

# PopED optimization loop across regimes ====

all_regimes_binary  <- list()
all_regimes_metrics <- list()

for (regime_name in names(regime_configs)) {
  target_n <- regime_configs[[regime_name]]
  message(sprintf("\n=== PopED Optimization: %s (Strict Total Samples = %d) ===", toupper(regime_name), target_n))
  
  # Select initial guess from candidate space
  init_idx <- unique(round(seq(1, length(discrete_space), length.out = target_n)))
  if (length(init_idx) < target_n) {
    init_idx <- sort(c(init_idx, setdiff(seq_along(discrete_space), init_idx)[1:(target_n - length(init_idx))]))
  }
  initial_xt <- discrete_space[init_idx]
  
  # When N = 4: fix Q (parameter 4) to ensure the 4x4 FIM has full rank
  if (target_n == 4L) {
    notfixed_bpop_setting <- c(1, 1, 1, 0, 1) # Estimates CL, V2, V3, KA; Q is fixed
  } else {
    notfixed_bpop_setting <- c(1, 1, 1, 1, 1) # Estimates all 5 parameters
  }
  
  poped_cohort_results <- foreach(
    c_id = cohort_ids_to_run,
    .packages = c("PopED", "dplyr"),
    .export = c("cohort_cache", "ff_2comp_oral", "sfg_model", "feps_model",
                "notfixed_bpop_setting", "initial_xt", "discrete_space", "target_n")
  ) %dopar% {
    c_entry <- cohort_cache[[as.character(c_id)]]
    
    poped_db <- PopED::create.poped.database(
      ff_fun         = ff_2comp_oral,
      fg_fun         = sfg_model,
      fError_fun     = feps_model,
      bpop           = c_entry$bpop,
      notfixed_bpop  = notfixed_bpop_setting,
      d              = c_entry$d,
      notfixed_d     = c(0, 0, 0, 0, 0),
      sigma          = c(0.05^2, 2^2),
      notfixed_sigma = c(0, 0),
      m              = 1,
      groupsize      = c_entry$n_subjects,
      xt             = matrix(initial_xt, nrow = 1),
      minxt          = matrix(min(discrete_space), nrow = 1, ncol = target_n),
      maxxt          = matrix(max(discrete_space), nrow = 1, ncol = target_n),
      discrete_xt    = list(discrete_space)
    )
    
    opt_res <- tryCatch({
      PopED::poped_optim(poped_db, opt_xt = TRUE, method = "LS", trace = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(opt_res)) {
      raw_opt <- opt_res$poped.db$design$xt[1, ]
    } else {
      raw_opt <- initial_xt
    }
    
    opt_times <- sort(unique(round(raw_opt, 6)))
    
    # Enforce exact target_n point count
    if (length(opt_times) > target_n) {
      opt_times <- opt_times[1:target_n]
    } else if (length(opt_times) < target_n) {
      missing_needed <- target_n - length(opt_times)
      fallback_pool  <- setdiff(discrete_space, opt_times)
      opt_times      <- sort(c(opt_times, fallback_pool[1:missing_needed]))
    }
    
    list(cohort_id = c_id, opt_times = opt_times)
  }
  
  regime_binary_list <- list()
  
  for (res in poped_cohort_results) {
    c_id         <- res$cohort_id
    final_design <- res$opt_times # Exactly 4, 7, or 15 points
    
    binary_row <- as.data.frame(t(as.integer(all_possible_times %in% final_design)))
    colnames(binary_row) <- paste0("T_", all_possible_times, "h")
    binary_row <- cbind(COHORT = c_id, DENSITY_REGIME = regime_name, N_SAMPLES = length(final_design), binary_row)
    
    regime_binary_list[[as.character(c_id)]] <- binary_row
    all_regimes_binary[[length(all_regimes_binary) + 1]] <- binary_row
    
    # Validation against ncappc
    c_entry    <- cohort_cache[[as.character(c_id)]]
    sparse_obs <- c_entry$data %>% dplyr::filter(round(TIME, 6) %in% round(final_design, 6))
    nca_cand   <- run_ncappc(sparse_obs %>% dplyr::select(ID, TIME, DV))
    val_met    <- calculate_standard_metrics(nca_cand[, 6], nca_cand[, 11], nca_cand[, 4],
                                             c_entry$true_cmax, c_entry$true_auc, c_entry$true_tmax)
    all_regimes_metrics[[length(all_regimes_metrics) + 1]] <- data.frame(DENSITY_REGIME = regime_name, COHORT = c_id, t(val_met))
  }
  
  write.csv(dplyr::bind_rows(regime_binary_list), file.path(results_dir, sprintf("best_schedules_poped_%s_binary.csv", regime_name)), row.names = FALSE)
  message(sprintf(">> Saved best_schedules_poped_%s_binary.csv (N = %d)", regime_name, target_n))
}

stopCluster(cl)

# Data exports and plots ====

master_binary_df  <- dplyr::bind_rows(all_regimes_binary)
master_metrics_df <- dplyr::bind_rows(all_regimes_metrics)

write.csv(master_binary_df,  file = file.path(results_dir, "all_best_schedules_poped_binary.csv"),   row.names = FALSE)
write.csv(master_metrics_df, file = file.path(results_dir, "selected_schedules_poped_nca_metrics.csv"), row.names = FALSE)

theme_set(theme_minimal(base_size = 12))

long_master <- master_binary_df %>%
  pivot_longer(cols = starts_with("T_"), names_to = "T_Label", values_to = "Selected") %>%
  mutate(Time = as.numeric(sub("^T_(.+)h$", "\\1", T_Label)),
         Status = factor(Selected, levels = c(0, 1), labels = c("Not Selected", "Selected")))

for (r_name in names(regime_configs)) {
  sub_binary <- long_master %>% dplyr::filter(DENSITY_REGIME == r_name)
  
  p_grid <- ggplot(sub_binary, aes(x = factor(Time), y = factor(COHORT), fill = Status)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Not Selected" = "#F0F2F5", "Selected" = "#6A1B9A")) +
    labs(title = sprintf("PopED D-Optimal Schedules (%s Regime: Strict N=%d)", toupper(r_name), regime_configs[[r_name]]),
         x = "Sampling Time (h)", y = "Cohort ID") +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
  
  ggsave(file.path(results_dir, sprintf("plot_final_schedule_only_poped_%s.png", r_name)), p_grid, width = 12, height = 5.5, dpi = 300)
}
message(">> PopED 96h pipeline completed successfully.")