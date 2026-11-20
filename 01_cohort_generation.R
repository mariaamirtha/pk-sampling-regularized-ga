# Title: Cohort Generation
# Purpose: Generate 1000 patient cohort then generate 10 subsets of 24 with replacement
# Author: Amirtha Maria

library(deSolve)
library(MASS)
library(dplyr)


# Simulate PK dataset and write directly to CSV ====

# set directory

data_dir <- "./data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

simulate_pk <- function(
    n_subjects = 24,
    pk_seed = 1568,
    error_seed = 98765,
    filename = "simulated_PK_data.csv"
) {
  
  # Typical PK parameters
  
  theta <- c(
    TVCL = 0.190,   # Clearance (L/h)
    TVV2 = 5.10,    # Central volume (L)
    TVV3 = 32.5,    # Peripheral volume (L)
    TVQ  = 0.336,   # Inter-compartmental clearance (L/h)
    TVKA = 1.201    # Absorption rate constant (1/h)
  )
  
  # Dose settings
  dose  <- 10000    # ug
  F1    <- 1
  ALAG1 <- 0.125    # h
  
  # Inter-individual variability
  omega <- c(
    OMEGA_CL = 0.04,
    OMEGA_V2 = 0.04,
    OMEGA_V3 = 0.04,
    OMEGA_Q  = 0.04,
    OMEGA_KA = 0.09
  )
  
  # Residual error
  SIGMA1_prop <- 0.05
  SIGMA2_add  <- 2
  
  # Lower limit of quantification
  LLOQ <- 75
  
  # Sampling times ====
  # Same grid as original GA simulation (Tsuchiwata and Tsuji 2023)
  
  times <- c(
    seq(0, 1, by = 0.125),
    seq(1.25, 2, by = 0.25),
    seq(2.5, 6, by = 0.5),
    seq(7, 12, by = 1),
    seq(14, 24, by = 2),
    seq(27, 48, by = 3),
    seq(54, 96, by = 6)
  )
  
  # Generate individual patient PK parameters ====

  
  set.seed(pk_seed)
  
  params <- MASS::mvrnorm(
    n = n_subjects,
    mu = log(theta),
    Sigma = diag(omega)
  )
  
  params <- exp(params)
  
  colnames(params) <- names(theta)
  
  # PK ODE model ====
  
  pkeq <- function(t, y, p) {
    
    CL <- p["TVCL"]
    V2 <- p["TVV2"]
    V3 <- p["TVV3"]
    Q  <- p["TVQ"]
    KA <- p["TVKA"]
    
    KE  <- CL / V2
    K23 <- Q / V2
    K32 <- Q / V3
    
    depot   <- y["depot"]
    central <- y["central"]
    periph  <- y["periph"]
    
    d_depot <- -KA * depot
    
    d_central <-
      KA * depot -
      (KE + K23) * central +
      K32 * periph
    
    d_periph <-
      K23 * central -
      K32 * periph
    
    list(
      c(
        depot   = d_depot,
        central = d_central,
        periph  = d_periph
      )
    )
  }
  
  # Dose event ====
  dose_events <- data.frame(
    var    = "depot",
    time   = ALAG1,
    value  = F1 * dose,
    method = "add"
  )
  
  initial_state <- c(
    depot   = 0,
    central = 0,
    periph  = 0
  )
  
  # Simulate each patient ====
  
  out_list <- lapply(
    seq_len(n_subjects),
    function(i) {
      
      out_i <- deSolve::ode(
        func = pkeq,
        times = times,
        y = initial_state,
        parms = params[i, ],
        events = list(data = dose_events)
      )
      
      out_i <- as.data.frame(out_i)
      
      # Noiseless concentration
      out_i$DV <- out_i$central / params[i, "TVV2"]
      
      # Subject ID
      out_i$ID <- i
      
      return(out_i)
    }
  )
  
  out <- dplyr::bind_rows(out_list)
  
  names(out) <- c(
    "TIME",
    "depot",
    "central",
    "periph",
    "DV",
    "ID"
  )
  
  # Save noiseless concentration ====
  out$ODV <- out$DV
  
  # Add residual error ====
  set.seed(error_seed)
  
  out$ERR1_prop <- rnorm(nrow(out))
  
  out$ERR2_add <- rnorm(nrow(out))
  
  out$DV <- ifelse(
    out$DV == 0,
    0,
    out$DV *
      (1 + out$ERR1_prop * SIGMA1_prop) +
      out$ERR2_add * SIGMA2_add
  )
  
  # BLQ handling ====
  # Same rule as original simulation (Tsuchiwata and Tsuji 2023)
  
  out <- out %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(TIME, .by_group = TRUE) %>%
    dplyr::mutate(
      
      # Mark concentration <= LLOQ as BLQ
      DV = ifelse(
        TIME != 0 & DV <= LLOQ,
        -99,
        DV
      ),
      
      # If previous observation was BLQ,
      # mark the next observation BLQ too
      DV = ifelse(
        TIME != 0 & dplyr::lag(DV) == -99,
        -99,
        DV
      ),
      
      # Replace BLQ marker with 0
      DV = ifelse(
        DV == -99,
        0,
        DV
      )
    ) %>%
    dplyr::ungroup()
  

  # Add patient-specific PK parameters ====
  
  param_df <- as.data.frame(params)
  
  param_df$ID <- seq_len(n_subjects)
  
  out <- dplyr::left_join(
    out,
    param_df,
    by = "ID"
  )
  
  # Final output dataset ====
  
  pk_output <- out %>%
    dplyr::select(
      ID,
      TIME,
      DV,
      ODV,
      TVCL,
      TVV2,
      TVV3,
      TVQ,
      TVKA
    ) %>%
    dplyr::arrange(
      ID,
      TIME
    )
  
  # Write CSV, print path, return dataset ====
  
  write.csv(
    pk_output,
    file = filename,
    row.names = FALSE
  )
  
  message(
    "Saved PK dataset to: ",
    normalizePath(filename, mustWork = FALSE)
  )
  
  return(pk_output)
}

# Optional quick sanity check (24 subjects) ====
# pk_data <- simulate_pk(
# n_subjects = 24,
#  pk_seed = 1568,
#  error_seed = 98765,
#  filename = file.path(data_dir,"simulated_PK_data.csv")
#)

# head(pk_data)

# dim(pk_data)

# file.exists("simulated_PK_data.csv")

# getwd()

# Generate 1,000 patient population ====

# Generate full cohort of 1000 subjects
pk_data_1000 <- simulate_pk(
  n_subjects = 1000,
  pk_seed    = 1568,
  error_seed = 98765,
  filename   = file.path(data_dir,"simulated_PK_data_1000.csv")
)

# Subset 10 cohort of 24 patients with replacement ====

sample_subcohorts <- function(
    full_data,
    n_cohorts      = 10,
    cohort_size    = 24,
    sampling_seed  = 2026,
    save_individual_csvs = FALSE
) {
  set.seed(sampling_seed)
  
  # Get all unique patient IDs from the base population
  all_ids <- unique(full_data$ID)
  
  # Sample 24 patient IDs for each of the 10 cohorts
  subcohort_list <- lapply(seq_len(n_cohorts), function(c_idx) {
    
    # Sample 24 patient IDs with replacement from the 1000-patient pool
    sampled_ids <- sample(all_ids, size = cohort_size, replace = TRUE)
    
    # Extract records for sampled patients while retaining duplicate draws
    cohort_df <- dplyr::bind_rows(
      lapply(seq_along(sampled_ids), function(rep_idx) {
        pid <- sampled_ids[rep_idx]
        
        full_data %>%
          dplyr::filter(ID == pid) %>%
          dplyr::mutate(
            ORIGINAL_ID = pid,       # Keep track of the source 1000-cohort ID
            NEW_ID      = rep_idx,   # Unique ID (1-24) within this sub-cohort
            COHORT      = c_idx      # Target sub-cohort index (1-10)
          )
      })
    )
    
    # Write individual sub-cohort CSV files
    if (save_individual_csvs) {
      write.csv(
        cohort_df,
        file = file.path(data_dir,sprintf("simulated_PK_cohort_%02d.csv", c_idx)),
        row.names = FALSE
      )
    }
    
    return(cohort_df)
  })
  
  # Combine all 10 sub-cohorts into a single master data frame
  combined_cohorts <- dplyr::bind_rows(subcohort_list)
  
  return(combined_cohorts)
}

# Run the sub-cohort sampler
subcohorts_df <- sample_subcohorts(
  full_data            = pk_data_1000,
  n_cohorts            = 10,
  cohort_size          = 24,
  sampling_seed        = 2026,
  save_individual_csvs = TRUE
)

# Verification checks ====

# Check total rows (Each cohort has 24 patients * time points per patient)
table(subcohorts_df$COHORT)

# View top of combined dataset
head(subcohorts_df)