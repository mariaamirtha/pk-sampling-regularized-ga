# Title: Manuscript and Appendix Figure Generator

# Purpose: This script ingests all the results from the data analysis script
# and generates all the figures in the manuscript and appendix with proper 
# scaling and IEEE formatted labels.

# Author: Amirtha Maria

# Input: 01_master_schedule_index.csv, 10_schedule_level_error_summary.csv,
# 21_primary_benchmark_tests.csv, 22_sensitivity1_factorial_swap_tests.csv,
# 23_sensitivity2_pchip15_vs_pchip49_tests.csv, 24_ablation_pure_ga_vs_tsuchiwata_tests.csv, 
# 25_friedman_method_omnibus_tests.csv, 26_friedman_sampling_intensity_omnibus_tests.csv

# Output: PDFs and PNGs of Figures 01-05 and Appendix Figures A1, B1-B5

# Clear workspace ====
rm(list = ls())
gc()

# Load packages ====
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "gridExtra", "pracma")
missing_pkgs  <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(pracma)
})

# Configure directories and ingest data ====

base_root <- "./results"
table_dir <- file.path(base_root, "tables")
plot_root <- file.path(base_root, "figures")

subdirs <- list(
  main     = file.path(plot_root, "01_main_text_figures"),
  appendix = file.path(plot_root, "02_appendix_figures")
)
for (d in subdirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

master_index     <- read.csv(file.path(table_dir, "01_master_schedule_index.csv"), stringsAsFactors = FALSE)
schedule_summary <- read.csv(file.path(table_dir, "10_schedule_level_error_summary.csv"), stringsAsFactors = FALSE)
primary_tests    <- read.csv(file.path(table_dir, "21_primary_benchmark_tests.csv"), stringsAsFactors = FALSE)
sens1_tests      <- read.csv(file.path(table_dir, "22_sensitivity1_factorial_swap_tests.csv"), stringsAsFactors = FALSE)
sens2_tests      <- read.csv(file.path(table_dir, "23_sensitivity2_pchip15_vs_pchip49_tests.csv"), stringsAsFactors = FALSE)
ablation_tests   <- read.csv(file.path(table_dir, "24_ablation_pure_ga_vs_tsuchiwata_tests.csv"), stringsAsFactors = FALSE)
friedman_method  <- read.csv(file.path(table_dir, "25_friedman_method_omnibus_tests.csv"), stringsAsFactors = FALSE)

rds_path       <- file.path(table_dir, "09_subject_errors_long.rds")
subject_errors <- if (file.exists(rds_path)) readRDS(rds_path) else NULL

# Establishing consistent color scheme for visualization, format for sub-labels ====
# and establishing file saving logistics

method_order <- c("Pure_GA", "PCHIP_49", "PopED", "Tsuchiwata_Paper_7pt", "Clinical", "PCHIP_15", "Tsuchiwata_GA")

palette_methods <- c(
  "Pure_GA"              = "#D55E00",
  "PCHIP_49"             = "#009E73",
  "PopED"                = "#0072B2",
  "Tsuchiwata_Paper_7pt" = "#56B4E9",
  "Clinical"             = "#222222",
  "PCHIP_15"             = "#6A3D9A",
  "Tsuchiwata_GA"        = "#CC79A7"
)

labels_methods <- c(
  "Pure_GA"              = "Pure GA",
  "PCHIP_49"             = "PCHIP-49 GA",
  "PopED"                = "PopED",
  "Tsuchiwata_Paper_7pt" = "Tsuchiwata 7pt",
  "Clinical"             = "Clinical 15pt",
  "PCHIP_15"             = "PCHIP-15 GA",
  "Tsuchiwata_GA"        = "Tsuchiwata GA"
)

metric_labels <- c(
  "Cmax"   = "Cmax",
  "AUC"    = "AUC (0-96h)",
  "Tmax"   = "Tmax",
  "CL_F"   = "CL/F",
  "Vz_F"   = "Vz/F",
  "T_half" = "t1/2"
)

regime_labels <- c(
  "sparse"      = "Sparse (N=4)",
  "semi-sparse" = "Semi-Sparse (N=7)",
  "intensive"   = "Intensive (N=15)"
)

format_p <- function(p) {
  if (!is.finite(p)) "NA" else if (p < 0.001) "<0.001" else formatC(p, digits = 3, format = "f")
}

# Main-text theme
plot_theme_ieee <- function(base_size = 11.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor      = element_blank(),
      strip.background      = element_rect(fill = "grey95", color = "grey80"),
      strip.text            = element_text(face = "bold", size = rel(0.95)),
      plot.title            = element_text(face = "bold", size = rel(1.08), color = "#0F172A"),
      plot.subtitle         = element_text(size = rel(0.88), color = "grey30", margin = margin(b = 6)),
      axis.title            = element_text(face = "bold", color = "#0F172A"),
      legend.position       = "bottom",
      legend.title          = element_text(face = "bold"),
      plot.caption          = element_text(hjust = 0.5, face = "bold", size = rel(1.08), color = "#0F172A", margin = margin(t = 8, b = 4)),
      plot.caption.position = "panel"
    )
}

# Appendix theme
plot_theme_app_ieee <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor      = element_blank(),
      strip.background      = element_rect(fill = "grey95", color = "grey80"),
      strip.text            = element_text(face = "bold", size = rel(0.95)),
      plot.title            = element_blank(),
      plot.subtitle         = element_blank(),
      axis.title            = element_text(face = "bold", color = "#0F172A"),
      legend.position       = "bottom",
      legend.title          = element_text(face = "bold"),
      plot.caption          = element_text(hjust = 0.5, face = "bold", size = rel(1.08), color = "#0F172A", margin = margin(t = 8, b = 4)),
      plot.caption.position = "panel"
    )
}

save_pub_plot <- function(p, dir_path, filename, w = 8.8, h = 10.4) {
  png_dir <- file.path(dir_path, "png")
  pdf_dir <- file.path(dir_path, "pdf")
  if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)
  if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
  
  ggsave(file.path(png_dir, paste0(filename, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(pdf_dir, paste0(filename, ".pdf")), p, width = w, height = h, bg = "white")
  message(">> Saved: ", filename, " (PNG + PDF)")
}

# Establish primary data set
# Preserves PopED Sparse row so x-axis slots stay aligned
df_primary_all <- schedule_summary %>%
  filter(
    Analysis_Type %in% c("Primary", "Benchmark"),
    Evaluation_Engine == "Direct_NCA",
    Schedule_Family %in% c("Pure_GA", "PCHIP_49", "PopED", "Tsuchiwata_Paper_7pt", "Clinical")
  ) %>%
  mutate(
    Regime_Label = factor(Sampling_Threshold, levels = names(regime_labels), labels = regime_labels),
    Family_Label = factor(Schedule_Family, levels = method_order, labels = labels_methods[method_order]),
    Is_Static    = Schedule_Family %in% c("Clinical", "Tsuchiwata_Paper_7pt")
  )

df_primary <- df_primary_all %>% filter(is.finite(MAPE))

# Figure 1: Kinetic anchor distribution semi-sparse (N=7) and sparse (N=4) comparison ====

message(">> Generating IEEE Figure 1...")

sched_idx <- read.csv(file.path(table_dir, "01_master_schedule_index.csv"), stringsAsFactors = FALSE)

parse_times <- function(s) as.numeric(strsplit(s, ";")[[1]])

# Extract modal schedule across the 10 primary Direct_NCA cohorts
get_modal_schedule <- function(df, family, threshold) {
  sub <- df[df$Schedule_Family == family &
              df$Sampling_Threshold == threshold &
              df$Evaluation_Engine == "Direct_NCA" &
              df$Analysis_Type == "Primary", ]
  tbl <- sort(table(sub$Schedule_Times), decreasing = TRUE)
  parse_times(names(tbl)[1])
}

# for families with 10 unique cohort schedules, pick the balanced cohort
# with the highest cross-cohort timepoint frequency overlap

get_consensus_rep_schedule <- function(df, family, threshold) {
  sub <- df[df$Schedule_Family == family &
              df$Sampling_Threshold == threshold &
              df$Evaluation_Engine == "Direct_NCA" &
              df$Analysis_Type == "Primary", ]
  all_t <- unlist(lapply(sub$Schedule_Times, parse_times))
  freq  <- table(all_t)
  scores <- sapply(sub$Schedule_Times, function(s) {
    ts <- parse_times(s)
    sum(freq[as.character(ts)])
  })
  
  # Tie-break by earliest Cohort_Match_ID
  best_idx <- which.max(scores)
  list(
    cohort = sub$Cohort_Match_ID[best_idx],
    times  = sort(parse_times(sub$Schedule_Times[best_idx]))
  )
}

# Structural two-compartment oral PK curve
dose_ug <- 10000; CL <- 0.190; V2 <- 5.10; V3 <- 32.5; Q <- 0.336; ka <- 1.201; ALAG1 <- 0.125
k10 <- CL / V2; k12 <- Q / V2; k21 <- Q / V3
b0  <- k10 + k12 + k21; b1 <- k10 * k21
disc <- sqrt(pmax(0, b0^2 - 4 * b1))
alpha <- (b0 + disc) / 2; beta <- (b0 - disc) / 2

A_coef <- (dose_ug * ka / V2) * (k21 - alpha) / ((ka - alpha) * (beta - alpha))
B_coef <- (dose_ug * ka / V2) * (k21 - beta)  / ((ka - beta)  * (alpha - beta))
C_coef <- (dose_ug * ka / V2) * (k21 - ka)    / ((alpha - ka) * (beta - ka))

pk_curve_table1 <- function(t) {
  t_eff <- pmax(0, t - ALAG1)
  y <- A_coef * exp(-alpha * t_eff) + B_coef * exp(-beta * t_eff) + C_coef * exp(-ka * t_eff)
  y[t < ALAG1] <- 0
  pmax(0, y)
}

t_dense        <- seq(0, 96, length.out = 1500)
df_curve       <- data.frame(Time = t_dense, Conc = pk_curve_table1(t_dense))
candidate_grid <- c(seq(0, 1, 0.125), seq(1.25, 2, 0.25), seq(2.5, 6, 0.5),
                    seq(7, 12, 1), seq(14, 24, 2), seq(27, 48, 3), seq(54, 96, 6))
df_candidates  <- data.frame(Time = candidate_grid, Conc = pk_curve_table1(candidate_grid))
y_max_lim      <- 1950

# Representative Semi-Sparse (N = 7) Pure GA Schedule

ga7_rep    <- get_consensus_rep_schedule(sched_idx, "Pure_GA", "semi-sparse")
ga_7_times <- ga7_rep$times  # c(2, 2.5, 3, 9, 22, 39, 84)

df_ga7 <- data.frame(
  Time    = ga_7_times,
  Conc    = pk_curve_table1(ga_7_times),
  Label   = paste0(ga_7_times, "h"),
  Label_X = c(2.0,  11.5, 12.0, 16.0, 24.0, 39.0, 84.0),
  Label_Y = c(1830, 1760, 1490, 1040, 540,  380,  320)
)

phases_ga <- data.frame(
  xmin = c(0, 3, 24), xmax = c(3, 24, 96), ymin = -Inf, ymax = Inf,
  Phase = factor(
    c("1. Absorption & Peak (0-3h)", "2. Distribution (3-24h)", "3. Terminal Elimination (24-96h)"),
    levels = c("1. Absorption & Peak (0-3h)", "2. Distribution (3-24h)", "3. Terminal Elimination (24-96h)")
  )
)

# (a) Semi-sparse Pure GA sampling point distribution

p_fig1a <- ggplot() +
  geom_rect(data = phases_ga, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Phase), alpha = 0.11) +
  scale_fill_manual(values = c("#0284C7", "#6366F1", "#0F172A"), name = "Kinetic Phase") +
  geom_point(data = df_candidates, aes(x = Time, y = 0), shape = 3, size = 1.5, color = "grey50") +
  geom_line(data = df_curve, aes(x = Time, y = Conc), color = "#0F172A", linewidth = 1.15) +
  geom_segment(data = df_ga7, aes(x = Time, xend = Time, y = 0, yend = Conc),
               linetype = "dashed", color = "#D55E00", linewidth = 0.7) +
  geom_segment(data = df_ga7, aes(x = Time, y = Conc, xend = Label_X, yend = Label_Y),
               color = "#D55E00", linewidth = 0.55, alpha = 0.85) +
  geom_point(data = df_ga7, aes(x = Time, y = Conc), shape = 21, fill = "#D55E00", color = "white", size = 3.8, stroke = 1.2) +
  geom_label(data = df_ga7, aes(x = Label_X, y = Label_Y, label = Label),
             size = 3.1, fontface = "bold", color = "#D55E00", fill = "white",
             label.size = 0.25, label.padding = unit(0.16, "lines")) +
  annotate("rect", xmin = 42, xmax = 94, ymin = 880, ymax = 1830,
           fill = "white", color = "grey70", linewidth = 0.65, alpha = 0.96) +
  annotate("text", x = 44, y = 1690, hjust = 0, size = 3.4, fontface = "bold",
           label = "Combinatorial Search Space (|T| = 49):") +
  annotate("text", x = 44, y = 1470, hjust = 0, size = 3.1, color = "grey25",
           label = "- Sparse (N = 4):       211,876 designs") +
  annotate("text", x = 44, y = 1250, hjust = 0, size = 3.1, fontface = "bold", color = "#D55E00",
           label = "- Semi-Sparse (N = 7):  85,900,584 designs") +
  annotate("text", x = 44, y = 1030, hjust = 0, size = 3.1, color = "grey25",
           label = "- Intensive (N = 15):    1.25 x 10^12 designs") +
  scale_x_continuous(breaks = c(0, 4, 12, 24, 48, 72, 96), limits = c(-2, 98), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, y_max_lim), expand = FALSE) +
  labs(title    = "Semi-Sparse Pure GA Kinetic Anchor Distribution (N = 7, Consensus Cohort 01)",
       subtitle = "Brackets peak & spans all 3 phases (family mean: 96.8% CL/F estimability; 1.85% Cmax MAPE)",
       x        = "Post-Dose Sampling Time (Hours)",
       y        = "Plasma Concentration (ng/mL)",
       caption  = "(a)") +
  plot_theme_ieee(base_size = 11.5) +
  theme(plot.margin = margin(t = 8, r = 16, b = 8, l = 8))

# (b) Modal Sparse (N = 4) PopED vs. Pure GA Comparison
poped_4_times <- get_modal_schedule(sched_idx, "PopED", "sparse")    # c(0.25, 4, 18, 48)
ga_4_times    <- get_modal_schedule(sched_idx, "Pure_GA", "sparse")  # c(2, 12, 27, 72)

poped_lbl <- "PopED D-Optimal Modal (N = 4: 0.0% CL/F Estimability)"
ga4_lbl   <- "Pure GA Modal Design (N = 4: 83.4% Modal / 68.2% Mean CL/F Estimability)"

df_poped4 <- data.frame(
  Time    = poped_4_times,
  Conc    = pk_curve_table1(poped_4_times),
  Method  = poped_lbl,
  Label   = paste0(poped_4_times, "h"),
  Label_X = c(6.5,  11.5, 21.5, 48.0),
  Label_Y = c(390,  1520, 630,  360)
)

df_ga4 <- data.frame(
  Time    = ga_4_times,
  Conc    = pk_curve_table1(ga_4_times),
  Method  = ga4_lbl,
  Label   = paste0(ga_4_times, "h"),
  Label_X = c(4.5,  15.0, 31.5, 72.0),
  Label_Y = c(1810, 930,  440,  330)
)

df_n4_compare <- bind_rows(df_poped4, df_ga4)

phases_poped <- data.frame(
  xmin = c(0, 24), xmax = c(24, 96), ymin = -Inf, ymax = Inf,
  Zone = factor(
    c("Absorption & Distribution (0-24h)", "Terminal Elimination Tail (24-96h)"),
    levels = c("Absorption & Distribution (0-24h)", "Terminal Elimination Tail (24-96h)")
  )
)

p_fig1b <- ggplot() +
  geom_rect(data = phases_poped, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Zone), alpha = 0.11) +
  scale_fill_manual(values = c("#0072B2", "#DC2626"), name = "Kinetic Domain") +
  geom_line(data = df_curve, aes(x = Time, y = Conc), color = "#0F172A", linewidth = 1.15) +
  geom_segment(data = df_n4_compare, aes(x = Time, xend = Time, y = 0, yend = Conc, color = Method, linetype = Method),
               linewidth = 0.7) +
  geom_segment(data = df_n4_compare, aes(x = Time, y = Conc, xend = Label_X, yend = Label_Y, color = Method),
               linewidth = 0.55, alpha = 0.85, show.legend = FALSE) +
  geom_point(data = df_n4_compare, aes(x = Time, y = Conc, color = Method, shape = Method),
             size = 3.8, stroke = 1.2, fill = "white") +
  geom_label(data = df_n4_compare, aes(x = Label_X, y = Label_Y, label = Label, color = Method),
             size = 3.0, fontface = "bold", fill = "white", label.size = 0.25,
             label.padding = unit(0.15, "lines"), show.legend = FALSE) +
  scale_color_manual(values = setNames(c("#0072B2", "#D55E00"), c(poped_lbl, ga4_lbl)),
                     name = "Sparse (N = 4) Schedule") +
  scale_shape_manual(values = setNames(c(21, 23), c(poped_lbl, ga4_lbl)),
                     name = "Sparse (N = 4) Schedule") +
  scale_linetype_manual(values = setNames(c("dashed", "solid"), c(poped_lbl, ga4_lbl)),
                        name = "Sparse (N = 4) Schedule") +
  guides(
    fill     = guide_legend(order = 1, nrow = 1),
    color    = guide_legend(order = 2, ncol = 1),
    shape    = guide_legend(order = 2, ncol = 1),
    linetype = guide_legend(order = 2, ncol = 1)
  ) +
  annotate("rect", xmin = 32, xmax = 95, ymin = 820, ymax = 1830,
           fill = "#FEF2F2", color = "#DC2626", linewidth = 0.65, alpha = 0.96) +
  annotate("text", x = 34, y = 1680, hjust = 0, size = 3.3, fontface = "bold", color = "#991B1B",
           label = "Sparse (N = 4) Head-to-Head Mechanism:") +
  annotate("text", x = 34, y = 1440, hjust = 0, size = 2.95, color = "#0072B2",
           label = "- PopED (0.25, 4, 18, 48h): 0.25h flank draw -> peak at 4h -> only 2 post-peak draws") +
  annotate("text", x = 34, y = 1210, hjust = 0, size = 2.95, color = "#D55E00",
           label = "- Pure GA (2, 12, 27, 72h): Skips flank -> 1 peak + 3 post-peak draws -> 83.4% est.") +
  annotate("text", x = 34, y = 980, hjust = 0, size = 2.9, fontface = "italic", color = "#7F1D1D",
           label = "- Regulatory NCA requires >= 3 post-peak points for log-linear lambda_z regression") +
  scale_x_continuous(breaks = c(0, 4, 12, 24, 48, 72, 96), limits = c(-2, 98), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, y_max_lim), expand = FALSE) +
  labs(title    = "Sparse (N = 4) Comparison: PopED Pre-Peak Allocation vs. Pure GA Post-Peak Anchoring",
       subtitle = "At identical sparse budgets (N = 4), PopED's 0.25h ka-sensitivity draw leaves only 2 post-peak points",
       x        = "Post-Dose Sampling Time (Hours)",
       y        = "Plasma Concentration (ng/mL)",
       caption  = "(b)") +
  plot_theme_ieee(base_size = 11.5) +
  theme(
    legend.box   = "vertical",
    legend.text  = element_text(size = 9.5),
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin  = margin(t = 8, r = 16, b = 10, l = 8)
  )

fig1_vertical <- arrangeGrob(p_fig1a, p_fig1b, nrow = 2, ncol = 1, heights = c(1.0, 1.08))
save_pub_plot(fig1_vertical, subdirs$main, "Figure_01_IEEE_Kinetic_Anchors_and_N4_Comparison", w = 8.8, h = 11.5)

# Figure 2: 6 Parameter primary benchmark and Sampling burden vs estimation error ====

message(">> Generating IEEE Figure 2...")

# (a) Out-of-Sample Exposure Recovery

df_fig2a_multi  <- df_primary %>% filter(Metric %in% c("AUC", "Cmax"), !Is_Static)
df_fig2a_static <- df_primary %>% filter(Metric %in% c("AUC", "Cmax"), Is_Static)

p_fig2a <- ggplot() +
  geom_boxplot(data = df_fig2a_multi, aes(x = Family_Label, y = MAPE, fill = Schedule_Family),
               alpha = 0.45, outlier.shape = NA, width = 0.55) +
  geom_jitter(data = df_fig2a_multi, aes(x = Family_Label, y = MAPE, color = Schedule_Family),
              width = 0.12, alpha = 0.75, size = 2.0) +
  stat_summary(data = df_fig2a_multi, aes(x = Family_Label, y = MAPE),
               fun = mean, geom = "point", shape = 23, size = 3.2, fill = "white", stroke = 1.0) +
  geom_point(data = df_fig2a_static, aes(x = Family_Label, y = MAPE, fill = Schedule_Family, color = Schedule_Family),
             shape = 23, size = 3.8, stroke = 1.2) +
  facet_grid(Metric ~ Regime_Label, scales = "free", labeller = labeller(Metric = metric_labels)) +
  scale_x_discrete(drop = TRUE, expand = expansion(add = c(0.55, 0.55))) +
  scale_fill_manual(values = palette_methods) +
  scale_color_manual(values = palette_methods) +
  labs(title    = "Out-of-Sample Exposure Recovery Across Sparsity Tiers (N = 1,000 Subjects)",
       subtitle = "White diamonds = 10-cohort family means; colored diamonds = static literature benchmarks",
       x        = NULL,
       y        = "Schedule-Level MAPE (%)",
       caption  = "(a)") +
  plot_theme_ieee(base_size = 11.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

# (b) Sampling Burden vs. Estimation Error

burden_df <- schedule_summary %>%
  filter(
    Analysis_Type == "Primary",
    Evaluation_Engine == "Direct_NCA",
    Schedule_Family %in% c("Pure_GA", "PopED", "PCHIP_49"),
    Metric %in% c("AUC", "Cmax"),
    is.finite(MAPE)
  ) %>%
  group_by(Schedule_Family, Number_of_Nominal_Samples, Metric) %>%
  summarise(Mean_MAPE = mean(MAPE, na.rm = TRUE), SD_MAPE = sd(MAPE, na.rm = TRUE), .groups = "drop")

p_fig2b <- ggplot(burden_df, aes(x = Number_of_Nominal_Samples, y = Mean_MAPE, color = Schedule_Family, group = Schedule_Family)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.3) +
  geom_errorbar(aes(ymin = pmax(0, Mean_MAPE - SD_MAPE), ymax = Mean_MAPE + SD_MAPE), width = 0.35, linewidth = 0.6) +
  facet_wrap(~ Metric, nrow = 1, scales = "free_y", labeller = labeller(Metric = metric_labels)) +
  scale_x_continuous(breaks = c(4, 7, 15), limits = c(3.5, 15.5)) +
  scale_color_manual(values = palette_methods, labels = labels_methods) +
  labs(title    = "Sampling Burden vs. Estimation Error: Diminishing Returns Past N = 7",
       subtitle = "Semi-sparse GA schedules (N = 7) achieve near-intensive exposure accuracy at <47% blood volume",
       x        = "Total Blood Draws per Subject (N)",
       y        = "Mean Schedule-Level MAPE (%)",
       color    = "Optimization Method",
       caption  = "(b)") +
  plot_theme_ieee(base_size = 11.5)

fig2_vertical <- arrangeGrob(p_fig2a, p_fig2b, nrow = 2, ncol = 1, heights = c(1.35, 1.0))
save_pub_plot(fig2_vertical, subdirs$main, "Figure_02_IEEE_Primary_Benchmark_and_Burden", w = 8.8, h = 10.8)

# Figure 3: Sensitivity analysis 1 (recovery engine comparison and mechanism visualization) ====

message(">> Generating IEEE Figure 3...")

# (a) t1/2 Error Inflation

df_fig3a <- schedule_summary %>%
  filter(
    Analysis_Type %in% c("Primary", "Sensitivity_1"),
    Sampling_Threshold == "semi-sparse",
    Schedule_Family %in% c("Pure_GA", "PCHIP_49"),
    Metric == "T_half",
    is.finite(MAPE)
  ) %>%
  mutate(
    Pairing = factor(
      case_when(
        Schedule_Family == "Pure_GA"  & Evaluation_Engine == "Direct_NCA"          ~ "Pure GA\n+ Direct NCA",
        Schedule_Family == "PCHIP_49" & Evaluation_Engine == "Direct_NCA"          ~ "PCHIP-49\n+ Direct NCA",
        Schedule_Family == "PCHIP_49" & Evaluation_Engine == "PCHIP_Reconstruction" ~ "PCHIP-49\n+ PCHIP Spline",
        Schedule_Family == "Pure_GA"  & Evaluation_Engine == "PCHIP_Reconstruction" ~ "Pure GA\n+ PCHIP Spline"
      ),
      levels = c("Pure GA\n+ Direct NCA", "PCHIP-49\n+ Direct NCA", "PCHIP-49\n+ PCHIP Spline", "Pure GA\n+ PCHIP Spline")
    ),
    Engine_Label = factor(ifelse(Evaluation_Engine == "Direct_NCA", "Direct Trapezoidal NCA", "PCHIP Spline Reconstruction"),
                          levels = c("Direct Trapezoidal NCA", "PCHIP Spline Reconstruction"))
  )

df_fig3a_sum <- df_fig3a %>%
  group_by(Pairing, Engine_Label) %>%
  summarise(Mean_MAPE = mean(MAPE), SE_MAPE = sd(MAPE) / sqrt(n()), .groups = "drop")

p_fig3a <- ggplot(df_fig3a_sum, aes(x = Pairing, y = Mean_MAPE, fill = Engine_Label)) +
  geom_col(width = 0.54, color = "#0F172A", linewidth = 0.5, alpha = 0.9) +
  geom_errorbar(aes(ymin = pmax(0, Mean_MAPE - SE_MAPE), ymax = Mean_MAPE + SE_MAPE),
                width = 0.16, linewidth = 0.65, color = "#0F172A") +
  geom_text(aes(label = sprintf("%.1f%%", Mean_MAPE), y = Mean_MAPE + SE_MAPE + 13),
            fontface = "bold", size = 3.6, color = "#0F172A") +
  annotate("rect", xmin = 2.55, xmax = 4.45, ymin = 88, ymax = 352,
           fill = NA, color = "#B91C1C", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 3.5, y = 338, size = 3.5, fontface = "bold", color = "#B91C1C",
           label = "Spline Artifact: +174.7% Error Surge (p_adj = 0.0117*)") +
  scale_fill_manual(values = c("Direct Trapezoidal NCA" = "#0072B2", "PCHIP Spline Reconstruction" = "#DC2626"),
                    name = "Post-Hoc Evaluation Engine") +
  coord_cartesian(ylim = c(0, 365)) +
  labs(title    = "2x2 Factorial Engine Swap: Terminal Half-Life (t1/2) Error (N = 7)",
       subtitle = "PCHIP spline reconstruction inflates t1/2 error (Mean +/- SE) regardless of discovery engine",
       x        = "Optimization & Post-Hoc Evaluation Pairing",
       y        = "Mean Absolute Percentage Error (MAPE %)",
       caption  = "(a)") +
  plot_theme_ieee(base_size = 11.5)

# (b) Spline Tail Distortion 
t_tail_knots <- c(8, 24, 48, 96)
y_tail_knots <- pk_curve_table1(t_tail_knots)
t_tail_grid  <- seq(48, 96, length.out = 300)
y_true_tail  <- pk_curve_table1(t_tail_grid)
y_pchip_lin  <- pracma::pchip(t_tail_knots, y_tail_knots, t_tail_grid)

# Compute exact 3-point terminal log-linear t1/2 on the 49-pt grid (84, 90, 96h)
t_reg       <- c(84, 90, 96)
y_reg_true  <- pk_curve_table1(t_reg)
y_reg_pchip <- pracma::pchip(t_tail_knots, y_tail_knots, t_reg)

t12_true  <- log(2) / (-coef(lm(log(y_reg_true)  ~ t_reg))[2])
t12_pchip <- log(2) / (-coef(lm(log(y_reg_pchip) ~ t_reg))[2])
t12_infl  <- (t12_pchip - t12_true) / t12_true * 100

df_tail_mech <- bind_rows(
  data.frame(Time = t_tail_grid, Conc = y_true_tail,
             Mode = "True 2-Compartment Terminal Decay (beta-Phase)"),
  data.frame(Time = t_tail_grid, Conc = y_pchip_lin,
             Mode = "PCHIP Cubic Polynomial in Linear Space")
)

# Reconstructed 49-pt grid points between 48h and 96h (6h steps) fed into ncappc
t_recon_pts <- seq(54, 90, by = 6)
df_recon_pts <- data.frame(
  Time = t_recon_pts,
  Conc = pracma::pchip(t_tail_knots, y_tail_knots, t_recon_pts)
)

df_knots_plot <- data.frame(
  Time    = c(48, 96),
  Conc    = pk_curve_table1(c(48, 96)),
  Label   = c("48h (Knot)", "96h (Knot: dC/dt = 0)"),
  Label_X = c(53.0, 88.0),
  Label_Y = c(121.0, 86.8)
)

p_fig3b <- ggplot() +
  # Highlight trailing 3-point ncappc log-linear regression window (84-96h)
  annotate("rect", xmin = 84, xmax = 96, ymin = 84, ymax = 125,
           fill = "#DC2626", alpha = 0.08) +
  annotate("text", x = 90, y = 105.5, size = 2.75, fontface = "italic", color = "#991B1B",
           label = "Trailing 3-pt\nlambda_z window") +
  geom_line(data = df_tail_mech, aes(x = Time, y = Conc, color = Mode, linetype = Mode),
            linewidth = 1.2) +
  # Show interpolated 49-pt grid points along the bowed PCHIP curve
  geom_point(data = df_recon_pts, aes(x = Time, y = Conc),
             shape = 21, fill = "white", color = "#DC2626", size = 2.4, stroke = 0.9) +
  geom_segment(data = df_knots_plot, aes(x = Time, y = Conc, xend = Label_X, yend = Label_Y),
               color = "#D55E00", linewidth = 0.55) +
  geom_point(data = df_knots_plot, aes(x = Time, y = Conc),
             shape = 21, fill = "#D55E00", color = "white", size = 4.0, stroke = 1.2) +
  geom_label(data = df_knots_plot, aes(x = Label_X, y = Label_Y, label = Label),
             fontface = "bold", size = 3.0, color = "#D55E00", fill = "white",
             label.size = 0.25, label.padding = unit(0.15, "lines")) +
  scale_color_manual(
    values = c("True 2-Compartment Terminal Decay (beta-Phase)" = "#0F172A",
               "PCHIP Cubic Polynomial in Linear Space"          = "#DC2626"),
    name = "Elimination Trajectory"
  ) +
  scale_linetype_manual(
    values = c("True 2-Compartment Terminal Decay (beta-Phase)" = "solid",
               "PCHIP Cubic Polynomial in Linear Space"          = "dashed"),
    name = "Elimination Trajectory"
  ) +
  scale_x_continuous(breaks = seq(48, 96, by = 6)) +
  scale_y_log10(breaks = c(85, 90, 95, 100, 110, 120)) +
  coord_cartesian(xlim = c(46, 98), ylim = c(84, 125), expand = FALSE) +
  annotate("rect", xmin = 59, xmax = 96.5, ymin = 110.5, ymax = 123.5,
           fill = "white", color = "grey70", linewidth = 0.6, alpha = 0.96) +
  annotate("text", x = 60.2, y = 121.0, hjust = 0, size = 3.05, fontface = "bold",
           label = "Mechanism of lambda_z Slope Corruption (84-96h):") +
  annotate("text", x = 60.2, y = 117.3, hjust = 0, size = 2.85, color = "grey25",
           label = sprintf("- True beta-phase tail is log-linear -> t1/2 = %.0f h", t12_true)) +
  annotate("text", x = 60.2, y = 113.7, hjust = 0, size = 2.85, color = "#B91C1C",
           label = sprintf("- PCHIP bows & clamps slope to 0 at 96h -> t1/2 = %.0f h (%+.0f%%)", t12_pchip, t12_infl)) +
  labs(title    = "Pharmacometric Mechanism in the Terminal Elimination Tail (Semi-Log Scale, 48-96h)",
       subtitle = "Linear-space PCHIP introduces cubic bowing and flattens into the final knot (open circles = 49-pt grid)",
       x        = "Post-Dose Sampling Time (Hours)",
       y        = "Plasma Concentration (ng/mL, Log10 Scale)",
       caption  = "(b)") +
  plot_theme_ieee(base_size = 11.5)

fig3_vertical <- arrangeGrob(p_fig3a, p_fig3b, nrow = 2, ncol = 1)
save_pub_plot(fig3_vertical, subdirs$main, "Figure_03_IEEE_Sensitivity1_Factorial_and_Mechanism", w = 8.8, h = 10.8)

# Figure 4: Sensitivity analysis 2 (PCHIP-49 VS PCHIP-15) ====

message(">> Generating IEEE Figure 4...")

# (a) PCHIP-49 vs. PCHIP-15 Parameter Invariance

df_sens2 <- schedule_summary %>%
  filter(
    Analysis_Type %in% c("Primary", "Sensitivity_2"),
    Evaluation_Engine == "Direct_NCA",
    Schedule_Family %in% c("PCHIP_49", "PCHIP_15"),
    is.finite(MAPE)
  ) %>%
  group_by(Schedule_Family, Sampling_Threshold, Metric) %>%
  summarise(Mean_MAPE = mean(MAPE), SE_MAPE = sd(MAPE) / sqrt(n()), .groups = "drop") %>%
  mutate(
    Family_Label = factor(Schedule_Family, levels = c("PCHIP_49", "PCHIP_15"),
                          labels = c("PCHIP-49 (Dense 49-pt Grid)", "PCHIP-15 (Clinical 15-pt Grid)")),
    Metric_Label = factor(Metric, levels = names(metric_labels), labels = metric_labels),
    Regime_Short = factor(Sampling_Threshold, levels = c("sparse", "semi-sparse", "intensive"),
                          labels = c("Sparse (N = 4)", "Semi-Sparse (N = 7)", "Intensive (N = 15)"))
  )

df_fig4a <- df_sens2 %>% filter(Sampling_Threshold == "semi-sparse")

p_fig4a <- ggplot(df_fig4a, aes(x = Metric_Label, y = Mean_MAPE, fill = Family_Label)) +
  geom_col(position = position_dodge(0.72), width = 0.62, color = "#0F172A", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, Mean_MAPE - SE_MAPE), ymax = Mean_MAPE + SE_MAPE),
                position = position_dodge(0.72), width = 0.18, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Mean_MAPE), y = Mean_MAPE + SE_MAPE + 2.8),
            position = position_dodge(0.72), size = 3.1, fontface = "bold") +
  annotate("rect", xmin = 0.55, xmax = 6.45, ymin = 76, ymax = 88,
           fill = "white", color = "grey50", linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 3.5, y = 82, size = 3.4, fontface = "bold", color = "#0F172A",
           label = "All 6 Endpoints Statistically Invariant at N = 7 (Holm-adjusted p_adj >= 0.703)") +
  scale_fill_manual(values = c("PCHIP-49 (Dense 49-pt Grid)" = "#009E73",
                               "PCHIP-15 (Clinical 15-pt Grid)" = "#6A3D9A"), name = "Surrogate Target Grid") +
  coord_cartesian(ylim = c(0, 92)) +
  labs(title    = "Semi-Sparse Clinical Target (N = 7): Complete Parameter Invariance",
       subtitle = "Restricting in-loop PCHIP reconstruction to the 15-point clinical grid does not degrade semi-sparse recovery",
       x        = "Pharmacokinetic Parameter",
       y        = "Mean Absolute Percentage Error (MAPE %)",
       caption  = "(a)") +
  plot_theme_ieee(base_size = 11.5)

# (b) PCHIP-49 vs. PCHIP-15 Peak Exposure Dynamics

df_fig4b <- df_sens2 %>% filter(Metric %in% c("Cmax", "Tmax"))

p_fig4b <- ggplot(df_fig4b, aes(x = Regime_Short, y = Mean_MAPE, fill = Family_Label)) +
  geom_col(position = position_dodge(0.72), width = 0.58, color = "#0F172A", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, Mean_MAPE - SE_MAPE), ymax = Mean_MAPE + SE_MAPE),
                position = position_dodge(0.72), width = 0.18, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Mean_MAPE), y = Mean_MAPE + SE_MAPE + (ifelse(Metric == "Cmax", 0.25, 1.3))),
            position = position_dodge(0.72), size = 3.2, fontface = "bold") +
  facet_wrap(~ Metric_Label, scales = "free_y") +
  scale_fill_manual(values = c("PCHIP-49 (Dense 49-pt Grid)" = "#009E73",
                               "PCHIP-15 (Clinical 15-pt Grid)" = "#6A3D9A"), name = "Surrogate Target Grid") +
  labs(title    = "Peak Exposure Dynamics (Cmax & Tmax) Across Sampling Regimes",
       subtitle = "Dense 49-pt target grid improves peak capture exclusively at rich sampling (N = 15, p_adj = 0.0117*)",
       x        = "Sampling Density Regime",
       y        = "Mean Absolute Percentage Error (MAPE %)",
       caption  = "(b)") +
  plot_theme_ieee(base_size = 11.5)

fig4_vertical <- arrangeGrob(p_fig4a, p_fig4b, nrow = 2, ncol = 1)
save_pub_plot(fig4_vertical, subdirs$main, "Figure_04_IEEE_Sensitivity2_Grid_Resolution", w = 8.8, h = 10.6)

# Objective function benchmark (N=7) ====

message(">> Generating IEEE Figure 5...")

# (a) Clearance Estimability Failure Semi-Sparse

df_fig5a <- schedule_summary %>%
  filter(
    Sampling_Threshold == "semi-sparse",
    Evaluation_Engine == "Direct_NCA",
    Schedule_Family %in% c("Pure_GA", "Tsuchiwata_GA", "Tsuchiwata_Paper_7pt"),
    Metric == "CL_F"
  ) %>%
  group_by(Schedule_Family) %>%
  summarise(Fail_Rate = mean(Nonestimable_pct), SE_Fail = sd(Nonestimable_pct) / sqrt(n()), .groups = "drop") %>%
  mutate(
    Arch_Label = factor(Schedule_Family,
                        levels = c("Pure_GA", "Tsuchiwata_GA", "Tsuchiwata_Paper_7pt"),
                        labels = c("Pure GA (Ours: Quadratic Barrier)", "Tsuchiwata GA (Step-Cliff Penalty)", "Tsuchiwata et al.\n(Published 7pt)"))
  )

p_fig5a <- ggplot(df_fig5a, aes(x = Arch_Label, y = Fail_Rate, fill = Schedule_Family)) +
  geom_col(width = 0.50, color = "#0F172A", linewidth = 0.5) +
  geom_errorbar(aes(ymin = pmax(0, Fail_Rate - SE_Fail), ymax = Fail_Rate + SE_Fail),
                width = 0.14, linewidth = 0.5, na.rm = TRUE) +
  geom_text(aes(label = sprintf("%.1f%% Failure", Fail_Rate),
                y = ifelse(is.na(SE_Fail), Fail_Rate + 0.9, Fail_Rate + SE_Fail + 0.9)),
            fontface = "bold", size = 3.5) +
  scale_fill_manual(values = c("Pure_GA" = "#0072B2", "Tsuchiwata_GA" = "#D55E00", "Tsuchiwata_Paper_7pt" = "#56B4E9")) +
  coord_cartesian(ylim = c(0, 14.5)) +
  labs(title    = "Clearance Estimability Failure Rate at Semi-Sparse Regime (N = 7)",
       subtitle = "Discontinuous 5% step-cliff penalty produces >3.2x higher rate of non-estimable clearance fits",
       x        = "Optimization Architecture",
       y        = "Non-Estimable Patients (%)",
       caption  = "(a)") +
  plot_theme_ieee(base_size = 11.5) +
  theme(legend.position = "none")

# (b) Parameter Recovery Trade-Off

df_fig5b <- schedule_summary %>%
  filter(
    Sampling_Threshold == "semi-sparse",
    Evaluation_Engine == "Direct_NCA",
    Schedule_Family %in% c("Pure_GA", "Tsuchiwata_GA"),
    Metric %in% c("Cmax", "Tmax", "CL_F", "AUC")
  ) %>%
  group_by(Schedule_Family, Metric) %>%
  summarise(Mean_MAPE = mean(MAPE), SE_MAPE = sd(MAPE) / sqrt(n()), .groups = "drop") %>%
  mutate(
    Metric_Order = factor(Metric, levels = c("Cmax", "Tmax", "CL_F", "AUC"),
                          labels = c("Cmax MAPE", "Tmax MAPE", "CL/F MAPE", "AUC MAPE")),
    Formulation  = factor(Schedule_Family, levels = c("Pure_GA", "Tsuchiwata_GA"),
                          labels = c("Pure GA (Quadratic Barrier, Ours)", "Tsuchiwata GA (Step-Cliff Baseline)"))
  )

p_fig5b <- ggplot(df_fig5b, aes(x = Metric_Order, y = Mean_MAPE, fill = Formulation)) +
  geom_col(position = position_dodge(0.70), width = 0.58, color = "#0F172A", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, Mean_MAPE - SE_MAPE), ymax = Mean_MAPE + SE_MAPE),
                position = position_dodge(0.70), width = 0.16, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Mean_MAPE), y = Mean_MAPE + SE_MAPE + 1.2),
            position = position_dodge(0.70), fontface = "bold", size = 3.4) +
  scale_fill_manual(values = c("Pure GA (Quadratic Barrier, Ours)" = "#0072B2",
                               "Tsuchiwata GA (Step-Cliff Baseline)" = "#D55E00"),
                    name = "Objective Formulation") +
  coord_cartesian(ylim = c(0, 32)) +
  labs(title    = "Parameter Recovery Trade-Off (N = 7)",
       subtitle = "Step-cliff over-indexes on peak timing (~7 min physical MAE diff) at the expense of CL/F & AUC",
       x        = "Pharmacokinetic Parameter",
       y        = "Mean Absolute Percentage Error (MAPE %)",
       caption  = "(b)") +
  plot_theme_ieee(base_size = 11.5)

fig5_vertical <- arrangeGrob(p_fig5a, p_fig5b, nrow = 2, ncol = 1, heights = c(1.0, 1.15))
save_pub_plot(fig5_vertical, subdirs$main, "Figure_05_IEEE_Objective_Function_Ablation", w = 8.8, h = 10.4)

# Appendix figures ====
# Figure A1: Cohort schedule heatmaps across schedule generation methods ====

message(">> Generating IEEE Appendix Figure A1...")

sched_for_heatmap <- master_index %>%
  filter(Evaluation_Engine == "Direct_NCA",
         Schedule_Family %in% c("Pure_GA", "PCHIP_49", "PopED", "Tsuchiwata_GA"),
         !is.na(Cohort_Match_ID))

heatmap_rows <- lapply(seq_len(nrow(sched_for_heatmap)), function(i) {
  r <- sched_for_heatmap[i, ]
  sel_t <- as.numeric(strsplit(r$Schedule_Times, ";")[[1]])
  data.frame(
    Schedule_Family    = r$Schedule_Family,
    Sampling_Threshold = r$Sampling_Threshold,
    Cohort             = factor(sprintf("%02d", r$Cohort_Match_ID), levels = sprintf("%02d", 10:1)),
    Time               = candidate_grid,
    Selected           = factor(ifelse(round(candidate_grid, 5) %in% round(sel_t, 5), "Selected", "Not Selected"),
                                levels = c("Not Selected", "Selected"))
  )
}) %>% bind_rows() %>%
  mutate(
    Family_Label = factor(Schedule_Family, levels = c("Pure_GA", "PCHIP_49", "PopED", "Tsuchiwata_GA"),
                          labels = c("Pure GA", "PCHIP-49 GA", "PopED (D-Optimal)", "Tsuchiwata GA")),
    Regime_Label = factor(Sampling_Threshold, levels = c("sparse", "semi-sparse", "intensive"),
                          labels = c("Sparse (N=4)", "Semi-Sparse (N=7)", "Intensive (N=15)"))
  )

p_app_a1 <- ggplot(heatmap_rows, aes(x = factor(Time), y = Cohort, fill = Selected)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(Family_Label + Regime_Label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Not Selected" = "#F1F5F9", "Selected" = "#D55E00"), name = "Sampling Point Status") +
  labs(x = "Candidate Sampling Timepoint (Hours)", y = "Training Cohort ID") +
  plot_theme_app_ieee(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5),
    axis.text.y   = element_text(size = 6.5),
    strip.text.y  = element_text(angle = 0, face = "bold", size = 7.5),
    panel.spacing = unit(0.22, "lines")
  )
save_pub_plot(p_app_a1, subdirs$appendix, "Appendix_Fig_A1_IEEE_Horizontal_Heatmaps", w = 9.2, h = 11.2)

# Figure B1: Master 6-metric PK parameter recovery boxplot ====

message(">> Generating IEEE Appendix Figure B1...")

p_app_b1 <- ggplot() +
  geom_boxplot(data = df_primary %>% filter(!Is_Static),
               aes(x = Family_Label, y = MAPE, fill = Schedule_Family),
               alpha = 0.45, outlier.shape = NA, width = 0.55) +
  geom_jitter(data = df_primary %>% filter(!Is_Static),
              aes(x = Family_Label, y = MAPE, color = Schedule_Family),
              width = 0.12, alpha = 0.75, size = 1.6) +
  stat_summary(data = df_primary %>% filter(!Is_Static),
               aes(x = Family_Label, y = MAPE),
               fun = mean, geom = "point", shape = 23, size = 2.8, fill = "white", stroke = 0.9) +
  geom_point(data = df_primary %>% filter(Is_Static),
             aes(x = Family_Label, y = MAPE, fill = Schedule_Family, color = Schedule_Family),
             shape = 23, size = 3.3, stroke = 1.1) +
  facet_grid(Metric ~ Regime_Label, scales = "free", labeller = labeller(Metric = metric_labels)) +
  scale_x_discrete(drop = TRUE, expand = expansion(add = c(0.55, 0.55))) +
  scale_fill_manual(values = palette_methods) +
  scale_color_manual(values = palette_methods) +
  labs(x = NULL, y = "Schedule-Level MAPE (%)") +
  plot_theme_app_ieee(base_size = 10.5) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_pub_plot(p_app_b1, subdirs$appendix, "Appendix_Fig_B1_IEEE_Master_6Metric", w = 9.0, h = 10.8)

# Figure B2: Elimination recovery boxplots ====

message(">> Generating IEEE Appendix Figure B2...")

build_elim_row_ieee <- function(metric_code, y_axis_label, y_max, y_breaks, header_line_y, badge_y, sub_label) {
  df_sub <- df_primary_all %>% filter(Metric == metric_code)
  df_fin <- df_sub %>% filter(is.finite(MAPE))
  
  df_badges <- df_sub %>%
    group_by(Regime_Label, Family_Label, Schedule_Family) %>%
    summarise(Mean_Fail = mean(Nonestimable_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Badge_Label = sprintf("%.1f%% fail", Mean_Fail),
      Badge_Fill  = ifelse(Mean_Fail >= 50, "#FEE2E2", ifelse(Mean_Fail >= 10, "#FEF3C7", "#F1F5F9")),
      Badge_Col   = ifelse(Mean_Fail >= 50, "#991B1B", ifelse(Mean_Fail >= 10, "#92400E", "#334155"))
    )
  
  df_na_notice <- df_badges %>% filter(Mean_Fail >= 99.9)
  
  ggplot() +
    geom_hline(yintercept = header_line_y, linetype = "dotted", color = "grey65", linewidth = 0.5) +
    geom_boxplot(data = df_fin %>% filter(!Is_Static),
                 aes(x = Family_Label, y = MAPE, fill = Schedule_Family),
                 alpha = 0.45, outlier.shape = NA, width = 0.52) +
    geom_jitter(data = df_fin %>% filter(!Is_Static),
                aes(x = Family_Label, y = MAPE, color = Schedule_Family),
                width = 0.12, alpha = 0.80, size = 1.9) +
    stat_summary(data = df_fin %>% filter(!Is_Static),
                 aes(x = Family_Label, y = MAPE),
                 fun = mean, geom = "point", shape = 23, size = 3.1, fill = "white", stroke = 1.0) +
    geom_point(data = df_fin %>% filter(Is_Static),
               aes(x = Family_Label, y = MAPE, fill = Schedule_Family, color = Schedule_Family),
               shape = 23, size = 3.6, stroke = 1.1) +
    geom_text(data = df_na_notice,
              aes(x = Family_Label, y = header_line_y * 0.48, label = "No Estimable\nProfiles\n(100% Failed)"),
              size = 2.8, fontface = "italic", color = "#991B1B", lineheight = 0.95) +
    geom_label(data = df_badges,
               aes(x = Family_Label, y = badge_y, label = Badge_Label),
               fill = df_badges$Badge_Fill, color = df_badges$Badge_Col,
               size = 2.75, fontface = "bold", label.size = 0.25,
               label.padding = unit(0.16, "lines")) +
    facet_wrap(~ Regime_Label, nrow = 1, scales = "free_x") +
    scale_x_discrete(drop = TRUE, expand = expansion(add = c(0.62, 0.62))) +
    scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = c(0, 0)) +
    scale_fill_manual(values = palette_methods) +
    scale_color_manual(values = palette_methods) +
    labs(x = NULL, y = y_axis_label, caption = sub_label) +
    plot_theme_app_ieee(base_size = 10.5) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 25, hjust = 1, size = 8.5, face = "bold")
    )
}

p_b2_row1_clf <- build_elim_row_ieee("CL_F", "CL/F MAPE (%)", 92, seq(0, 80, by = 20), 80, 86, "(a)")
p_b2_row2_thalf <- build_elim_row_ieee("T_half", "t1/2 MAPE (%)", 198, seq(0, 160, by = 40), 172, 185, "(b)")

df_tmax_h <- df_primary_all %>% filter(Metric == "Tmax", is.finite(MAE))

p_b2_row3_tmax <- ggplot() +
  geom_boxplot(data = df_tmax_h %>% filter(!Is_Static),
               aes(x = Family_Label, y = MAE, fill = Schedule_Family),
               alpha = 0.45, outlier.shape = NA, width = 0.52) +
  geom_jitter(data = df_tmax_h %>% filter(!Is_Static),
              aes(x = Family_Label, y = MAE, color = Schedule_Family),
              width = 0.12, alpha = 0.80, size = 1.9) +
  stat_summary(data = df_tmax_h %>% filter(!Is_Static),
               aes(x = Family_Label, y = MAE),
               fun = mean, geom = "point", shape = 23, size = 3.1, fill = "white", stroke = 1.0) +
  geom_point(data = df_tmax_h %>% filter(Is_Static),
             aes(x = Family_Label, y = MAE, fill = Schedule_Family, color = Schedule_Family),
             shape = 23, size = 3.6, stroke = 1.1) +
  facet_wrap(~ Regime_Label, nrow = 1, scales = "free_x") +
  scale_x_discrete(drop = TRUE, expand = expansion(add = c(0.62, 0.62))) +
  scale_y_continuous(breaks = seq(0, 2.5, by = 0.5), limits = c(0, 2.75), expand = c(0, 0)) +
  scale_fill_manual(values = palette_methods) +
  scale_color_manual(values = palette_methods) +
  labs(x = NULL, y = "Tmax MAE (Hours)", caption = "(c)") +
  plot_theme_app_ieee(base_size = 10.5) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 25, hjust = 1, size = 8.5, face = "bold")
  )

app_b2_rescaled <- arrangeGrob(p_b2_row1_clf, p_b2_row2_thalf, p_b2_row3_tmax, nrow = 3, ncol = 1, heights = c(1.05, 1.05, 0.95))
save_pub_plot(app_b2_rescaled, subdirs$appendix, "Appendix_Fig_B2_IEEE_Rescaled_Elimination_and_Tmax", w = 9.2, h = 11.4)

# Figure B3: Friedman omnibus and KS separation heatmaps ====

message(">> Generating IEEE Appendix Figure B3...")

# (a) Nonparametric Friedman Rank-Sum ANOVA Omnibus

df_fm <- friedman_method %>% filter(is.finite(p_adj_Holm)) %>%
  mutate(Metric_Display = factor(Metric, levels = rev(names(metric_labels)), labels = rev(metric_labels)),
         Regime_Label   = factor(Sampling_Threshold, levels = names(regime_labels), labels = regime_labels),
         NegLog10P      = -log10(pmax(p_adj_Holm, 1e-5)),
         Cell_Text      = paste0("Q = ", round(Q_Statistic, 2), "\np = ", vapply(p_adj_Holm, format_p, character(1))))

p_app_b3a <- ggplot(df_fm, aes(x = Regime_Label, y = Metric_Display, fill = NegLog10P)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Cell_Text), size = 3.0, lineheight = 0.95) +
  facet_wrap(~ Outcome, nrow = 1) +
  scale_fill_gradient(low = "white", high = "#B40426", name = "-log10(p_Holm)") +
  labs(x = NULL, y = NULL, caption = "(a)") +
  plot_theme_app_ieee(base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold"), legend.position = "right")

# (b) KS Test Heatmap (Pure GA vs. PopED)

df_ks_poped <- primary_tests %>% filter(Group_B == "PopED", Outcome == "MAPE", is.finite(KS_D)) %>%
  mutate(Metric_Display = factor(Metric, levels = rev(names(metric_labels)), labels = rev(metric_labels)),
         Regime_Label   = factor(Regime, levels = names(regime_labels), labels = regime_labels),
         Cell_Label     = paste0("D = ", formatC(KS_D, digits = 2, format = "f"), " | d = ", formatC(Mean_Delta, digits = 1, format = "f"), "%\np = ", vapply(Perm_P_Holm, format_p, character(1))))


p_app_b3b <- ggplot(df_ks_poped, aes(x = Regime_Label, y = Metric_Display, fill = KS_D)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Cell_Label), size = 2.9, lineheight = 0.95) +
  scale_fill_gradient(low = "white", high = "#3B4CC0", limits = c(0, 1), name = "KS D Stat") +
  labs(x = NULL, y = NULL, caption = "(b)") +
  plot_theme_app_ieee(base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold"), legend.position = "right")

# (c) KS Test Heatmap (Pure GA vs. PCHIP-49)

df_ks_pchip <- primary_tests %>% filter(Group_B == "PCHIP_49", Outcome == "MAPE", is.finite(KS_D)) %>%
  mutate(Metric_Display = factor(Metric, levels = rev(names(metric_labels)), labels = rev(metric_labels)),
         Regime_Label   = factor(Regime, levels = names(regime_labels), labels = regime_labels),
         Cell_Label     = paste0("D = ", formatC(KS_D, digits = 2, format = "f"), " | d = ", formatC(Mean_Delta, digits = 1, format = "f"), "%\np = ", vapply(Perm_P_Holm, format_p, character(1))))

p_app_b3c <- ggplot(df_ks_pchip, aes(x = Regime_Label, y = Metric_Display, fill = KS_D)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Cell_Label), size = 2.9, lineheight = 0.95) +
  scale_fill_gradient(low = "white", high = "#009E73", limits = c(0, 1), name = "KS D Stat") +
  labs(x = "Sampling Density Regime", y = NULL, caption = "(c)") +
  plot_theme_app_ieee(base_size = 10.5) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(face = "bold"), legend.position = "right")

save_pub_plot(arrangeGrob(p_app_b3a, p_app_b3b, p_app_b3c, nrow = 3, ncol = 1), subdirs$appendix, "Appendix_Fig_B3_IEEE_Heatmaps", w = 8.8, h = 11.4)

# Figure B4: Empirical cumulative error distributions and log10 patient APE distributions ====

message(">> Generating IEEE Appendix Figure B4...")

# (a) Empirical Cumulative Error Distribution

df_ecdf <- schedule_summary %>%
  filter(Sampling_Threshold == "semi-sparse", Evaluation_Engine == "Direct_NCA",
         Schedule_Family %in% c("Pure_GA", "PopED", "PCHIP_49"),
         Metric %in% c("AUC", "Cmax", "CL_F"), is.finite(MAPE)) %>%
  mutate(Family_Label   = factor(Schedule_Family, levels = method_order, labels = labels_methods[method_order]),
         Metric_Display = factor(Metric, levels = names(metric_labels), labels = metric_labels))

p_app_b4a <- ggplot(df_ecdf, aes(x = MAPE, color = Family_Label)) +
  stat_ecdf(linewidth = 1.05) +
  geom_rug(alpha = 0.55, sides = "b") +
  facet_wrap(~ Metric_Display, scales = "free_x") +
  scale_color_manual(values = setNames(palette_methods, labels_methods)) +
  labs(x = "Schedule-Level MAPE (%)", y = "Cumulative Proportion", color = "Method", caption = "(a)") +
  plot_theme_app_ieee(base_size = 10.5)

# (b) Unaggregated Patient-Level Dispersion

if (!is.null(subject_errors)) {
  pat_df <- subject_errors %>%
    filter(
      Evaluation_Engine == "Direct_NCA",
      (Schedule_Family %in% c("Pure_GA", "PopED") & Sampling_Threshold == "semi-sparse") |
        (Schedule_Family == "Clinical" & Sampling_Threshold == "intensive"),
      Metric %in% c("AUC", "Cmax"),
      is.finite(APE)
    ) %>%
    mutate(
      Family_Label = factor(
        Schedule_Family,
        levels = c("Pure_GA", "PopED", "Clinical"),
        labels = c("Pure GA\n(N = 7)", "PopED\n(N = 7)", "Clinical 15pt\n(N = 15)")
      ),
      Metric_Display = factor(Metric, levels = names(metric_labels), labels = metric_labels)
    )
  
  p_app_b4b <- ggplot(pat_df, aes(x = Family_Label, y = APE + 0.1, fill = Schedule_Family)) +
    geom_boxplot(outlier.size = 0.7, outlier.alpha = 0.25, width = 0.52, alpha = 0.6) +
    facet_wrap(~ Metric_Display, scales = "free_y") +
    scale_fill_manual(values = palette_methods) +
    scale_y_log10(labels = label_comma()) +
    labs(x = NULL, y = "Subject-Level Abs. Percentage Error + 0.1% (Log10 Scale)", caption = "(b)") +
    plot_theme_app_ieee(base_size = 10.5) +
    theme(legend.position = "none")
  
  save_pub_plot(
    arrangeGrob(p_app_b4a, p_app_b4b, nrow = 2, ncol = 1),
    subdirs$appendix,
    "Appendix_Fig_B4_IEEE_ECDF_and_Log_Patient_APE",
    w = 8.8, h = 10.2
  )
}

# Figure B5: Pareto frontier and Effect size forest plot ====

message(">> Generating IEEE Appendix Figure B5...")

df_tradeoff <- schedule_summary %>%
  filter(Evaluation_Engine == "Direct_NCA", Metric %in% c("CL_F", "Vz_F", "T_half"),
         Schedule_Family %in% c("Pure_GA", "PopED", "PCHIP_49", "Tsuchiwata_GA", "Tsuchiwata_Paper_7pt"),
         is.finite(MAPE), is.finite(Nonestimable_pct)) %>%
  mutate(Regime_Label = factor(Sampling_Threshold, levels = names(regime_labels), labels = regime_labels))

p_app_b5a <- ggplot(df_tradeoff, aes(x = Nonestimable_pct, y = MAPE, color = Schedule_Family, shape = Schedule_Family)) +
  geom_point(size = 2.8, alpha = 0.85) +
  facet_grid(Metric ~ Regime_Label, scales = "free_y", labeller = labeller(Metric = metric_labels)) +
  scale_color_manual(values = palette_methods, labels = labels_methods, name = "Method") +
  scale_shape_manual(values = c("Pure_GA" = 16, "PCHIP_49" = 17, "PopED" = 15, "Tsuchiwata_GA" = 18, "Tsuchiwata_Paper_7pt" = 8),
                     labels = labels_methods, name = "Method") +
  labs(x = "Patients Failing Terminal Regression (%)", y = "MAPE among Estimable Patients (%)", caption = "(a)") +
  plot_theme_app_ieee(base_size = 10.5)

df_forest <- primary_tests %>%
  filter(Group_B == "PopED", Outcome == "MAPE", is.finite(Mean_Delta)) %>%
  mutate(Metric_Display = factor(Metric, levels = rev(names(metric_labels)), labels = rev(metric_labels)),
         Regime_Label   = factor(Regime, levels = names(regime_labels), labels = regime_labels))

p_app_b5b <- ggplot(df_forest, aes(x = Mean_Delta, y = Metric_Display, xmin = CI_Lower, xmax = CI_Upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(height = 0.22, linewidth = 0.7, color = "#222222") +
  geom_point(size = 2.8, color = "#D55E00") +
  facet_wrap(~ Regime_Label, scales = "free_x") +
  labs(x = "Paired Mean Difference in MAPE (Percentage Points)", y = NULL, caption = "(b)") +
  plot_theme_app_ieee(base_size = 10.5)

save_pub_plot(arrangeGrob(p_app_b5a, p_app_b5b, nrow = 2, ncol = 1, heights = c(1.35, 1.0)), subdirs$appendix, "Appendix_Fig_B5_IEEE_Pareto_Forest", w = 8.8, h = 10.6)

message("All IEEE-Formatted Main (1-5) & Appendix (A1, B1-B5) Figures Saved In:")
message(normalizePath(plot_root))

