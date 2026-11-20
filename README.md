# Direct Empirical Optimization of Clinical Pharmacokinetic Sampling Schedules Using Regularized Genetic Algorithms

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.3.0-276DC3.svg)](https://www.r-project.org/)

Official R codebase, Monte Carlo simulation datasets, schedule optimization pipelines, nonparametric statistical suites, and publication figure generators for:

> **Direct Empirical Optimization of Clinical Pharmacokinetic Sampling Schedules Using Regularized Genetic Algorithms**  
> Amirtha Maria — *University of California, San Francisco (UCSF)*

## Overview

Regulatory pharmacokinetic (PK) exposure and bioequivalence assessments rely on empirical **Non-Compartmental Analysis (NCA)**, which conventionally requires 15–20 blood draws per patient. While analytical D-optimal frameworks (e.g., **PopED**) maximize Fisher Information Matrix ($\det \text{FIM}$) sensitivities, they act as indirect parametric proxies that over-allocate draws to steep early absorption/distribution flanks (e.g., $0.25\text{ h}$) and leave too few post-peak observations for log-linear $\lambda_z$ regression under sparse ($N=4$) and semi-sparse ($N=7$) phlebotomy budgets.

This repository implements a **quadratically regularized binary Genetic Algorithm (`Pure GA`)** that optimizes sampling schedules directly against empirical NCA recovery metrics using:
- **In-loop vectorized trapezoidal integration** and empirical peak extraction ($C_{\max}$, $T_{\max}$, $AUC_{0\text{-}96}$).
- **A continuous one-sided quadratic tolerance barrier** replacing discontinuous "lethal chromosome" step-cliffs.
- **Precision-weighted loss formulation** (75% RMSPE / 25% MPE) with $T_{\max}$ attenuation ($w_{T_{\max}} = 0.025$).
- **Unconstrained multi-seed frequency consensus** for benchmark validation—demonstrating organic terminal-tail discovery—paired with an optional `enforce_terminal_anchor` configuration flag (`FALSE` by default) for prospective clinical protocol deployment.

---

## Repository Structure

All scripts use relative paths anchored at the repository root (`./data`, `./schedules`, `./results`). Pre-computed intermediate schedules and summary tables are included so any stage of the pipeline can be reproduced independently.

```text
pk-sampling-regularized-ga/
├── README.md
├── LICENSE
├── 01_cohort_generation.R          # Simulates N=1,000 PK subjects & 10 bootstrap cohorts (K=24)
├── 02_run_pure_ga.R                # Direct empirical Pure GA optimizer (N = 4, 7, 15)
├── 03_run_pchip_ga_49.R            # PCHIP-assisted GA on dense 49-pt grid (N = 4, 7, 15)
├── 04_popED_all_frequencies.R      # PopED D-optimal FIM benchmark generator (N = 4, 7, 15)
├── 05_run_pchip_ga_15.R            # PCHIP-assisted GA on discrete 15-pt clinical grid (N = 4, 7, 15)
├── 06_run_tsuchiwata_ga.R          # Tsuchiwata & Tsuji (2023) step-cliff objective benchmark (N = 7)
├── 07_data_analysis.R              # Out-of-sample NCA evaluation (N=1,000) and permutation/Friedman tests
├── 08_ieee_fig_generator.R         # Generates publication Main (Figs 1–5) & Appendix (Figs A1, B1–B5) plots
├── data/                           # Simulated 1,000-patient population & 10 training cohort CSVs
├── schedules/                      # Archived binary schedule matrices & selection frequency summaries
|── results/
│   ├── tables/                     # Master schedule index, error summaries, and statistical test tables
│   └── figures/
│       ├── 01_main_text_figures/   # Main manuscript figures (png/ and pdf/ subfolders)
│       └── 02_appendix_figures/    # Appendix figures (png/ and pdf/ subfolders)
└── paper/
    └── maria_2026_regularized_ga_pk_sampling.pdf   # IEEE-formatted manuscript
```
---

## Requirements and Dependencies

* R (($$\ge 4.3.0$$)
* Required CRAN packages (automatically checked and installed at the top of each script):
```r
install.packages(c(
  "deSolve", "MASS", "dplyr", "tidyr", "ggplot2", "scales",
  "GA", "doParallel", "foreach", "ncappc", "pracma", "gridExtra", "PopED"
))
```
---
## Reproducing the Pipeline

Set your R working directory to the root of this repository (`setwd("path/to/pk-sampling-regularized-ga")`) before running any script. 

### Option A: Regenerate Figures Only
Uses the pre-computed statistical tables in `results/tables/`:
```r
source("08_ieee_fig_generator.R")
```
*Note: Appendix Figure B.4(b) requires `results/tables/09_subject_errors_long.rds`. If that file is not stored locally, run Option B first to regenerate it.*

### Option B: Re-Run Out-of-Sample Evaluation and Statistical Tests
Evaluates all schedules in `schedules/` across all $N=1000$ virtual patients in `data/` via `ncappc` and executes exact sign-flip permutation tests ($2^{10} = 1{,}024$ permutations), 20,000-replicate bootstrap CIs, paired permutation KS tests, and Friedman ANOVAs:
```r
source("07_data_analysis.R")
source("08_ieee_fig_generator.R")
```

### Option C: Full End-to-End Execution from Simulation
Runs Monte Carlo population generation, multi-core GA/PCHIP/PopED schedule discovery across all 10 cohorts, out-of-sample benchmarking, and figure rendering:
```r
source("01_cohort_generation.R")
source("02_run_pure_ga.R")
source("03_run_pchip_ga_49.R")
source("04_popED_all_frequencies.R")
source("05_run_pchip_ga_15.R")
source("06_run_tsuchiwata_ga.R")
source("07_data_analysis.R")
source("08_ieee_fig_generator.R")
```
### Reproducibility and Prospective Use Notes

The exact converged schedules evaluated in the manuscript (including cohort-matched PopED benchmarks) are archived in `schedules/` for deterministic reproduction via Option B; `04_popED_all_frequencies.R` is configured with nominal population priors as the reference PopED script.

In `02_run_pure_ga.R`, `03_run_pchip_ga_49.R`, and `05_run_pchip_ga_15.R`, `enforce_terminal_anchor` defaults to `FALSE` to match the unconstrained manuscript benchmarks and may be set to `TRUE` for prospective clinical protocol generation (Section V.E of manuscript).

---
## Key Findings

1. **Direct Empirical Optimization Outperforms D-Optimal Proxies:** At the semi-sparse clinical target ($N=7$), Pure GA reduces blood draws by 53.33% relative to standard 15-point clinical schedules while achieving 1.85% $C_{\max}$ MAPE, 11.39% $AUC_{0\text{-}96}$ MAPE (8.29% median; $-2.41\%$ mean bias), and 96.83% apparent clearance ($CL/F$) estimability—outperforming PopED (11.71% $C_{\max}$ MAPE, $p_{\text{adj}} = 0.0117$; 73.83% estimability at $N=7$ and 0.00% at $N=4$, where Pure GA rescues 68.19% estimability). With zero hard-coded terminal-window rules, Pure GA organically discovered $72\text{-}96\text{ h}$ draws in 70% ($N=4$) and 80% ($N=7$) of cohorts.
2. **The Spline Paradox (($2 \times 2$ Factorial Swap):** Linear-space PCHIP compromises terminal recovery twice: upstream, in-loop linear extrapolation masked truncation error and caused both PCHIP-49 and PCHIP-15 to drop to 0/10 terminal draws ($72\text{-}96\text{ h}$) at $N=4$; downstream, post-hoc PCHIP reconstruction introduced cubic bowing and zero-slope clamping into exponential decay tails, inflating $t_{1/2}$ MAPE by $+174.68\%$ ($p_{\text{adj}} = 0.0117$) regardless of whether PCHIP was embedded inside the optimization loop.
3. **Target Grid Invariance:** Restricting in-loop PCHIP reconstruction to a discrete 15-point clinical grid (PCHIP-15) matches the dense 49-point simulation grid (PCHIP-49) across all six endpoints at both $N=4$ ($p_{\text{adj}} \ge 0.2344$) and $N=7$ ($p_{\text{adj}} \ge 0.7031$).
4. **Regularized Barrier vs. Step-Cliff Penalty:** Under identical unconstrained consensus aggregation, both Pure GA and the re-implemented Tsuchiwata & Tsuji (2023) step-cliff GA spontaneously placed a $72\text{-}96\text{ h}$ draw in 8/10 semi-sparse cohorts. However, by preventing peak-chasing ($3.70$ vs. $2.50$ draws past $4\text{ h}$), Pure GA's quadratically regularized loss reduced non-estimable clearance failures $>3$-fold (3.17% vs. 10.42%) while improving both $CL/F$ and $AUC_{0\text{-}96}$ accuracy.
---
## Citation

If you use this codebase, simulation data, or methodology in your research, please cite:

```bibtex
@mastersthesis{maria2026pksampling,
  author       = {Maria, Amirtha},
  title        = {Direct Empirical Optimization of Clinical Pharmacokinetic Sampling Schedules Using Regularized Genetic Algorithms},
  year         = {2026},
  institution  = {University of California, San Francisco},
  url    = {[https://github.com/mariaamirtha/pk-sampling-regularized-ga](https://github.com/mariaamirtha/pk-sampling-regularized-ga)}
```
---
## Contact

For questions about the code or reproducibility, please open a [GitHub Issue](https://github.com/mariaamirtha/pk-sampling-regularized-ga/issues). For research or professional inquiries, feel free to connect on [LinkedIn](https://www.linkedin.com/in/amirtha-maria-42182a1a6).

---
## License
This project is released under the MIT License.
