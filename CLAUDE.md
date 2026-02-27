# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

An R package (`probe`) implementing two Bayesian sparse high-dimensional linear regression algorithms:

- **PROBE** (PaRtitiOned empirical Bayes ECM): for homoscedastic models
- **H-PROBE**: extends PROBE to handle heteroscedastic (non-constant variance) models

The package is published on CRAN and GitHub at `alexmclain/PROBE` (subdir `probe`). It requires compilation of C++ code via `Rcpp`/`RcppArmadillo`.

## Repository Layout

```
probe/              # The R package source
  R/                # R source files
  src/              # C++ source files (Rcpp/RcppArmadillo)
  data/             # .rda datasets for examples and tests
  man/              # Roxygen2-generated documentation
  DESCRIPTION       # Package metadata (version 1.2)
Drug Response Example/   # Reproducible real-data example scripts
Simulation functions/    # Simulation study scripts
probe-manual.pdf    # Full PDF manual
```

## Key R Functions

| Function | File | Role |
|---|---|---|
| `probe()` | `R/probe_wrapper.R` | Main entry point for PROBE (all-at-once variant) |
| `hprobe()` | `R/hprobe.R` | Main entry point for H-PROBE |
| `predict_probe_func()` | `R/predict_probe_func.R` | Predictions + credible/prediction intervals from PROBE |
| `predict_hprobe_func()` | `R/predict_hprobe_func.R` | Predictions + intervals from H-PROBE |
| `e_step_func()` | `R/e_step_func.R` | Empirical Bayes E-step (estimates γ posteriors via local FDR) |
| `m_step_regression()` | `R/m_step_func.R` | M-step for α₀ and σ² (calibration model) |

## Algorithm Architecture

Both algorithms share the same PX-ECM loop structure:

1. **M-step (β update)**: Calls C++ functions (`PROBE_cpp0_5_6`, `PROBE_cpp0_5_6_covs` for PROBE; `*_h` variants for H-PROBE) that do column-wise linear regression. First iteration uses marginal regressions; subsequent iterations use the partitioned update.
2. **Damping step**: New estimates are blended with old via `fact = 1/(count+1)` to stabilize convergence.
3. **E-step**: `e_step_func()` computes t-statistics, estimates π₀ (null proportion) via Storey's method, fits a kernel density to the marginals, and returns local FDR values → γ (posterior inclusion probabilities).
4. **W update**: `m_update_func()` computes E[W] and E[W²] using the current γ and β̃.
5. **Calibration model**: `m_step_regression()` regresses Y on the aggregated signal W̃ to estimate α₀ (intercept/scale) and σ².
6. **Convergence**: Checked via a χ² statistic on the change in W̃ between iterations.

**H-PROBE difference**: The variance model uses `V %*% ω` (log-linear) with `Σ_y = diag(exp(-V %*% ω))`. The ω coefficients are estimated by BFGS optimization inside `m_step_regression.h()`. All C++ calls receive `Sigma_y_inv` instead of a scalar `sigma2`.

### C++ Files (`src/`)

| File | Function |
|---|---|
| `PROBE_cpp0_5_6.cpp` | All-at-once M-step, no covariates Z |
| `PROBE_cpp0_5_6_covs.cpp` | All-at-once M-step, with Z covariates |
| `PROBE_cpp0_5_6_h.cpp` | H-PROBE M-step, no Z |
| `PROBE_cpp0_5_6_covs_h.cpp` | H-PROBE M-step, with Z |
| `lm_by_col.cpp` / `lm_by_col_h.cpp` | Initial marginal regressions |
| `lm_w_covs_by_col.cpp` / `*_h.cpp` | Initial marginal regressions with Z |
| `PROBE_one_cpp.cpp` | One-at-a-time variant (alternative) |
| `mvm.cpp` | Matrix-vector multiply utility (MVM) |
| `row_sum.cpp` | Row sums utility |
| `inv_sigma.cpp` | Matrix inversion utility (inv_cpp) |

## Common Development Tasks

### Install/reinstall the package locally

```r
# From the repo root
devtools::install("probe")
# Or from the probe/ subdirectory:
R CMD build probe
R CMD INSTALL probe_1.2.tar.gz
```

### Recompile C++ after editing src/ files

```r
# From inside the probe/ directory
devtools::clean_dll()
devtools::load_all()
```

### Rebuild documentation (after editing roxygen comments)

```r
devtools::document("probe")
```

### Run the package examples

```r
library(probe)
# PROBE example
data(Sim_data); data(Sim_data_test)
attach(Sim_data); attach(Sim_data_test)
full_res <- probe(Y = Y, X = X, Y_test = Y_test, X_test = X_test,
                  alpha = 0.05, plot_ind = TRUE)
pred_res <- predict_probe_func(full_res, X = X_test)

# H-PROBE example
data(h_Sim_data); attach(h_sim_data)
res <- hprobe(Y = Y, X = X, V = V)
pred_res <- predict_hprobe_func(res, X_test, V = V_test)
```

### Check the package

```bash
R CMD check probe
```

## Key Parameters

- `adj`: Bandwidth multiplier for the kernel density in the E-step (default 5). Raise if you see "bandwidth may be too narrow" warnings.
- `ep`: Convergence threshold (default 0.1); compared to the χ² convergence criterion.
- `Z`: Optional non-sparse covariates not subject to variable selection.
- `V`: (H-PROBE only) Design matrix for the variance model.
- `signal` / `eta_i`: Only for simulation diagnostics/plotting; do not affect estimates.

## Output Structure (from `probe()` / `hprobe()`)

| Field | Description |
|---|---|
| `beta_ast_hat` | MAP estimates of regression coefficients (β* = γ × β̃ × α₁) |
| `gamma_hat` | Posterior inclusion probabilities |
| `sigma2_est` | MAP estimate of residual variance |
| `E_step` | Full E-step output (lfdr, p_vals, T_vals, pi0, etc.) |
| `Calb_mod` | Calibration model output (coef, VCV, res_data) |
| `count` | Number of iterations until convergence |
| `conv` | Convergence flag: 0=converged, 1=did not converge, 2=optimization failed |
| `omega` | (H-PROBE only) Variance model coefficients |
