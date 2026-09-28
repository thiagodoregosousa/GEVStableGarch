# GEVStableGarch

An R package for ARMA-GARCH and ARMA-APARCH models with **stable**,
**GEV** (Generalized Extreme Value) and **GAt** (Generalized Asymmetric t,
Paolella 1997) innovations: estimation, simulation, stationarity checks and
forecasting (including Value-at-Risk). These three families are not offered
by fGarch, rugarch or tsgarch. Useful for financial time series with
volatility clustering, heavy tails and skewness.

The innovation distribution is a pluggable object: a user defined density
gets the same fitting, simulation, stationarity and forecasting machinery,
and the constructor checks that everything it declares is consistent, so
fitting never fails silently.

Version 2.0.0 is a rewrite of the version 1.1 released on CRAN in 2015
(archived by CRAN in 2020). The interface changed; see `NEWS.md`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("thiagodoregosousa/GEVStableGarch", build_vignettes = TRUE)
```

## Usage

```r
library(GEVStableGarch)

# AR(1)-APARCH(1,1) with stable innovations
spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())

# Simulate from known parameters, then fit
model <- gs_model(spec, mu = 0.0005, ar = 0.05, omega = 0.0002, alpha = 0.05,
                  gamma = 0.3, beta = 0.85, delta = 1.2,
                  dist_par = c(stable_alpha = 1.75, stable_beta = 0))
x <- gs_sim(model, n = 1500, seed = 1)$y
fit <- gs_fit(x, spec)
fit

# Persistence (below 1: stationary scale level)
gs_stationarity(fit)

# One step location, scale and quantiles (VaR is minus the quantile)
predict(fit, n_ahead = 1, level = c(0.01, 0.05))

# Probability integral transform of the residuals (e.g. for copulas)
u <- residuals(fit, type = "pit")
```

See `vignette("GEVStableGarch")` for forecasting, rolling windows and user
defined distributions.

## Functions

| Function | Description |
|---|---|
| `gs_spec()` | Specifies an ARMA(m,n)-GARCH/APARCH(p,q) model and its innovation distribution |
| `gs_model()` | Attaches parameter values to a specification, with admissibility checks |
| `gs_fit()` | Maximum likelihood estimation, warm starts, standard errors, stationarity report |
| `gs_sim()` | Simulates a series from a model |
| `predict()` | Forecasts location, mean, scale and quantiles; exact one step, simulated beyond |
| `gs_stationarity()` | Persistence of the APARCH recursion |
| `gs_aparch_moment()` | $E(\|z\| - \gamma z)^\delta$ of the innovations (closed form or numerical) |
| `gs_stable()`, `gs_gev()`, `gs_gat()` | Built in innovation distributions |
| `gs_dist()`, `gs_check_dist()` | User defined innovation distributions and their self check |
| `dgat()`, `pgat()`, `qgat()`, `rgat()` | Density, distribution, quantile and random generation of the GAt distribution |

Methods for fitted models: `coef()`, `vcov()`, `logLik()`, `AIC()`, `BIC()`,
`residuals(type = c("raw", "standardized", "pit"))`, `sigma()`, `fitted()`,
`print()`.

## Package layout

- `R/`: package source (one file per concern: distributions, parameters,
  filters, likelihood, fit, simulation, forecasting)
- `tests/testthat/`: unit tests, including cross checks against fGarch
- `vignettes/`: worked example
- `dev/`: development material excluded from the build (restructuring plan
  notes, legacy code of version 1.1, benchmarks, exploratory scripts)
- `PLAN.md`: restructuring plan and decisions

## Reference

Mittnik, S., Paolella, M. S., Rachev, S. T. (2002). *Stationarity of stable
power-GARCH processes*. Journal of Econometrics, 106, 97-107.

Paolella, M. S. (1997). *Tail estimation and conditional modeling of
heteroskedastic time series*. PhD thesis, University of Kiel.

Wuertz, D., Chalabi, Y., Luksan, L. (2006). *Parameter Estimation of ARMA
Models with GARCH/APARCH Errors: An R and SPlus Software Implementation*.
Journal of Statistical Software.

Zhao, X., Scarrott, C.J., Oxley, L., Reale, M. (2011). *GARCH dependence in
extreme value models with Bayesian inference*. Mathematics and Computers in
Simulation, 81(7), 1430-1440.
