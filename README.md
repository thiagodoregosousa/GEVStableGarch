# GEVStableGarch

An R package for ARMA-GARCH and ARMA-APARCH estimation with GEV
(Generalized Extreme Value), stable, GAt, and skew Student's t conditional
distributions, in addition to the usual normal, Student's t, skew-t (from
`fGarch`) and GED innovations. Useful for modeling financial time series
with volatility clustering, fat tails, and skewness — e.g. for volatility
forecasting and Value-at-Risk estimation.

The package was previously on CRAN (last release: version 1.1, 2015) but
was archived by CRAN in 2020. The top level of this repo is being
reorganized into a standard `R/`-based package layout (see the note below),
so `remotes::install_github()` does not yet work directly against it. In
the meantime, install the last CRAN release straight from the CRAN
archive:

## Installation

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/GEVStableGarch/GEVStableGarch_1.1.tar.gz",
  repos = NULL, type = "source"
)
```

## Usage

```r
library(GEVStableGarch)

# Fit an ARMA(1,1)-GARCH(1,1) model with GEV conditional distribution
# to the dem2gbp exchange-rate returns (from the fGarch package)
data(dem2gbp)
x <- dem2gbp[, 1]

fit <- gsFit(data = x, formula = ~garch(1, 1), cond.dist = "gev")
fit
```

## Functions

*The package is currently being refactored, so the function names below may
still change. This table reflects what is exported today.*

| Function | Description |
|---|---|
| `gsFit()` | Maximum likelihood estimation of ARMA-GARCH/APARCH models, optionally restricted to the stationarity region |
| `gsSpec()` | Specifies an ARMA-GARCH/APARCH model with GEV or stable conditional distribution |
| `gsSim()` | Simulates a time series from an ARMA-GARCH/APARCH model |
| `gsSelect()` | Fits models over a range of orders and selects the best one by AIC/AICc/BIC |
| `gsMomentAparch()` | Computes the moment expression used to check APARCH stationarity conditions |
| `dgat()`, `pgat()`, `qgat()`, `rgat()` | Density, distribution, quantile, and random generation for the Generalized Asymmetric t (GAt) distribution |
| `dskstd()`, `pskstd()`, `qskstd()`, `rskstd()` | Density, distribution, quantile, and random generation for the skew Student's t distribution (Fernandez & Steel, 1998) |

`show` methods are also provided for the `GEVSTABLEGARCH` and
`GEVSTABLEGARCHSPEC` S4 classes returned by `gsFit()` and `gsSpec()`.

## Package layout

- Top-level `GEVStableGarch-*.R`, `class-*.R`, `dist-*.R`, `gevStableGarch-Spec*.R`
  and `methods-show.R` — package source. This is a legacy layout (no `R/`
  folder, no root `DESCRIPTION`/`NAMESPACE`) that predates the current
  refactor; the last buildable package tree is the one under `CRAN_versions/`
- `CRAN_versions/` — snapshots of the package as submitted to CRAN (each with
  its own `DESCRIPTION`, `NAMESPACE`, `R/`, `man/`), including the checked
  `.tar.gz` and `R CMD check` logs
- `tests/` — unit test scripts and saved test-run results
- `docs/` — changelog and supporting notes (e.g. on the filtering process
  used during estimation)
- `LICENSE` — GPL (>= 2)

## Reference

Wuertz, D., Chalabi, Y., Luksan, L. (2009). *Parameter Estimation of ARMA
Models with GARCH/APARCH Errors: An R and SPlus Software Implementation*.
Journal of Statistical Software.

Zhao, X., Scarrott, C.J., Oxley, L., Reale, M. (2011). *GARCH dependence in
extreme value models with Bayesian inference*. Mathematics and Computers in
Simulation, 81(7), 1430-1440.
