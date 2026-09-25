# GEVStableGarch restructuring plan

Status date: 2026-09-25. Branch: `package-hygiene`.

## 1. Goal

A minimal R package that fits and simulates ARMA(m,n)-GARCH/APARCH(p,q) models with three innovation families that no mainstream GARCH package offers: **stable (S0 only)**, **GEV** and **GAT** (Paolella). Everything else is either delegated to other packages or left out.

Principles:

* Only stable, GEV and GAT are public distributions. Normal stays internal, for tests only.
* Stable density and random numbers come from `libstable4u` (faster). No `stabledist`, no `stable`, one parametrization (S0).
* Plain lists with S3 classes instead of S4 classes. Roxygen generates NAMESPACE and man; testthat runs the tests. No hand-written Rd.
* The innovation distribution is a pluggable object, so users can supply their own density and reuse the whole ARMA-APARCH machinery.

## 2. Where the code is today

The January 2026 work (commits `3e8b986` to `dbb4383`) started the migration but stopped halfway:

* **Done**: package layout (`R/`, `tests/testthat/`), roxygen for `dist_gat.R`, GAT tests and Monte Carlo benchmark, `libstable4u` swapped into `.armaGarchDist`, new list based `GEVStableGarch_model_spec()` and `GEVStableGarch_set_params()` with tests.
* **Broken**: `dbb4383` deleted `gsSpec`, `.getFormula`, `.getOrder`, `gsSelect`, the S4 classes and show methods, but `gsFit` still calls `.getFormula()` and `new("GEVSTABLEGARCH")`. The package does not load a working fit today.
* **Half written**: `gsSim()` mixes the new list API with old S4 slots (`spec@distribution`), uses undefined `h`, `m`, `rseed`, `n.start`.
* **Bugs found while reading**:
  * `GEVStableGarch_set_params()` checks `class(spec) != "GEVStableGarch_model"`, should be `"GEVStableGarch_spec"`. Every call errors.
  * `GEVStableGarch_model_spec()` accepts `"norm"` but rejects `"gat"`, and the check is duplicated.
  * `R/GEVStableGarch_fit.R` ends with a top level `GSgarch.dstable <<- ...` using `stabledist`, a side effect at load time.
  * `.armaGarchDist` still has a `std` branch calling `dstd` (fGarch).
  * `R/GEVStableGarch_stationarity_aparch.R` defines `.gevMomentAparch` twice and still carries std, sstd, skstd, ged moments using fGarch/fExtremes functions.
  * `.filterGarch11Fit` prints on every call.
  * testthat files still call the deleted `gsSpec()` and old `gsSim(spec=)`; `test_sim.R` is a placeholder.
  * `DESCRIPTION`: depends on `fGarch`, `fExtremes`, `stabledist`, `skewt`, `timeDate`, `timeSeries`; `libstable4u` is missing; BugReports URL points to an old account.
  * `man/` is gitignored, so `remotes::install_github()` would ship no help pages.
* Open TODO items carried over: GAT stationarity region unknown (simulations explode), `gamma = 0` issue in APARCH estimation, uninformative Hessian error, `gsFit` too large.

Reference material now lives in `dev/` (excluded from the build): old R files recovered from `fccc2dc` in `dev/legacy/R`, old Rd files in `dev/legacy/man`, old test cases in `dev/legacy/test_cases`, the handwritten filtering notes, GAT benchmark, and 2025/2026 exploratory scripts (fast stable via `libstable4u`, ABEV3 fits compared with rugarch).

## 3. What is worth doing and what is not

Worth doing (brings value no other package gives):

* **Stable S0, GEV and GAT innovations in ARMA-APARCH.** fGarch, rugarch, tseries and tsgarch offer none of these three. This is the reason the package exists.
* **Stationarity conditions for APARCH under these innovations** (`E(|z| - gamma z)^delta`, `gsMomentAparch`) and **stationarity constrained estimation** (`sqp.restriction`). Few implementations exist, especially for stable innovations (Mittnik, Paolella, Rachev 2002).
* **Pluggable innovation distribution.** A small object (`log density`, `random generator`, `parameter names`, `bounds`, `start`, optional `APARCH moment`) lets users fit ARMA-APARCH with any density. rugarch and tsgarch have closed distribution lists. This is cheap to build because the three built-in families would use the same interface, and it is the main feature that makes the package more than a copy.
* **Fast stable likelihood via `libstable4u`.** Stable density evaluation dominates fitting time, so this is where speed comes from.

Not worth doing:

* **Normal, t, skew t, GED and similar innovations** as public options. fGarch and rugarch already do it well; keep normal internal for cross-checks.
* **Several stable parametrizations** (S1, S2) and the optional Nolan `stable` package.
* **The "matrix filtering" as a selling point.** The handwritten notes build the lagged vectors and call `stats::filter(..., "recursive")`. That is the same vectorized filter fGarch uses (Wuertz, Chalabi, Luksan 2006, the `.filterGarch11Fit` code even credits it). Keep it as the internal implementation, since it is correct and fast in R, but do not advertise it as new. A C++ recursion would only matter after the density is fast, and is out of scope for now.
* **Formula parsing** (`~arma(1,1)+aparch(1,1)`) and **`gsSelect`** (AIC grid search). The list spec replaces the formula; model selection is a three line loop for users. Can come back later if asked.
* **S4 classes and show methods.** Replace by S3 `print`, `coef`, `logLik`, `residuals`, `sigma`, `AIC`/`BIC` via `logLik`.
* **Forecasting.** Useful (stable prediction per Brockwell and Mittnik/Paolella) but a separate phase after the core is stable.

## 4. Target design

Public API (names follow the `snake_case` decision in TODO; prefix `gs_` proposed, see open question 1):

```r
spec  <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = FALSE,
                 include_mean = TRUE, dist = gs_stable())
model <- gs_model(spec, mean = 0, ar = 0.1, omega = 0.05,
                  alpha = 0.1, beta = 0.85, dist_par = list(alpha = 1.8, beta = 0))
x     <- gs_sim(model, n = 2000, burnin = 1000)
fit   <- gs_fit(x$y, spec, algorithm = c("sqp", "sqp_stationary", "nlminb"))
coef(fit); logLik(fit); residuals(fit); sigma(fit)
gs_stationarity(model)          # value < 1 means stationary
```

Distribution objects (`R/dist_*.R`):

```r
gs_stable()  gs_gev()  gs_gat()          # built in
gs_dist(name, log_density, random, par_names, lower, upper, start,
        aparch_moment = NULL)           # user defined, same structure
```

Each object stores `log_density(z, par)` for standardized innovations, so the likelihood is `sum(log_density(z / sigma, par) - log(sigma))` for every family. `aparch_moment(gamma, delta, par)` is optional; without it `gs_stationarity` falls back to numerical integration of the density.

Internal pieces:

* `filter_arma()`, `filter_aparch()`: the existing vectorized filters, cleaned.
* `neg_loglik(par, data, spec)`: unpacks a named parameter vector built from the spec (replaces the index arithmetic in `gsFit`).
* `start_values(data, spec)`: from `.getStart`.
* `gs_fit()` split into: build parameter vector and bounds, optimize, compute Hessian and standard errors, assemble S3 object.

Dependencies after cleanup: `Imports: Rsolnp, libstable4u, stats`; `Suggests: testthat, fGarch, rugarch`. GEV density written directly (already done in `.armaGarchDist`), GEV random numbers via inverse CDF.

## 5. Phases

Each phase is one `feat/*` branch merged into `master` with passing `devtools::check()`.

0. **package-hygiene** (this branch). Restore lost work, move legacy material to `dev/`, write this plan. Remaining: fix DESCRIPTION dependencies, remove the load time side effect, commit `man/`, add GitHub Actions `R-CMD-check`. Goal: `devtools::load_all()` and `devtools::check()` run (tests may be skipped).
1. **feat/dist-objects.** `gs_stable`, `gs_gev`, `gs_gat`, `gs_dist`, internal normal. Tests: densities integrate to 1, random generator matches density (KS test), `libstable4u` S0 matches reference values.
2. **feat/spec-model.** Rename and fix `model_spec`/`set_params` into `gs_spec`/`gs_model`, parameter vector packing and unpacking. Tests from `test_model_spec.R`.
3. **feat/sim.** `gs_sim` with burn-in and seed. Tests: moments and ACF of simulated normal GARCH(1,1) against theory; no NaN for stationary parameters.
4. **feat/fit.** Filters, `neg_loglik`, `start_values`, `gs_fit`, S3 methods. Tests: normal GARCH(1,1) on `dem2gbp` matches fGarch; simulate and fit round trip for stable, GEV, GAT on (m,n,p,q) grid (skip on CRAN); edge cases (stable alpha near 2 vs normal, GAT large nu).
5. **feat/stationarity.** Port `gsMomentAparch` for the three families, numerical fallback for user distributions, `sqp_stationary` algorithm. Investigate GAT stationarity region.
6. **docs/release.** README example, one vignette (custom distribution example), NEWS, version 2.0.0, CRAN submission.

## 6. Branch cleanup

After this branch merges into `master`, delete `develop`, `feature/repository-organization`, `feature/organize-tests` and `doc/readme_formatting` (all fully contained in `master` by then). Work model: `master` always passes check; short lived `feat/*`, `fix/*`, `docs/*`, `chore/*` branches.

## 7. Open questions

1. Function prefix: `gs_` (short, keeps the old `gs` identity) or `gevstable_`? Old exported names (`gsFit`, `gsSim`, `gsSpec`) break either way; CRAN version 1.1 users would need a major version bump.
2. Keep `sqp.restriction` style stationarity constrained fitting as a first class option, or just report stationarity after fitting?
3. Should forecasting enter the scope of version 2.0?
4. Is `rugarch`'s custom distribution route worth a quick test to confirm it cannot host these densities? (It is believed not to support user densities in `ugarchfit`.)
