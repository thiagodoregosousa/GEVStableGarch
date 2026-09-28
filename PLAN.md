# GEVStableGarch restructuring plan

Status date: 2026-09-28. Branch: `package-hygiene`. Revision 4 (author decisions of 2026-09-28 applied on top of review rounds 1 and 2).

## 1. Goal

A minimal R package that fits, simulates and forecasts ARMA(m,n)-GARCH/APARCH(p,q) models with three innovation families that no mainstream GARCH package offers: **stable (S0 only)**, **GEV** and **GAT** (Paolella). Everything else is delegated to other packages or left out.

Principles:

* Only stable, GEV and GAT are built in. Normal is shipped as the worked example of a user defined distribution, and used as the test oracle against fGarch.
* Stable density and random numbers come from `libstable4u` (faster). No `stabledist` or `stable` at run time, one parametrization (S0).
* Plain lists with S3 classes instead of S4 classes. Roxygen generates NAMESPACE and man; testthat runs the tests. No hand written Rd.
* The innovation distribution is a pluggable object. Its value is not the density hook alone (any `optim` script can do that) but everything that comes with it: simulation, forecasting, stationarity checks, standard errors and S3 methods, for any density. A user distribution must declare everything the machinery needs; incomplete or inconsistent definitions fail loudly at construction, never silently during fitting.

## 2. Where the code is today

The January 2026 work (commits `3e8b986` to `dbb4383`) started the migration but stopped halfway.

**Done**: package layout (`R/`, `tests/testthat/`), roxygen for `dist_gat.R`, GAT tests and Monte Carlo benchmark, `libstable4u` swapped into `.armaGarchDist`, list based `GEVStableGarch_model_spec()` and `GEVStableGarch_set_params()` with tests.

**Broken**: `dbb4383` deleted `gsSpec`, `.getFormula`, `.getOrder`, `gsSelect`, the S4 classes and show methods, but `gsFit` still calls `.getFormula()` (`fit.R:118`) and `new("GEVSTABLEGARCH")` (`fit.R:540`). No working fit today.

**Half written**: `gsSim()` mixes the new list API with old S4 slots (`spec@distribution`); `h`, `m`, `mu`, `ar`, `ma`, `omega`, `rseed`, `n.start` and `rgev` are undefined; `burnin` is ignored; it still uses `stabledist::rstable`.

**Bugs found** (verified by review):

* `GEVStableGarch_set_params()` checks `class(spec) != "GEVStableGarch_model"`, should be `"GEVStableGarch_spec"`; every call errors (`model_spec.R:127`).
* `GEVStableGarch_model_spec()` accepts `"norm"` but rejects `"gat"`; the check is duplicated (`model_spec.R:30-36`).
* Top level `GSgarch.dstable <<- ...` at load time using `stabledist` and `stable` (`fit.R:554-568`).
* `.armaGarchDist` still has a `std` branch calling `dstd` from fGarch (`armaGarchDist.R:91`); the GEV branch never checks the support `1 + xi z > 0` and returns NaN outside it (`armaGarchDist.R:124-127`).
* `.armaDist` calls `stable::dstable.quick` (not on CRAN) and drops zero densities from the likelihood, which biases it (`armaDist.R:110-118`).
* `stationarity_aparch.R` defines `.gevMomentAparch` twice (lines 389 and 580), keeps std, sstd, skstd, ged moments, and the S0 stable moment integrates `stable::dstable.quick` (line 531).
* Debug `print()` calls in `filter.R:177`, `dist_gat.R:277`, `get_start.R:264`, and several in `fit.R`.
* Parameter bounds allow the APARCH power `delta` up to 2 for stable (`get_start.R:230`) although `E|z|^delta` is finite only for `delta < alpha`.
* testthat files call the deleted `gsSpec()` and old `gsSim(spec=)`; they need `stabledist`; `test_sim.R` is a placeholder.
* `DESCRIPTION`: depends on `fGarch`, `fExtremes`, `stabledist`, `skewt`, `timeDate`, `timeSeries`, `R >= 2.15`; `libstable4u` missing; BugReports URL points to an old account; `Date` 2015.
* `man/` is gitignored, so `remotes::install_github()` would ship no help pages.

**Duplicated and scattered code** (to be removed in the restructuring, section 5):

* The parameter vector index arithmetic (`parm[(2+m+n+1):(3+m+n+p-1)]` etc.) is written three times in `gsFit` (LLH, stationarity function, output naming) and again implicitly in `.getStart`.
* Per distribution `if (cond.dist == ...)` ladders are repeated in `.armaGarchDist`, `.armaDist`, `.getStart` (bounds and starts), `.stationarityAparch`, `gsFit` (names, shape length) and `gsSim`. Adding a family means editing six files.
* Three APARCH filters (`.filterAparch`, `.filterAparchForLoop`, `.filterGarch11Fit`) and dummy order tricks (`m = 1, a = 0` for no AR; `q = 1, beta = 0` for ARCH) instead of handling zero orders directly.
* The Hessian fallback block is copied twice in `gsFit` (lines 383-424); input validation is repeated per argument in `set_params`.
* Output renaming and `outindex` subsetting (`fit.R:476-490`) duplicate the unpacking logic.
* `gat.fit` (`dist_gat.R`) is a second, independent solnp fitting path; `.armaDist` is a separate pure ARMA likelihood.
* `sigma = abs(h)^(1/delta)` (`filter.R:151,236`) silently hides negative `h` instead of penalizing it.
* The distribution object design (section 5) removes most of this: each family is declared once and every other function reads from it.

Open TODO items carried over: GAT stationarity region (likely explained by `delta >= nu * d`, see section 3), `gamma = 0` issue in APARCH estimation, uninformative Hessian error.

Reference material is in `dev/` (excluded from the build): old R files recovered from `fccc2dc` in `dev/legacy/R`, old Rd in `dev/legacy/man`, old test cases in `dev/legacy/test_cases`, the handwritten filtering notes, the GAT benchmark, and 2025/2026 exploratory scripts (fast stable via `libstable4u`, ABEV3 fits compared with rugarch).

## 3. Statistical points the design must respect

* **Innovations are location scale, not standardized.** Stable has infinite variance; stable with `beta != 0` in S0 and GEV (location 0, scale 1, mean `(Gamma(1 - xi) - 1) / xi`) have nonzero mean; GAT has infinite variance when `nu * d <= 2`. The likelihood `sum(log f(z / sigma) - log(sigma))` is valid for all three, but `mu` is a location, not the conditional mean, and `sigma_t` is a conditional scale, not a standard deviation. Document this in `gs_fit` and in `sigma()`. Decision for the GEV mean: keep it as a location (simpler, matches the thesis), and report the implied conditional mean `mu + sigma_t E[z]` as a derived quantity.
* **Mean equation in intercept form, as in Wuertz et al. (2006) and fGarch:** `x_t = mu + sum a_i x_{t-i} + sum b_j e_{t-j} + e_t`, so `mu` is the intercept and results can be compared with fGarch directly when debugging. The old `.filterArma` used the demeaned form `x_t - mu = sum a_i (x_{t-i} - mu) + ...`, which does not match its own header comment; this is fixed in `filter_arma()`. Papers needing the unconditional location recompute it as `mu / (1 - sum a_i)`.
* **Moment existence bounds `delta`.** `E(|z| - gamma z)^delta` is finite only for `delta < alpha` (stable), `delta < 1 / xi` when `xi > 0` (GEV), `delta < nu * d` (GAT). Each distribution object carries a `max_power(par)` function. Because this couples `delta` with distribution parameters, box bounds cannot enforce it: it is enforced for every algorithm through the penalty in `neg_loglik` (with a margin). This probably explains the exploding GAT simulations.
* **Stationarity is reported, not imposed.** Constrained fitting (`sqp.restriction`) is dropped to future work: an estimate pinned at the constraint boundary (persistence equal to `1 - tol`) is not an interior maximum, its standard errors are invalid, and it hides the real message of the data (near integrated volatility). Instead the estimation uses plain parameter intervals, and after fitting `gs_fit` computes the stationarity measure once and warns when it is above or close to 1.
* **What "GARCH" means per family.** In GARCH mode `delta` is fixed by the family, not always 2: the code uses `delta = 1` for stable (absolute value GARCH, Taylor/Schwert, `fit.R:220-227`) and 2 otherwise. Each distribution object declares `default_delta`; the fixed value must satisfy `max_power` (GAT needs `nu * d > 2`, GEV needs `xi < 0.5` when `delta = 2`). Document that with stable, "GARCH" is a scale recursion, not a variance recursion.
* **Stable S0 moments.** The closed forms in the code (Mittnik, Paolella, Rachev 2002) are for S1. S0 is S1 shifted by `beta * tan(pi * alpha / 2)`, and `E(|z + c| - gamma (z + c))^delta` has no closed form for `c != 0`, so **the closed form covers S0 only at `beta = 0`**. For `beta != 0` the moment is numerical integration of the `libstable4u` density. Since the moment is now computed only once per fit (stationarity report) and once per forecast, plain numerical integration is fast enough; no cached grid is needed.
* **GEV support.** The likelihood must return `+Inf` (penalty) when any `1 + xi z_t / sigma_t <= 0`.
* **APARCH identifiability.** `gamma_i` is not identified when `alpha_i` goes to 0, and `delta` trades off against `alpha`. Flag these cases instead of reporting meaningless standard errors; consider profiling over `delta`.
* **Forecasting with infinite variance.** With stable innovations there is no forecast variance; forecasts are the conditional location/mean, the conditional scale `sigma_{t+k}` and predictive quantiles. The scale recursion for `k >= 2` uses `kappa = E(|z| - gamma z)^delta` (the same moment as the stationarity report). The one step predictive distribution is exact (a shifted and scaled stable). For `k >= 2` the predictive distribution is not stable (mixture over random scales), so quantiles come from simulating paths from the fitted state. The point forecast of the mean needs `E[z]` (for S0, `-beta tan(pi alpha / 2)`, finite because `alpha > 1`) and `E[sigma_{t+k}]`, which is exact from the recursion only when `delta = 1` (the stable GARCH default); for other `delta` it is taken from the simulated paths.

## 4. What is worth doing and what is not

Worth doing:

* **Stable S0, GEV and GAT innovations in ARMA-APARCH.** fGarch (norm, snorm, ged, sged, std, sstd, snig), rugarch (adds nig, ghyp, jsu, ghst), tsgarch, MSGARCH, GAS and bayesGARCH offer none of the three, and rugarch has no user density hook in `ugarchfit`. This is the reason the package exists.
* **Stationarity conditions** for APARCH under these innovations (`gsMomentAparch`), reported after fitting. Rare in software, especially for stable.
* **Forecasting** (scale, mean and predictive quantiles/VaR), validated first for stable innovations, since the author needs it now.
* **Pluggable innovation distribution**, sold for the machinery that comes with it (see section 1).
* **Fast stable likelihood via `libstable4u`**, since the stable density dominates fitting time. Verified with a benchmark, not assumed.

Not worth doing:

* **Normal, t, skew t, GED** as built in options (fGarch and rugarch cover them).
* **S1 and S2 stable parametrizations** and the optional Nolan `stable` package.
* **"Matrix filtering" as a selling point.** The handwritten notes vectorize the recursion with lagged vectors and `stats::filter(..., "recursive")`, which is the same method as fGarch's default `llh = "filter"` path (Wuertz, Chalabi, Luksan 2006). Keep it as the internal implementation, do not advertise it as new. A C++ recursion is out of scope.
* **Formula parsing and `gsSelect`.** The list spec replaces the formula; model selection is a short user loop. Can come back later.
* **S4 classes and show methods.** Replaced by S3 methods.
* **Stationarity constrained fitting** (`sqp.restriction`). Future work, see section 3.

## 5. Target design

Public API (prefix `gs_`, decided):

```r
spec  <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = FALSE,
                 include_mean = TRUE, dist = gs_stable())
model <- gs_model(spec, mean = 0, ar = 0.1, omega = 0.05,
                  alpha = 0.1, beta = 0.85, dist_par = c(alpha = 1.8, beta = 0))
x     <- gs_sim(model, n = 2000, burnin = 1000, seed = 1)
fit   <- gs_fit(x$y, spec, algorithm = c("sqp", "nlminb"),
                start = NULL)   # optional start values, e.g. coef(previous_fit) in rolling windows
coef(fit); vcov(fit); logLik(fit); sigma(fit); AIC(fit)
residuals(fit, type = c("raw", "standardized", "pit"))
                                # e_t; e_t / sigma_t; F(e_t / sigma_t) under the fitted innovation
gs_stationarity(fit)            # value < 1 means stationary; also printed by gs_fit
fc    <- predict(fit, n_ahead = 10, level = c(0.01, 0.05), n_sim = 10000)
                                # location/mean, scale, quantiles per horizon
```

Distribution object (one file per family, `R/dist_stable.R`, `R/dist_gev.R`, `R/dist_gat.R`):

```r
gs_dist(name,
        log_density,     # function(z, par): vector of log f(z), location 0 scale 1
        random,          # function(n, par)
        cdf, quantile,   # for PIT residual diagnostics, VaR and forecast quantiles
        par_names, lower, upper,
        start,           # function(z): start values from standardized data
        default_delta,   # delta used in GARCH mode (1 for stable, 2 otherwise)
        max_power,       # function(par): supremum of delta with finite E|z|^delta
        mean,            # function(par): E[z]; must return NA explicitly where it does not exist
        valid = NULL,    # function(par): extra non box constraints, TRUE/FALSE
        aparch_moment = NULL)   # function(gamma, delta, par); NULL means numerical integration
gs_stable(); gs_gev(); gs_gat()   # built in, created with gs_dist()
```

Contract: `z` is a numeric vector, `par` a named numeric vector; functions are vectorized in `z`. Stable start values come from McCulloch quantile estimates on standardized data instead of the fixed 1.9 and 0.

**No silent failures with user distributions.** Every field above is required except `valid` and `aparch_moment`, which have well defined fallbacks. In particular `max_power` and `mean` have no defaults: the user must state when moments exist (returning `NA` for a nonexistent mean is allowed, leaving it out is not). `gs_dist()` then runs a self check at the start values and at the bounds, and stops with an error naming the failing field when:

* `par_names`, `lower`, `upper` and `start` do not have matching lengths and names, or `start` is outside `(lower, upper)`;
* `log_density` is not finite on a grid, or `exp(log_density)` does not integrate to 1;
* `cdf` is not consistent with the density, or `quantile` is not its inverse;
* draws from `random` fail a KS test against `cdf`;
* `max_power(start) <= default_delta` (GARCH mode would have no finite moment);
* `mean(start)` disagrees with the numerical integral, when finite;
* `aparch_moment`, when given, disagrees with numerical integration.

At fit time, if the optimizer ends on a penalty region, `gs_fit` reports which condition was hit (bounds, `valid`, `max_power`, support, negative `h`) instead of returning `1e99` silently. Built in families go through the same checks in the test suite.

Internal structure, replacing the duplicated code of section 2:

* `R/params.R`: one pair `pack_params(model)` / `unpack_params(vec, spec)` with named vectors, used by likelihood, stationarity, start values, bounds and output. Removes all index arithmetic.
* `R/params.R` also owns output naming and subsetting (replaces `outindex`).
* `R/filter.R`: `filter_arma()` and `filter_aparch()` that accept zero orders directly (no dummy orders); one implementation each. Pure ARMA (`p = q = 0`) is not supported, as today, so `.armaDist` is dropped. Negative `h` is a penalty, not `abs(h)`. The initialization of the lags (`mean(|z|^delta)` today) is a documented argument.
* `R/likelihood.R`: `neg_loglik(par, data, spec)` = unpack, filter, `sum(dist$log_density(z / sigma) - log(sigma))`, penalty on invalid values (bounds, `valid`, `max_power`, GEV support, negative `h`).
* `R/start.R`: `start_values(data, spec)` and bounds, reading distribution specific parts from the object.
* `R/fit.R`: `gs_fit()` split into small steps (scale, bounds, optimize, Hessian, assemble). The data is rescaled by `sd(x)` internally (as fGarch does) and parameters back transformed, so finite difference steps are meaningful. Hessian by `numDeriv::hessian` or `optimHess` with relative steps; if any evaluation during differencing hits a penalty, the affected standard errors are NA with a warning; parameters on bounds are flagged.
* `R/stationarity.R`: moments and `gs_stationarity`, used by start values, the post fit report, simulation and forecasting.
* `R/forecast.R`: `predict.gs_fit()`: ARMA mean recursion, APARCH scale recursion with `kappa`, exact one step quantiles, simulated paths from the fitted state (reusing the `gs_sim` recursion) for multi step quantiles.
* `R/sim.R`, `R/methods.R` (S3), `R/validate.R` (shared argument checks).
* `gat.fit` (iid GAT fitting) is either removed or rewritten as a thin call on the same likelihood with zero orders for the variance, so there is one fitting path.

Dependencies: `Imports: Rsolnp, libstable4u (>= current version), numDeriv, stats`; `Suggests: testthat, fGarch, stabledist` (the latter two only as test oracles). GEV density and random numbers written directly.

Licence: code derived from fGarch keeps the package GPL (>= 2); keep Wuertz as contributor and copyright holder.

## 6. Phases

Each phase is one short lived branch merged into `master` with passing `devtools::check()` and GitHub Actions.

0. **chore/package-hygiene** (this branch). Done: restore lost work, move legacy material to `dev/`, this plan. Remaining: fix DESCRIPTION (dependencies, R version, date, URL), remove load time side effect and `print()` calls, stop ignoring `man/`, modernize `inst/CITATION` (`bibentry`), move `inst/doc/bibliography.txt` out of the reserved vignette folder, unexport the broken `gsSim` (or `\dontrun` its example), add GitHub Actions `R-CMD-check`. Goal: `devtools::check()` passes with the fit still broken; old tests marked `skip()` until replaced. Nothing works for users yet at this point, which is expected.
1. **feat/dist-objects.** `gs_dist` with its self checks, `gs_stable`, `gs_gev`, `gs_gat`, normal as example, **plus the moment functions** (closed forms where valid, numerical integration bounded by `max_power`), since start values, the stationarity report, simulation and forecasting all need them. Tests: densities integrate to 1; random generators match densities (KS); `libstable4u` against `stabledist` in the tails (`|x|` up to 1e4) and for alpha near 1 and 2; GEV support; `max_power`; S0 moment equals S1 closed form at `beta = 0`; incomplete or inconsistent user distributions are rejected with a clear error.
2. **feat/spec-likelihood.** `gs_spec`, `gs_model`, pack and unpack, filters, `neg_loglik`. Tests against fGarch on `dem2gbp` (normal GARCH(1,1)): first compare the log likelihood **at fGarch's own estimates** with matching initialization and scaling; then run plain `optim` started from fGarch's estimates and compare with relative tolerance about 1e-3.
3. **feat/fit.** `start_values`, `gs_fit`, scaling, Hessian strategy, S3 methods, post fit stationarity report and penalty diagnostics. Test: fGarch cross check for normal GARCH(1,1) and ARMA(1,1)-APARCH(1,1).
4. **feat/sim.** `gs_sim` with burn in, seed and stationarity check. Tests: simulate and fit round trips for stable, GEV, GAT over a small (m,n,p,q) grid (skip on CRAN); edge cases (stable alpha near 2 vs normal, GAT large nu); investigation of the GAT stationarity region with `max_power`.
5. **feat/forecast.** `predict.gs_fit()`, validated for **stable** in this release; GEV and GAT run through the same generic code but are marked experimental until validated. Validation for stable: (a) the recursion reproduces fGarch `predict` for the normal example (checks the ARMA and scale recursions); (b) one step quantiles equal the analytic stable quantiles; (c) multi step point forecasts from the recursion agree with the simulated paths within Monte Carlo error (exact case `delta = 1`); (d) on long simulated stable ARMA-GARCH series, rolling one step VaR at 1% and 5% passes the Kupiec coverage test.
6. **docs/release.** README, one vignette (custom distribution with the normal example, plus a stable forecasting example), `NEWS.md` replacing `ChangeLog`, examples and tests use `dem2gbp` from fGarch (Suggests), no data set shipped, benchmark script (`dev/bench/`, libstable4u vs stabledist inside a full fit), version 2.0.0, CRAN submission. If 1.1 is still on CRAN, add `.Defunct()` stubs for `gsFit`, `gsSim`, `gsSpec`, `gsSelect` pointing to the new names.

## 7. Branches

All branches are kept (`develop`, `feature/repository-organization`, `feature/organize-tests`, `doc/readme_formatting`, `package-hygiene`) as a record of the modification history. Work model: `master` always passes check; each phase gets its own `feat/*`, `fix/*`, `docs/*` or `chore/*` branch, merged with `--no-ff` so the history of each change stays visible.

## 8. Decisions (2026-09-28)

1. Prefix `gs_`. Old exported names break; version 2.0.0.
2. `mu` is a location parameter for every family; the conditional mean is a derived quantity.
3. No stationarity constrained fitting for now; stationarity is reported after fitting (section 3). Constrained fitting is future work.
4. Forecasting enters now (Phase 5), validated for stable; GEV and GAT experimental until validated.
5. No data set shipped; tests and examples use `dem2gbp` from fGarch (Suggests).
6. Mean equation in intercept form (Wuertz et al. 2006, fGarch), see section 3.
7. Multi step forecasting stays in Phase 5 as planned; one step is exact, multi step by simulated paths.
8. Needs from the AR-APARCH plus vine copula VaR paper (Klin, Fiorucci, Otiniano, Maluf, Sousa) are covered: stable AR(1)-APARCH(1,1) fit with `delta` estimated, warm starts for rolling windows (`start`), one step `mu_{t+1}` and `sigma_{t+1}` from `predict`, standardized and PIT residuals for the copula step, `quantile` of the fitted innovation to map copula scenarios back to returns.

## 9. Review log

* Round 1: confirmed the bug list and the `stats::filter` equivalence with fGarch; added moment existence (`max_power`), S0 vs S1 moments, GEV location and support, Hessian strategy, fGarch cross check before simulation, libstable4u tail tests, API stubs, benchmark, licence.
* Round 2: moved moments into Phase 1 (start values, constrained fitting and simulation depend on them); limited the S0 closed form to `beta = 0`; `max_power` enforced by penalty for every algorithm and `default_delta` per family; internal rescaling and penalty aware Hessian; pinned the fGarch comparison (initialization, scaling, tolerance); extra `gs_dist` fields (`cdf`, `quantile`, `mean`, `valid`, `start` as function); hidden release work (CITATION, data, NEWS); extra duplication (`gat.fit`, `.armaDist`, `outindex`, `abs(h)`).
* Author decisions (2026-09-28): constrained fitting dropped to future work; strict self checks for user distributions; forecasting added as Phase 5, validated for stable; all branches kept; open questions closed (section 8).
