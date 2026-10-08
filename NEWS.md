# GEVStableGarch 2.0.0

Complete rewrite. The version 1.1 interface (`gsFit`, `gsSpec`, `gsSim`,
`gsSelect`, S4 classes) is replaced; the package was archived on CRAN in 2020,
so no compatibility layer is provided.

## Scope

* Innovations: stable (S0 parametrization only, via `libstable4u`), GEV and
  GAt. Normal, t, skew t and GED innovations are left to fGarch and rugarch;
  S1 and S2 stable parametrizations and the `stable` package are dropped.
* Any other distribution can be plugged in with `gs_dist()`, which requires
  every field the machinery needs (including when moments exist) and checks
  them against each other (`gs_check_dist()`).

## New interface

* `gs_spec()`, `gs_model()`: specification and parameters as plain lists.
* `gs_fit()`: maximum likelihood on internally rescaled data, `sqp` (Rsolnp)
  or `nlminb`, warm starts through `start`, optional Hessian, standard errors
  set to `NA` with a warning when a parameter is on a bound or the Hessian is
  not usable, persistence warning. It records convergence diagnostics in
  `$diagnostics` (scaled gradient, Hessian reciprocal condition number and
  definiteness, optimizer iterations) and warns when the gradient is not near
  zero, the Hessian is ill conditioned, or the iteration cap is hit.
* `gs_bootstrap()`: parametric bootstrap standard errors and percentile
  confidence intervals, recommended over the Hessian errors near a boundary
  and for GEV innovations, whose support depends on the shape parameter.
* `gs_sim()`, `predict()`: simulation and forecasting (exact one step,
  APARCH recursion and simulated paths beyond).
* `gs_stationarity()`, `gs_aparch_moment()`.
* S3 methods: `coef`, `vcov`, `logLik`, `residuals` (raw, standardized, pit),
  `sigma`, `fitted`, `print`.

## Behaviour changes

* The mean equation is in intercept form, as in Wuertz et al. (2006) and
  fGarch: `mu` is the intercept, the unconditional location is
  `mu / (1 - sum(ar))`.
* The recursion starts from `mean(|e|^delta)`, which is scale equivariant.
  Starting from fGarch's own initial value, the likelihood reproduces fGarch
  exactly for GARCH and APARCH models.
* Stationarity is reported after fitting, not imposed (the old
  `sqp.restriction` algorithm is not available).
* The APARCH power `delta` is no longer bounded below by 1 for stable
  innovations; `delta < stable_alpha` is enforced instead. Version 1.1
  estimates at `delta = 1` were on that bound.

## Bug fixes

* `libstable4u` 1.0.5 returns about half the stable density for points within
  roughly 1e-5 of `zeta = -beta tan(pi alpha / 2)`. The resulting jumps of the
  log likelihood (about 0.6) broke optimization and standard errors; the
  density is now interpolated in a small window around `zeta`.
* `pgat()` kept only a few digits next to zero; it now uses the exact upper
  beta tail.
* The out-of-bounds penalty is now graded: it adds a term proportional to the
  size of the constraint violation (for example `delta` above the innovation
  moment bound, or a non-stationary AR polynomial) instead of a flat constant.
  The objective slopes back toward the feasible region, so the optimizer no
  longer stalls on a flat plateau.

## Performance

* The stable density from `libstable4u` is about 90 times faster than
  `stabledist` (2000 points: 31 ms against 2.9 s). A stable AR(1)-APARCH(1,1)
  fit on 2000 observations takes about 20 seconds with standard errors and
  3 seconds from a warm start (`dev/bench/bench_stable_density.R`).
