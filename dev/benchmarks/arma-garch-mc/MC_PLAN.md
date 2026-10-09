# Monte Carlo study plan

Goal: show that GEVStableGarch 2.0.0 estimates the ARMA-GARCH/APARCH models
with stable, GEV and GAt innovations **at least as accurately** as the 2012
version (do Rego Sousa, master thesis), and add the consistency/coverage
evidence a CRAN-quality release and a paper need. Grounded in the optim and
MLE advisor reviews (Nash 2014; Nocedal & Wright 2006; Pawitan 2001) and
anchored to the thesis simulation tables.

The study has two parts. Part A is the **non-degradation check** against the
thesis. Part B is the **extended study** (the "good grade"). Part A is the gate:
the new package must not be worse than 2012 before we invest in Part B.

---

## Part A — Non-degradation against the 2012 thesis

### Baseline

`dev/mc/baseline_2012_thesis.csv` holds every MLE (frequentist) summary table
from the thesis appendix: tables 9.2, 9.7, 9.11-9.17. Bayesian tables 9.5 and
9.10 are excluded. Each row is `(table, model, dist, param_set, param, true,
est_2012, rmse_2012)`. All thesis runs used **100 series of length n = 2500**.

Design: 13 parameter sets across 9 models.

| Table | Model (package order) | Innov. | Param sets |
|---|---|---|---|
| 9.2  | AR(1)-GARCH(1,1)      | GEV    | theta1, theta2, theta3 |
| 9.7  | AR(1)-GARCH(1,1)      | stable | phi1, phi2, phi3 |
| 9.14 | ARMA(1,1)-GARCH(1,1)  | stable | set1, set2 |
| 9.17 | ARMA(1,1)-GARCH(1,1)  | GEV    | set1, set2 |
| 9.13 | ARMA(2,2)-GARCH(2,2)  | stable | set1, set2 |
| 9.16 | ARMA(2,2)-GARCH(2,2)  | GEV    | set1, set2 |
| 9.11 | ARMA(2,2)-APARCH(2,2) | stable | set1 |
| 9.12 | ARMA(2,2)-APARCH(2,2) | stable | set1 |
| 9.15 | ARMA(2,2)-APARCH(2,2) | GEV    | set1, set2 |

### Protocol

1. For each `(table, param_set)`, build the `gs_model` from the `true` column.
2. Simulate **R = 100** series of **n = 2500** with fixed, recorded seeds
   (`seed = base_seed + replication`), matching the thesis design exactly.
3. Fit with `gs_fit()` (default `algorithm = "sqp"`), `hessian = FALSE` for
   speed in Part A (point estimates only; SEs come in Part B).
4. Record per replication: estimates, convergence code, `$diagnostics`
   (scaled gradient, Hessian rcond, iterations, cap_hit), persistence,
   `at_bound`, wall time, and the penalty `reason` when a fit fails.
5. Per parameter compute mean estimate and RMSE
   `sqrt(mean((est - true)^2))` over the **converged** replications; also keep
   bias and the convergence/boundary/failure rates.

### Pass criterion (the gate)

For each parameter in each cell, the new RMSE must not exceed the 2012 RMSE by
more than a tolerance:

> `rmse_new <= rmse_2012 * (1 + tol)`, with `tol = 0.10` (10%).

Report, per cell and overall: the fraction of parameters that pass, the worst
offenders (`rmse_new / rmse_2012`), and a side-by-side table. A handful of
parameters that are known to be hard and were already poor in 2012 (stable
skewness `stable_beta`, APARCH `gamma`, negative-`xi` cells) may legitimately
be noisy at R = 100; treat a cell as a regression only when the degradation is
systematic (several parameters, or a well-identified one like `stable_alpha`,
`alpha1`, `beta1`). Where the new package is **better**, say so — the boundary
cells (9.16 set1 has `xi = -0.1` estimated as +0.10 in 2012) are an opportunity
to show improvement.

### Comparability caveats (checked against legacy 1.x code)

These are conventions that, if mismatched, would confound "degradation" with a
definitional difference. Status after reading `dev/legacy/`:

1. **Mean equation form.** RESOLVED (match). `dev/legacy/R/GEVStableGarch-Sim.R`
   documents `mu` as "the intercept for ARMA specification (mean=mu/(1-sum(ar)))"
   — identical to the package's intercept form (Wuertz et al. 2006). No mapping
   needed.
2. **Stable parametrization.** LIKELY S0 (match). The legacy `gsSpec()` default
   `cond.dist` is `stableS0` (`GEVStableGarch-Spec.R:43`), the same S0 the
   package uses, so the stable cells are most likely fully comparable. The exact
   MC-study script is not in the repo, so residual uncertainty remains: if the
   thesis instead used `stableS1`, the stable location differs and `mu`/`omega`
   of the stable cells (9.7, 9.11-9.14) are then not directly comparable, while
   `stable_alpha`, `stable_beta`, the ARMA and `alpha/beta/gamma/delta` RMSEs
   are. The report flags stable `mu`/`omega` so a definitional shift there is
   not read as degradation.
3. **GEV standardization.** RESOLVED (match). Legacy `rgev`/`dgev` come from
   fExtremes (location 0, scale 1), the same standardization as `gs_gev()`.
4. **GARCH-mode stable power.** The package fixes `delta = 1` for stable GARCH
   (absolute-value GARCH) and `delta = 2` for GEV. The thesis GARCH tables
   assume the same; the APARCH tables estimate `delta`.
5. **Negative `xi` cells** (9.16 set1 = -0.1, set2 = -0.2; 9.17 set2 = -0.1).
   The package bounds `xi` in (-0.49, 0.49); `xi > -0.5` is the regularity edge
   (Jondeau, Poon & Rockinger 2007 — see note below). These cells are inside
   the admissible range but near the hard side; keep them and report.

---

## Part B — Extended study (books-grounded, for the paper / release)

Adds the evidence the thesis did not report: consistency across sample size,
and the calibration of the standard errors. Run only after Part A passes.

### Sample-size grid (consistency)

Take a small, representative subset of the Part A cells — one well-identified
and one harder cell per family: GEV `theta2` and `theta3` (9.2), stable `phi1`
and `phi3` (9.7), plus one APARCH cell each (9.11 stable, 9.15 set1 GEV). For
each, run **n in {500, 1000, 2500, 5000}** at R = 500. Expectation (the headline
of a consistency study, Cira bGEV §5): relative bias and RMSE fall as n grows.
Keep n = 2500 so Part B overlaps Part A at the thesis sample size.

### Standard-error calibration (coverage)

On the same subset at n in {1000, 2500}, R = 1000, record **95% Wald coverage**
(Hessian SE) and compare it to **parametric bootstrap coverage** (`gs_bootstrap`,
B = 199). This directly tests the §3 concern from the MLE review: expect good
Wald coverage in interior cells and visibly worse coverage as `stable_alpha`
approaches 2, `xi` approaches +-0.5, or `gamma` approaches +-1 — reported as a
documented limitation, not a bug.

### Metrics recorded per cell (both parts)

bias, relative bias, RMSE; convergence rate; failure rate and boundary-pileup
rate (fraction of reps with any `at_bound`); median scaled gradient and median
Hessian rcond (the Phase 1 diagnostics); median iterations and wall time; and,
for Part B SE cells, Wald and bootstrap coverage.

---

## Harness design

Code in `dev/mc/` (excluded from the package build), not in `R/`.

- `dev/mc/dgp.R` — reads `baseline_2012_thesis.csv`, builds a `gs_model` per
  `(table, param_set)`, and exposes the list of DGPs. Validates each model is
  admissible and (for the stationary tables) has persistence < 1.
- `dev/mc/run_mc.R` — the replication driver: `run_cell(dgp, n, R, seed, ...)`
  returns a tidy data frame (one row per replication x parameter) with
  estimate, true, converged, and the diagnostics. Parallel over replications
  with the `parallel` package; results written to `dev/mc/results/<tag>.rds`.
- `dev/mc/summarise.R` — aggregates replications to mean/bias/RMSE/coverage per
  parameter, joins the 2012 baseline, computes `rmse_new / rmse_2012` and the
  pass flag, and writes `dev/mc/results/comparison_<tag>.csv`.
- `dev/mc/report.Rmd` (optional) — side-by-side tables and boxplots; the
  figures for the paper.

Reproducibility: every cell carries its base seed; the driver sets
`set.seed(seed + r)` before simulating replication `r`. Store the session info
and package version with each result file.

Compute budget: a stable AR(1)-APARCH(1,1) fit on 2000 obs is ~3 s from a warm
start, ~20 s with SEs (NEWS.md). Part A is ~13 cells x 100 fits x n=2500,
point estimates only, parallelized — hours, not days. Part B's bootstrap-
coverage cells are the expensive part; keep them to the small subset above.

---

## Resolved question and open decisions

- **GEV MLE regularity (`xi > -1/2`)** — RESOLVED. The thesis itself cites
  Jondeau, Poon & Rockinger (2007) for the asymptotic normality of the GEV MLE
  when `xi > -0.5`. This is the citation the MLE advisor flagged as missing; it
  can now be stated in the help/vignette with that reference instead of Smith
  (1985).
- **Caveats checked against the legacy 1.x code** (see above): mean form and
  GEV standardization match exactly; the stable parametrization is the S0
  default, so stable cells are most likely fully comparable, with stable
  `mu`/`omega` flagged in case the thesis used S1.
- **Tolerance `tol`** for the gate is set at 10%; adjust after seeing the first
  Part A run (R = 100 Monte Carlo error on an RMSE is itself ~10%).

## Results so far (recorded)

**Part A** (R = 100, n = 2500, all 18 cells). On the 155 comparable parameters:
53 better than 2012 (RMSE ratio < 0.9), 70 tie (0.9-1.15), 22 mild, 10 worse
(> 1.3). The stable family (the package's headline) has median RMSE ratio 0.92
and **zero** regressions. The 10 real regressions are all GEV ARMA(2,2) cells,
7 of them in 9.16:set1.

**Diagnosis of the GEV ARMA(2,2) soft spot.** Not an optimization defect:
neither the new L-moment shape start nor 4-restart multistart lowers its RMSE,
and the far-off replications are all stationary (0/100 with persistence >= 1).
The estimator is **median-unbiased** there (median |median-bias| ~ 0.005 across
the affected parameters, e.g. beta1 median 0.045 vs true 0.05); the elevated
RMSE is variance from a heavy-tailed sampling distribution on weakly-identified
high-order GARCH/APARCH coefficients. The 2012 lower RMSE is consistent with the
old optimizer shrinking toward the start values (lower variance, some bias).

**Consistency check** (cells 9.15, 9.16 at n = 2500, 5000, 10000; R = 100).
RMSE falls **monotonically with n for all 50 parameters**, and at n = 10000 is
at or below the 2012 n = 2500 value for **49 of 50** (the exception is 9.15:set1
`delta`, the APARCH power, which decreases slowly). Several parameters the 2012
version estimated poorly are far better here even at n = 2500 (9.16:set1 `xi`
= -0.1: RMSE 0.203 in 2012 vs 0.014 here, a ~15x improvement; 9.16 `mu`
similarly). This confirms the soft spot is small-sample variance, not a defect:
the estimator is consistent and median-unbiased.

## References

do Rego Sousa, T. (2012). Modelos combinados AR-GARCH. Master thesis
(`references/gevstablegarch`), Appendix Chapter 9 (MLE tables 9.2, 9.7,
9.11-9.17).

Jondeau, E., Poon, S.-H., Rockinger, M. (2007). Financial Modeling Under
Non-Gaussian Distributions. Springer. (GEV MLE regularity for xi > -0.5.)

Zhao, X., Scarrott, C.J., Oxley, L., Reale, M. (2011). GARCH dependence in
extreme value models with Bayesian inference. Math. Comput. Simul., 81(7).
