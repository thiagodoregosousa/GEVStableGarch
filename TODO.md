# TODO after the 2.0.0 restructuring (2026-09-28)

Details of each finding in PLAN.md section 9 and NEWS.md.

## Investigate GEV standardization vs the 2012 thesis (Monte Carlo)

- The Part B study (dev/benchmarks/arma-garch-mc) found that the GEV
  AR(1)-GARCH(1,1) cells with a large ARCH coefficient (Table 9.2 theta2
  and theta3, alpha1 = 0.80 / 0.50) are near-non-stationary in 2.0.0:
  persistence 3.2 / 3.6, because GEV is standardized to location 0 scale
  1 (so E z^2 ~ 3), not unit variance. Point estimates stay consistent
  (RMSE falls with n), but Hessian standard errors are mostly NA and the
  strict convergence gate drops with n, as expected for an
  explosive-variance model.
- The 2012 thesis treated these same cells as ordinary stationary GARCH,
  which only works if its GEV innovations were unit-variance
  standardized. So the standardization convention likely differs for
  these high-alpha GEV cells. Part A RMSE still matched, so the point
  estimates are comparable, but confirm the convention before the paper
  leans on theta2/theta3.
- Action: read the legacy 1.x GEV fit/sim path (dev/legacy) and the
  thesis to settle whether GEV was unit-variance or scale-1
  standardized; decide whether to offer a unit-variance GEV option or to
  document the scale-1 convention and drop the explosive cells from
  comparisons.

## Report the libstable4u density bug upstream

- Where: <https://github.com/swihart/libstable4u/issues> (maintainer
  Bruce Swihart, same as <https://github.com/swihart/stable>, which is
  the `stable` package, not libstable4u).

- Nobody has reported it yet. The only open issue, \#6 (Dec 2025,
  “libstable \| Wolfram Mathematica: Output difference”), reports MLE
  estimates of alpha differing from Mathematica, without replies. Our
  bug may be one cause of it: mention it there too.

- Bug: libstable4u 1.0.5 `stable_pdf(x, c(alpha, beta, 1, 0), 0)`
  returns about half the density for x within ~1e-5 of zeta =
  -beta*tan(pi*alpha/2); exact at zeta and beyond ~1e-4 (band up to
  ~1e-2 when alpha is close to 1). The cdf is fine (error ~1e-5).

- Reproducible example for the issue:

  ``` r

  a <- 1.746193; b <- -0.341291; zeta <- -b * tan(pi * a / 2)
  z <- zeta + c(-1e-4, -1e-5, 0, 1e-5, 1e-4)
  libstable4u::stable_pdf(z, c(a, b, 1, 0), 0L) / stabledist::dstable(z, a, b, 1, 0, pm = 0)
  # 1.000 0.529 1.000 0.528 1.000
  ```

- Effect: jumps of ~0.6 in the log likelihood, broken Hessians, possibly
  biased MLE. Workaround in the package: `.stable_pdf()` in
  R/dist_stable.R. Remove it once fixed upstream (the test “stable
  density is correct next to zeta” guards it).

## VaR paper (Klin, Fiorucci, Otiniano, Maluf, Sousa)

- Refit the stable AR(1)-APARCH(1,1) windows with version 2.0.0: the old
  code bounded delta below by 1 (hence delta ~1.0001 in Table 3) and the
  libstable4u bug above may have distorted those fits.
- Mean equation is now intercept form (Wuertz/fGarch): the paper’s mu
  equals the package mu.
- Write the estimation algorithm section (red TODOs in the paper) from
  the gs_fit docs.

## Package

- Push master and the phase branches to GitHub; check that GitHub
  Actions builds the vignette.
- Decide on a CRAN submission (version 2.0.0).
- Future work: stationarity constrained fitting, GAt stationarity
  region, formula interface or model selection helper if users ask.

# USING DEVTOOLS

- devtools::load_all() load your package (when run inside its folder)
- devtools::test() to run all tests
- devtools::test_active_file() run test on active file
- devtools::document() automatic generation of .Rd files in the man
  folder

# Change log

## 202601

- allows only S0 parametrization for stable (cleaner code), and replaced
  dependence from stabledist by libstable4u.
- It allows exactly only stable, gev and gat distributions. Normal just
  for testing.
- gat-tests finished. Found out that estimation of nu is more difficult.
- gat-dist comments are now in roxygen format. That should be done with
  all the other functions for automation.
- created benchmark and test_that by using the amazing general framework
  provided by simDesign package.
- decided on minimal structure of lists for model specification and
  parameters
- tested model and spec, plus documentation

# Tasks by priority

## Now

- create simulation function using model spec and parameters.
- test sim function, try to get conditions for stationarity and simulate
  only that.
- Build minimum test set that should go into pipeline (see tests
  folder):
  1.  OK garch(1,1) with normal dem2… with expected values from fGarch
      package
  2.  simulate and estimate using my package 1.1) OK test garch11 for
      all conditional distributions 1.2)
      model(m\>=0,n\>=0,p\>=1,q\>=0) + distributions (stableS0, GEV,
      GAT): test model combinations of (m,n,p,q) for m,n,q in \[0,2\]
      and p in \[1,2\]
  3.  testing extreme cases: simulate arma(1,1)-garch(1,1) with normal
      innovations. Estimate the same model using normal innovations,
      stable (alpha should be close to 2) and GAT (nu should be close to
      infinty)
  4.  test models and distribution
  5.  simulate arma(1,1)-garch(1,1), fit and compare results. Do this
      for both stable, GEV and gat
  6.  simulate arma(1,1)-garch(1,1) with normal and estimate with fGarch
      and rugarch. Same for stable conditional with alpha close to 2,
      and GAT with d=2, theta=1 should be t-student with v degrees of
      freedom. Do the oposite direction, simulate with fGarch and
      rugarch and estimate with mine.

## Latter

- variable notation to \_ instead of you.me.plot which seems method
  dispatch in R.
- Simulation with libstable4u  
- It is not clear when model with conditional GAT will be stationary, it
  is exploding for several parameters
- Printing and computing hessian mixed inside gsFit function. gsFit is
  super big function. e.g., c(mu = 0, omega = 0.4, alpha1 = 0.2, beta1 =
  0.2, skew = 0.4, shape1 = 2, shape2 = 1)
- My pdf Filtering Process for estimation (PDF DOC) is missing et = zt
  \* ht in the equation
- Error message for computing std using hessian, not informative. users
  need mathematical reasons to investigage better the output of the
  function.
- enforce stationarity using sqp.restriction algorithm must be tested
  with others datasets.
- fix gamma = 0 for estimation in aparch.

# Functions in this package

.armaDist Calculates the likelihood function for a vector of points (z)
\# according to the specified distribution.

.armaGarchDist similar to .armaDist but density for z/hh

.filterEQUATION. filter data according to a specified model by EQUATION

Fit - Fits ARMA-GARCH or ARMA-APARCH model. Returns an object of class
GEVSTABLEGARCH

getFormula - get model parameters from string, like,
‘arma(1,1)-garch(2,2)’

gsSelect (previous or related name = GSgarch.FitAIC) (returns model
parameters with minimum AIC). It needs getOrder to get list of different
parameters to fit and decide which one is the best.

getStart - starting parameters for estimation, relies on preliminary
arima fitting

otherUsefulCodes - for testing ?

Sim - model simulation

Spec - specifies model and returns an instance of class
GEVSTABLEGARCHSPEC

stationarityAparch - compute a value that should be \< 1 for the model
to be stationary.

# Naming convention:

user functions: gsFit, gsSelect, gsMomentAparch variables: cond.dist,
arma.order, garch.llh constans: TOLG, TOLSTABLE, ARMA.ORDER internal
functions: .getStart, .getFormula .

# Future modifications on package:

- Prediction methods using the results of Brockwell for stable
  prediction. See paper from Parameter “Estimation of ARMA Models with
  GARCH/APARCH Errors An R and SPlus Software Implementation”
