# Fit an ARMA-GARCH/APARCH model by maximum likelihood

Estimates the model of
[`gs_spec()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_spec.md)
by conditional maximum likelihood. The data is rescaled internally to
unit standard deviation (as in fGarch) and the estimates are transformed
back. Parameters move in plain intervals; stationarity is not imposed
but reported after fitting (see
[`gs_stationarity()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_stationarity.md)),
with a warning when the persistence is 1 or more.

## Usage

``` r
gs_fit(
  data,
  spec,
  algorithm = c("sqp", "nlminb"),
  start = NULL,
  hessian = TRUE,
  restarts = 0L,
  control = list()
)
```

## Arguments

- data:

  Numeric vector, the time series (for example log returns).

- spec:

  A `gs_spec` object.

- algorithm:

  `"sqp"` (sequential quadratic programming via
  [`Rsolnp::solnp()`](https://rdrr.io/pkg/Rsolnp/man/solnp.html),
  Galanos and Ye 2025) or `"nlminb"`.

- start:

  Optional named vector of start values on the original scale, for
  example `coef(previous_fit)` in rolling windows. Missing names are
  filled with the default start values.

- hessian:

  Logical, compute standard errors. Set `FALSE` to save time when only
  point estimates are needed.

- restarts:

  Integer, number of additional optimizer runs from perturbed start
  values; the fit with the highest log likelihood is kept. A few
  restarts (for example 3 to 5) make the estimate robust to the local
  optima of multimodal surfaces, such as ARMA-GARCH models with
  heavy-tailed innovations, at a proportional increase in time. The
  default 0 keeps the single run from the data-driven start.

- control:

  List passed to the optimizer.

## Value

An object of class `gs_fit` with methods
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`sigma()`](https://rdrr.io/r/stats/sigma.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`print()`](https://rdrr.io/r/base/print.html).

## Details

Standard errors come from a numerical Hessian (central differences with
steps adapted to the curvature). When a parameter sits on a bound, or a
differentiation step leaves the parameter space, the affected standard
errors are `NA` and a warning explains why.

These Hessian standard errors are asymptotic and rely on the usual
regularity conditions. They become unreliable *near* a boundary of the
parameter space, even when not exactly on it: the stable tail index
approaching 2 (where the skewness is no longer identified), the GEV
shape approaching \\\pm 0.5\\, or an APARCH asymmetry approaching \\\pm
1\\. The GEV support depends on its shape parameter, so for GEV
innovations the usual maximum likelihood confidence intervals are not
guaranteed valid; prefer the parametric bootstrap of
[`gs_bootstrap()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_bootstrap.md)
there and in the near-boundary cases above. The fit records convergence
diagnostics in `$diagnostics` (scaled gradient, Hessian reciprocal
condition number and definiteness, optimizer iterations), which flag
these situations.

## References

Galanos, A. and Ye, Y. (2025). Rsolnp: General Non-Linear Optimization.
R package version 2.0.1. <https://CRAN.R-project.org/package=Rsolnp>

## Examples

``` r
# \donttest{
if (requireNamespace("fGarch", quietly = TRUE)) {
  data("dem2gbp", package = "fGarch")
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
  fit <- gs_fit(dem2gbp[, 1], spec)
  fit
}
#> ARMA(1,0)-GARCH(1,1) with stable innovations
#>                Estimate Std. Error t value  Pr(>|t|)    
#> mu            0.0141412  0.0074983  1.8859 0.0593073 .  
#> ar1           0.0214631  0.0217327  0.9876 0.3233504    
#> omega         0.0038210  0.0016143  2.3670 0.0179324 *  
#> alpha1        0.0770510  0.0118725  6.4899 8.591e-11 ***
#> beta1         0.8911481  0.0179740 49.5797 < 2.2e-16 ***
#> stable_alpha  1.7462593  0.0331827 52.6256 < 2.2e-16 ***
#> stable_beta  -0.3748142  0.0989306 -3.7887 0.0001515 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log likelihood: -987.6780   AIC: 1989.3560   BIC: 2028.4707   n = 1974
#> Persistence: 0.9930
#> Convergence: scaled gradient 3.5e-05, Hessian rcond 6.7e-05
# }
```
