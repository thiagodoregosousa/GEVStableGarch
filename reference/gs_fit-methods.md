# Methods for fitted models

Accessors for `gs_fit` objects.
[`residuals()`](https://rdrr.io/r/stats/residuals.html) returns, by
`type`: `"raw"` the residuals \\e_t\\, `"standardized"` the innovations
\\e_t / \sigma_t\\, `"pit"` their probability integral transform \\F(e_t
/ \sigma_t)\\ under the fitted innovation distribution (uniform under a
correct model; the input of copula models and of KS diagnostics).
[`sigma()`](https://rdrr.io/r/stats/sigma.html) returns the conditional
scale \\\sigma_t\\ (not a standard deviation when the innovation
variance is infinite),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) the conditional
location \\x_t - e_t\\.

## Usage

``` r
# S3 method for class 'gs_fit'
coef(object, ...)

# S3 method for class 'gs_fit'
vcov(object, ...)

# S3 method for class 'gs_fit'
logLik(object, ...)

# S3 method for class 'gs_fit'
nobs(object, ...)

# S3 method for class 'gs_fit'
residuals(object, type = c("raw", "standardized", "pit"), ...)

# S3 method for class 'gs_fit'
sigma(object, ...)

# S3 method for class 'gs_fit'
fitted(object, ...)

# S3 method for class 'gs_fit'
print(x, ...)

# S3 method for class 'gs_fit'
summary(object, ...)
```

## Arguments

- object, x:

  A `gs_fit` object.

- ...:

  Unused.

- type:

  Type of residuals.

## Value

- [`coef()`](https://rdrr.io/r/stats/coef.html): named numeric vector of
  the estimated parameters.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html): numeric matrix, the
  estimated covariance matrix of the parameters (rows and columns in the
  order of [`coef()`](https://rdrr.io/r/stats/coef.html)).

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): an object of class
  `"logLik"`, the maximized log-likelihood, with attributes `df` (number
  of parameters) and `nobs`.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html): integer, the number of
  observations used in the fit.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): numeric
  vector with one value per observation, as described by `type`.

- [`sigma()`](https://rdrr.io/r/stats/sigma.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): numeric
  vectors with one value per observation, the conditional scale and the
  conditional location.

- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html): print the
  coefficient table (estimates, standard errors, t and p values),
  information criteria, persistence and convergence diagnostics, and
  return `x` (`object`) invisibly.
