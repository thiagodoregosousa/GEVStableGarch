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
