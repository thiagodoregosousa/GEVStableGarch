# Set the parameters of an ARMA-GARCH/APARCH model

Attaches parameter values to a
[`gs_spec()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_spec.md),
checking lengths, admissible ranges and moment conditions. Parameters
can be given one by one or as a named vector `par` (for example
`coef(fit)`).

## Usage

``` r
gs_model(
  spec,
  mu = NULL,
  ar = NULL,
  ma = NULL,
  omega = NULL,
  alpha = NULL,
  gamma = NULL,
  beta = NULL,
  delta = NULL,
  dist_par = NULL,
  par = NULL
)
```

## Arguments

- spec:

  A `gs_spec` object.

- mu:

  Intercept (only when `spec$include_mean`).

- ar, ma:

  AR and MA coefficients, lengths m and n.

- omega:

  Positive constant of the scale recursion.

- alpha, beta:

  Nonnegative coefficients, lengths p and q.

- gamma:

  Asymmetry coefficients in (-1, 1), length p (APARCH only).

- delta:

  Positive power (APARCH only).

- dist_par:

  Distribution parameters, named by `spec$dist$par_names`.

- par:

  Alternatively, a named vector with all free parameters.

## Value

An object of class `gs_model`.

## Examples

``` r
spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.85,
         dist_par = c(stable_alpha = 1.8, stable_beta = 0))
#> ARMA(1,0)-GARCH(1,1) with stable innovations
#>           mu          ar1        omega       alpha1        beta1 stable_alpha 
#>         0.00         0.10         0.05         0.10         0.85         1.80 
#>  stable_beta 
#>         0.00 
```
