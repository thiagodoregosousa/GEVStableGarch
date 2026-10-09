# Simulate an ARMA-GARCH/APARCH model

Simulates the model equations of
[`gs_spec()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_spec.md)
with innovations drawn from the model's distribution. The recursion
starts at the unconditional location and scale level and the first
`burnin` values are discarded.

## Usage

``` r
gs_sim(model, n = 1000, burnin = 1000, seed = NULL)
```

## Arguments

- model:

  A `gs_model` object.

- n:

  Length of the returned series.

- burnin:

  Number of initial values discarded.

- seed:

  Optional seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Value

A data frame with columns `y` (series), `sigma` (conditional scale) and
`z` (innovations).

## Examples

``` r
spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
model <- gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.8,
                  dist_par = c(stable_alpha = 1.8, stable_beta = 0))
x <- gs_sim(model, n = 500, seed = 1)
plot(x$y, type = "l")
```
