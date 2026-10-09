# Persistence of the APARCH scale recursion

Computes \\\sum_j \beta_j + \sum_i \alpha_i \kappa_i\\, with \\\kappa_i
= E(\|z\| - \gamma_i z)^\delta\\ from
[`gs_aparch_moment()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_aparch_moment.md).
A value below 1 means \\E(\sigma_t^\delta)\\ is finite (Ding, Granger
and Engle, 1993; Mittnik, Paolella and Rachev, 2002). It is `Inf` when
\\\delta\\ is not below the moment bound of the innovations.

## Usage

``` r
gs_stationarity(object)
```

## Arguments

- object:

  A `gs_model` or `gs_fit` object.

## Value

The persistence, a single number.

## Examples

``` r
spec <- gs_spec(garch = c(1, 1), dist = gs_stable())
m <- gs_model(spec, mu = 0, omega = 0.05, alpha = 0.1, beta = 0.8,
              dist_par = c(stable_alpha = 1.8, stable_beta = 0))
gs_stationarity(m)
#> [1] 0.9268715
```
