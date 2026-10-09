# APARCH moment of the innovation distribution

Computes \\\kappa = E(\|z\| - \gamma z)^\delta\\ for standardized
innovations \\z\\, the quantity that enters the APARCH stationarity
condition and the multi step scale forecast. A closed form is used when
the distribution provides one, otherwise numerical integration.

## Usage

``` r
gs_aparch_moment(dist, par, gamma, delta)
```

## Arguments

- dist:

  A `gs_dist` object.

- par:

  Distribution parameters, named by `dist$par_names`.

- gamma:

  Asymmetry parameter in (-1, 1).

- delta:

  Power, positive.

## Value

\\\kappa\\, or `Inf` when \\\delta\\ is not below `dist$max_power(par)`.

## References

Ding, Z., Granger, C. W. J. and Engle, R. F. (1993). A long memory
property of stock market returns and a new model. Journal of Empirical
Finance, 1, 83-106.

Mittnik, S., Paolella, M. S. and Rachev, S. T. (2002). Stationarity of
stable power-GARCH processes. Journal of Econometrics, 106, 97-107.

## Examples

``` r
gs_aparch_moment(gs_stable(), c(stable_alpha = 1.7, stable_beta = 0),
                 gamma = 0.3, delta = 1)
#> [1] 1.370884
```
