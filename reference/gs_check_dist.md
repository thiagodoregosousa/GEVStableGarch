# Self check of an innovation distribution

Verifies that the pieces of a
[`gs_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_dist.md)
object are consistent at a parameter value: names and bounds, finite and
vectorized log density, density integrating to one, `quantile` inverting
`cdf`, `cdf` matching the density, `random` matching `cdf`
(Kolmogorov-Smirnov), moment conditions (`max_power` above
`default_delta`, `mean` and `aparch_moment` against numerical
integration). It stops with an error naming the failing field.

## Usage

``` r
gs_check_dist(dist, par = dist$start, tol = 0.001)
```

## Arguments

- dist:

  A `gs_dist` object.

- par:

  Parameter value where the checks are done (default `dist$start`).

- tol:

  Tolerance for the numerical comparisons.

## Value

`dist`, invisibly.
