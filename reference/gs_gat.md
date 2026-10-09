# GAt innovations

Generalized asymmetric t distribution of Paolella (1997) with location
0, scale 1 and parameters `nu` (tail), `d` (peakedness) and `xi`
(asymmetry, 1 is symmetric), see
[`dgat()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/GAT.md).
Moments \\E\|z\|^\delta\\ exist for \\\delta \< \nu d\\; the GARCH mode
power is \\\delta = 2\\.

## Usage

``` r
gs_gat(
  lower = c(nu = 0.05, d = 0.1, xi = 0.05),
  upper = c(nu = 100, d = 50, xi = 20),
  start = c(nu = 2, d = 4, xi = 1)
)
```

## Arguments

- lower, upper, start:

  Optional named vectors replacing the default bounds and start value.

## Value

A `gs_dist` object.

## Examples

``` r
d <- gs_gat()
d$max_power(c(nu = 2, d = 4, xi = 1))
#> [1] 8
```
