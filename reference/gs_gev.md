# GEV innovations

Generalized extreme value distribution with location 0, scale 1 and
shape `xi`, \\F(z) = \exp\\-(1 + \xi z)^{-1/\xi}\\\\ on \\1 + \xi z \>
0\\ (Gumbel when \\\xi = 0\\). The shape is restricted to (-0.5, 0.5):
above -0.5 for regular likelihood properties, below 0.5 for a finite
variance, which the GARCH mode power \\\delta = 2\\ requires.

## Usage

``` r
gs_gev(lower = c(xi = -0.49), upper = c(xi = 0.49), start = c(xi = 0.01))
```

## Arguments

- lower, upper, start:

  Optional named vectors replacing the default bounds and start value.

## Value

A `gs_dist` object.

## Details

The innovations are not centered: \\E\[z\] = (\Gamma(1-\xi) - 1)/\xi\\,
so the model parameter `mu` is a location, not the conditional mean.

## Examples

``` r
d <- gs_gev()
d$mean(c(xi = 0.1))
#> [1] 0.686287
```
