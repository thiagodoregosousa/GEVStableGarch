# Stable innovations (S0 parametrization)

Standardized stable distribution \\S(\alpha, \beta, 1, 0; 0)\\ in
Nolan's S0 parametrization, computed with `libstable4u`. Parameters:
`stable_alpha` (tail index, in (1, 2)) and `stable_beta` (skewness, in
(-1, 1)).

## Usage

``` r
gs_stable(
  lower = c(stable_alpha = 1.01, stable_beta = -0.99),
  upper = c(stable_alpha = 1.99, stable_beta = 0.99),
  start = c(stable_alpha = 1.8, stable_beta = 0)
)
```

## Arguments

- lower, upper, start:

  Optional named vectors replacing the default bounds and start value.

## Value

A `gs_dist` object.

## Details

The variance is infinite, so \\\sigma_t\\ is a conditional scale, not a
standard deviation. In GARCH mode the power is \\\delta = 1\\ (absolute
value GARCH), since \\E\|z\|^\delta\\ is finite only for \\\delta \<
\alpha\\. The mean is \\-\beta \tan(\pi\alpha/2)\\.

## References

Nolan, J. P. (2020). Univariate Stable Distributions. Springer.

Royuela-del-Val, J., Simmross-Wattenberg, F. and Alberola-Lopez, C.
(2017). libstable: Fast, Parallel, and High-Precision Computation of
alpha-Stable Distributions in R, C/C++, and MATLAB. Journal of
Statistical Software, 78(1), 1-25.

## Examples

``` r
d <- gs_stable()
d$log_density(c(-1, 0, 1), c(stable_alpha = 1.7, stable_beta = 0.2))
#> [1] -1.552620 -1.259304 -1.560797
```
