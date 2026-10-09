# Specify an ARMA-GARCH/APARCH model

The model, in the intercept form of Wuertz et al. (2006) and fGarch, is
\$\$x_t = \mu + \sum\_{i=1}^m a_i x\_{t-i} + \sum\_{j=1}^n b_j
e\_{t-j} + e_t, \qquad e_t = \sigma_t z_t,\$\$ \$\$\sigma_t^\delta =
\omega + \sum\_{i=1}^p \alpha_i (\|e\_{t-i}\| - \gamma_i
e\_{t-i})^\delta + \sum\_{j=1}^q \beta_j \sigma\_{t-j}^\delta,\$\$ with
\\z_t\\ i.i.d. from `dist`. In GARCH mode (`aparch = FALSE`) the
asymmetry \\\gamma_i\\ is 0 and \\\delta\\ is fixed at
`dist$default_delta` (1 for stable, 2 otherwise). `mu` is an intercept:
the unconditional location is \\\mu / (1 - \sum a_i)\\. Since the
innovations are location scale (not standardized), \\\sigma_t\\ is a
conditional scale.

## Usage

``` r
gs_spec(
  arma = c(0, 0),
  garch = c(1, 1),
  aparch = FALSE,
  include_mean = TRUE,
  dist = gs_stable()
)
```

## Arguments

- arma:

  Integer vector `c(m, n)`, ARMA orders.

- garch:

  Integer vector `c(p, q)`, number of \\\alpha\\ (p \>= 1) and \\\beta\\
  (q \>= 0) terms.

- aparch:

  Logical, estimate the asymmetry \\\gamma\\ and the power \\\delta\\.

- include_mean:

  Logical, include the intercept \\\mu\\.

- dist:

  Innovation distribution, a `gs_dist` object.

## Value

An object of class `gs_spec`.

## References

Wuertz, D., Chalabi, Y. and Luksan, L. (2006). Parameter estimation of
ARMA models with GARCH/APARCH errors: an R and SPlus software
implementation. Journal of Statistical Software.

## Examples

``` r
gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
#> ARMA(1,0)-APARCH(1,1) with stable innovations
#> Parameters: mu, ar1, omega, alpha1, gamma1, beta1, delta, stable_alpha, stable_beta 
```
