# Innovation distribution for ARMA-GARCH/APARCH models

Builds the object that describes the innovation distribution \\z_t\\ of
the model. The built in families
[`gs_stable()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_stable.md),
[`gs_gev()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gev.md)
and
[`gs_gat()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gat.md)
are created with this constructor, and a user defined distribution gets
exactly the same fitting, simulation, stationarity and forecasting
machinery.

## Usage

``` r
gs_dist(
  name,
  log_density,
  random,
  cdf,
  quantile,
  par_names,
  lower,
  upper,
  start,
  default_delta,
  max_power,
  mean,
  valid = NULL,
  aparch_moment = NULL,
  start_fun = NULL,
  check = TRUE
)
```

## Arguments

- name:

  Character, name of the distribution.

- log_density:

  `function(z, par)` returning \\\log f(z)\\, vectorized in `z`, and
  `-Inf` outside the support.

- random:

  `function(n, par)` returning `n` draws.

- cdf:

  `function(q, par)`, distribution function.

- quantile:

  `function(p, par)`, quantile function.

- par_names:

  Character vector with the parameter names. They must not clash with
  model parameter names (`mu`, `ar1`, `omega`, `alpha1`, ...).

- lower, upper, start:

  Numeric vectors named by `par_names`: box bounds used in estimation
  and the default start value (strictly inside).

- default_delta:

  Power \\\delta\\ used when the model is a GARCH (not APARCH) model,
  e.g. 2 for a variance recursion, 1 for a scale recursion.

- max_power:

  `function(par)` returning the supremum of the powers \\\delta\\ with
  \\E\|z\|^\delta \< \infty\\ (`Inf` if all moments exist).

- mean:

  `function(par)` returning \\E\[z\]\\; it must return `NA` explicitly
  where the mean does not exist.

- valid:

  Optional `function(par)` returning `TRUE`/`FALSE` for constraints that
  are not box bounds.

- aparch_moment:

  Optional `function(gamma, delta, par)` returning \\E(\|z\| - \gamma
  z)^\delta\\ in closed form. It may return `NA` where no closed form is
  available; numerical integration is then used.

- start_fun:

  Optional `function(z)` returning data driven start values from
  standardized data.

- check:

  Logical, run
  [`gs_check_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_check_dist.md)
  at `start`.

## Value

An object of class `gs_dist`.

## Details

All functions work on the standardized innovation (location 0, scale 1).
`z` is a numeric vector and `par` a numeric vector named by `par_names`.

Every field is required except `valid`, `aparch_moment` and `start_fun`,
which have well defined fallbacks. In particular `max_power` and `mean`
have no defaults: the user must state when moments exist, so that
fitting never fails silently. Unless `check = FALSE`,
[`gs_check_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_check_dist.md)
runs a self check at `start` and stops with an error naming the failing
field.

## See also

[`gs_check_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_check_dist.md),
[`gs_stable()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_stable.md),
[`gs_gev()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gev.md),
[`gs_gat()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gat.md)

## Examples

``` r
# Normal innovations written as a user defined distribution
norm <- gs_dist(
  name = "norm",
  log_density = function(z, par) dnorm(z, log = TRUE),
  random = function(n, par) rnorm(n),
  cdf = function(q, par) pnorm(q),
  quantile = function(p, par) qnorm(p),
  par_names = character(0),
  lower = numeric(0), upper = numeric(0), start = numeric(0),
  default_delta = 2,
  max_power = function(par) Inf,
  mean = function(par) 0)
norm
#> Innovation distribution: norm 
#> No parameters
#> GARCH mode delta: 2 
```
