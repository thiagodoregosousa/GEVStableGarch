# A new innovation distribution: bimodal GEV from the bgev package

The built-in innovations (stable, GEV, GAt) are created with the same
constructor,
[`gs_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_dist.md),
that is offered to users. This vignette shows the point of that design:
a distribution from **another CRAN package** is dropped in as the
innovation of an ARMA-GARCH/APARCH model and immediately gets the full
machinery — estimation, simulation, stationarity and forecasting — with
no changes to the package.

We use the **bimodal GEV** (BGEV) of Otiniano, Oliveira, Vila and Dias
(2025), provided by the CRAN package
[`bgev`](https://CRAN.R-project.org/package=bgev). It extends the GEV
with a second shape parameter `delta` that reshapes the body and the
tail.

## Wrapping bgev as an innovation distribution

An innovation is standardized to location 0 and scale 1; the conditional
scale `sigma_t` of the model carries the scale. So we fix bgev’s
`mu = 0`, `sigma = 1` and expose its two shape parameters. The only
bookkeeping:

- bgev’s second shape is called `delta`, which would clash with the
  APARCH power `delta`, so we rename the innovation parameters `bgev_xi`
  and `bgev_delta`;
- `max_power` must say where moments stop existing. The BGEV paper
  (Remark 1) gives, for `xi > 0`, `E|z|^k < Inf` if and only if
  `k < (delta + 1) / xi`;
- `mean` must return `E[z]`. BGEV is not centered, exactly like the
  built-in GEV, and the machinery accounts for it; we compute the mean
  numerically.

``` r

bgev_mean <- function(xi, delta)
  integrate(function(p) bgev::qbgev(p, 0, 1, xi, delta), 0, 1,
            subdivisions = 400L, rel.tol = 1e-6)$value

bgev_innov <- gs_dist(
  name = "bgev",
  log_density = function(z, par) {
    d <- bgev::dbgev(z, 0, 1, par[["bgev_xi"]], par[["bgev_delta"]])
    ifelse(d > 0, log(d), -Inf)
  },
  random   = function(n, par) bgev::rbgev(n, 0, 1, par[["bgev_xi"]], par[["bgev_delta"]]),
  cdf      = function(q, par) bgev::pbgev(q, 0, 1, par[["bgev_xi"]], par[["bgev_delta"]]),
  quantile = function(p, par) bgev::qbgev(p, 0, 1, par[["bgev_xi"]], par[["bgev_delta"]]),
  par_names = c("bgev_xi", "bgev_delta"),
  lower = c(bgev_xi = 0.02, bgev_delta = 0.5),
  upper = c(bgev_xi = 0.45, bgev_delta = 5),
  start = c(bgev_xi = 0.20, bgev_delta = 2),
  default_delta = 2,                                   # GARCH mode: variance recursion
  max_power = function(par) (par[["bgev_delta"]] + 1) / par[["bgev_xi"]],
  mean = function(par) bgev_mean(par[["bgev_xi"]], par[["bgev_delta"]]))
```

The constructor runs a self-check so that a distribution that declares
something inconsistent (densities that do not integrate to one, a wrong
`max_power`, a `quantile` that does not invert the `cdf`, …) fails
loudly rather than silently biasing the fit:

``` r

gs_check_dist(bgev_innov)
```

## Fit, simulate, forecast

From here it is an ordinary model. We simulate a GARCH(1,1) with BGEV
innovations and estimate it back.

``` r

spec  <- gs_spec(garch = c(1, 1), dist = bgev_innov)
model <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.6,
                  dist_par = c(bgev_xi = 0.1, bgev_delta = 4))
gs_stationarity(model)
#> [1] 0.6923729

x   <- gs_sim(model, n = 3000, seed = 7)$y
fit <- gs_fit(x, spec)
fit
#> GARCH(1,1) with bgev innovations
#>              Estimate Std. Error t value  Pr(>|t|)    
#> mu         -0.0025125  0.0022469 -1.1182    0.2635    
#> omega       0.1040038  0.0184086  5.6497 1.607e-08 ***
#> alpha1      0.0933730  0.0123417  7.5656 3.860e-14 ***
#> beta1       0.5950794  0.0619475  9.6062 < 2.2e-16 ***
#> bgev_xi     0.1508719  0.0204074  7.3930 1.435e-13 ***
#> bgev_delta  4.1736280  0.0867042 48.1364 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log likelihood: 553.2916   AIC: -1094.5833   BIC: -1058.5451   n = 3000
#> Persistence: 0.6822
#> Convergence: scaled gradient 1.5e-05, Hessian rcond 2.0e-04
```

The estimates sit close to the true values (`omega = 0.1`,
`alpha1 = 0.1`, `beta1 = 0.6`, `bgev_xi = 0.1`, `bgev_delta = 4`) — note
the shape parameters are recovered from the default start
`(bgev_xi = 0.2, bgev_delta = 2)`, far from the truth. Forecasting and
Value-at-Risk work unchanged; with a non-centered innovation the
conditional mean uses the declared `mean`:

``` r

predict(fit, n_ahead = 5, level = c(0.01, 0.05), n_sim = 5000, seed = 1)
#>   horizon     location      mean mean_method     sigma     q_0.01     q_0.05
#> 1       1 -0.002512544 0.1785613       exact 0.5685197 -0.6061889 -0.5722561
#> 2       2           NA 0.1786121   simulated 0.5696493 -0.6062665 -0.5739388
#> 3       3           NA 0.1782490   simulated 0.5704186 -0.6119872 -0.5738215
#> 4       4           NA 0.1700287   simulated 0.5709428 -0.6111660 -0.5780694
#> 5       5           NA 0.1762631   simulated 0.5713002 -0.6112368 -0.5761988
```

And the probability integral transform of the residuals is uniform under
a correct model, which is the usual goodness-of-fit check:

``` r

u <- residuals(fit, type = "pit")
round(quantile(u, c(0.1, 0.5, 0.9)), 3)
#>   10%   50%   90% 
#> 0.100 0.491 0.898
```

That is the whole contribution on display: the pluggable interface turns
any documented density — here a third party’s, reused rather than
reimplemented — into a fully functional heavy-tailed ARMA-GARCH/APARCH
model.

## Reference

Otiniano, C. E. G., Oliveira, B. S., Vila, R., Dias, M. (2025). *The
bimodal GEV distribution*. (CRAN package `bgev`.)
