# Stable, GEV and GAt ARMA-APARCH models

GEVStableGarch fits, simulates and forecasts ARMA-GARCH/APARCH models
whose innovations follow a stable, GEV or GAt distribution, or any
distribution you define. This vignette walks through the workflow with
GEV innovations and then builds a user defined distribution; stable and
GAt work the same way.

## The model

In the intercept form of Wuertz et al. (2006),

``` math
x_t = \mu + \sum_{i=1}^m a_i x_{t-i} + \sum_{j=1}^n b_j e_{t-j} + e_t, \qquad e_t = \sigma_t z_t,
```
``` math
\sigma_t^\delta = \omega + \sum_{i=1}^p \alpha_i (|e_{t-i}| - \gamma_i e_{t-i})^\delta + \sum_{j=1}^q \beta_j \sigma_{t-j}^\delta,
```

with $`z_t`$ i.i.d. from the innovation distribution, standardized to
location 0 and scale 1. The families differ in which moments exist:
stable innovations have infinite variance, so $`\sigma_t`$ is a
conditional *scale* and $`E|z|^\delta`$ is finite only for
$`\delta < \alpha`$; GEV has finite variance when $`\xi < 1/2`$. In
GARCH mode (`aparch = FALSE`) the power is fixed at 1 for stable
innovations and 2 for GEV and GAt.

## GEV AR(1)-APARCH(1,1)

A specification says which model to fit; a model adds parameter values.
The same code fits stable or GAt innovations by passing
[`gs_stable()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_stable.md)
or
[`gs_gat()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gat.md)
instead of
[`gs_gev()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_gev.md).

``` r

spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_gev())
spec
#> ARMA(1,0)-APARCH(1,1) with gev innovations
#> Parameters: mu, ar1, omega, alpha1, gamma1, beta1, delta, xi
model <- gs_model(spec, mu = 0.0005, ar = 0.05, omega = 0.02, alpha = 0.12,
                  gamma = 0.5, beta = 0.6, delta = 1.0, dist_par = c(xi = 0.25))
gs_stationarity(model)
#> [1] 0.7013768
```

The persistence $`\alpha\, E(|z| - \gamma z)^\delta + \beta`$ is below
one, so the model has a stationary scale level. Simulate and fit:

``` r

x <- gs_sim(model, n = 4000, seed = 7)$y
fit <- gs_fit(x, spec)
fit
#> ARMA(1,0)-APARCH(1,1) with gev innovations
#>          Estimate Std. Error t value  Pr(>|t|)    
#> mu     -0.0012890  0.0011505 -1.1203   0.26257    
#> ar1     0.0488366  0.0093522  5.2219 1.771e-07 ***
#> omega   0.0232768  0.0100046  2.3266   0.01999 *  
#> alpha1  0.1240942  0.0181085  6.8528 7.242e-12 ***
#> gamma1  0.5150681  0.0658449  7.8224 5.181e-15 ***
#> beta1   0.6106628  0.0505859 12.0718 < 2.2e-16 ***
#> delta   0.9264252  0.1643601  5.6366 1.735e-08 ***
#> xi      0.2698621  0.0133310 20.2432 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log likelihood: 3937.9401   AIC: -7859.8802   BIC: -7809.5278   n = 4000
#> Persistence: 0.7137
#> Convergence: scaled gradient 1.0e-04, Hessian rcond 9.3e-05
```

The estimates recover the true values — note the fit starts from
data-driven defaults far from the truth (the leverage `gamma` from 0,
the power `delta` from 2, `beta` from 0.8), so this is a genuine
recovery, not a nudge from a warm start.

Three kinds of residuals are available. The raw residuals are
$`e_t = x_t - \hat\mu - \sum \hat a_i x_{t-i} - \sum \hat b_j e_{t-j}`$.
Dividing them by the fitted scale gives the standardized residuals
$`\hat z_t = e_t / \hat\sigma_t`$, the estimates of the innovations: if
the model is right they are i.i.d. from the fitted GEV law, so they
should show no remaining autocorrelation, in themselves or in their
absolute values (Ljung-Box tests).

Applying the fitted distribution function gives the probability integral
transform $`\hat u_t = F(\hat z_t)`$. Under a correct model $`\hat u_t`$
is i.i.d. uniform on (0, 1), which gives a goodness of fit test of the
innovation distribution (Kolmogorov-Smirnov below), and these are the
pseudo observations that copula models take as input.

``` r

z <- residuals(fit, type = "standardized")
Box.test(abs(z), lag = 20, type = "Ljung-Box")$p.value
#> [1] 0.190353
u <- residuals(fit, type = "pit")
ks.test(u, "punif")$p.value
#> [1] 0.9696844
```

## Forecasting and value at risk

One step ahead the forecast is exact: location, scale and quantiles of a
shifted and scaled GEV law. The value at risk at level $`p`$ is minus
the $`p`$ quantile.

``` r

predict(fit, n_ahead = 1, level = c(0.01, 0.05))
#>   horizon      location       mean mean_method      sigma      q_0.01
#> 1       1 -0.0001446314 0.05702881       exact 0.06103953 -0.07654194
#>        q_0.05
#> 1 -0.05811188
```

Beyond one step the scale follows the APARCH recursion and the quantiles
come from simulated paths.

``` r

predict(fit, n_ahead = 5, level = 0.01, n_sim = 5000, seed = 2)
#>   horizon      location       mean mean_method      sigma      q_0.01
#> 1       1 -0.0001446314 0.05702881       exact 0.06103953 -0.07654194
#> 2       2            NA 0.06164169   simulated 0.06263350 -0.07906787
#> 3       3            NA 0.06349711   simulated 0.06377301 -0.08195583
#> 4       4            NA 0.06331053   simulated 0.06458724 -0.08083717
#> 5       5            NA 0.06382599   simulated 0.06516885 -0.08273912
```

A rolling window backtest evaluates one step forecasts out of sample.
The window has a fixed length $`w`$. At each origin $`t`$ the model is
fitted to the last $`w`$ observations $`x_{t-w+1}, \dots, x_t`$ and the
one step VaR for $`t + 1`$ is computed. The window then moves forward by
one: the *observed* value $`x_{t+1}`$ (not the forecast) enters on the
right and $`x_{t-w+1}`$ leaves on the left, the model is refitted and
the next VaR computed. Comparing each VaR with the observation it was
made for gives the violations used in the Kupiec and Christoffersen
tests.

Consecutive windows share $`w - 1`$ observations, so the estimates
change little from one origin to the next. Starting each fit from the
previous estimates and skipping the Hessian makes each refit several
times faster:

``` r

w <- 1000
origins <- w:(length(x) - 1)
var_1 <- numeric(length(origins)); prev <- NULL
for (k in seq_along(origins)) {
  t <- origins[k]
  fit_t <- gs_fit(x[(t - w + 1):t], spec, start = prev, hessian = FALSE)
  prev <- coef(fit_t)
  var_1[k] <- -predict(fit_t, level = 0.01)$q_0.01    # VaR for x[t + 1]
}
violations <- x[origins + 1] < -var_1
mean(violations)                                      # should be close to 0.01
```

## A user defined innovation distribution

Any distribution can be plugged in with
[`gs_dist()`](https://thiagodoregosousa.github.io/GEVStableGarch/reference/gs_dist.md).
Every field that the machinery needs is required, including when moments
exist (`max_power`) and the mean (`mean`, which must return `NA`
explicitly if it does not exist). The constructor checks the pieces
against each other and stops with an error naming the faulty field, so a
mistake cannot silently distort a fit. Here the normal distribution:

``` r

norm <- gs_dist(
  name = "norm",
  log_density = function(z, par) dnorm(z, log = TRUE),
  random = function(n, par) rnorm(n),
  cdf = function(q, par) pnorm(q),
  quantile = function(p, par) qnorm(p),
  par_names = character(0), lower = numeric(0), upper = numeric(0), start = numeric(0),
  default_delta = 2,
  max_power = function(par) Inf,
  mean = function(par) 0)
```

A wrong declaration is caught immediately:

``` r

gs_dist(
  name = "norm",
  log_density = function(z, par) dnorm(z, log = TRUE),
  random = function(n, par) rnorm(n),
  cdf = function(q, par) pnorm(q),
  quantile = function(p, par) qnorm(p),
  par_names = character(0), lower = numeric(0), upper = numeric(0), start = numeric(0),
  default_delta = 2,
  max_power = function(par) Inf,
  mean = function(par) 0.5)
#> Error:
#> ! gs_dist 'norm': `mean` (0.5) disagrees with numerical integration (0)
```

The new distribution gets the same fitting, stationarity and
forecasting:

``` r

spec_n <- gs_spec(garch = c(1, 1), dist = norm)
m_n <- gs_model(spec_n, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.8)
fit_n <- gs_fit(gs_sim(m_n, n = 1000, seed = 3)$y, spec_n)
coef(fit_n)
#>          mu       omega      alpha1       beta1 
#> -0.03224600  0.09284506  0.13943092  0.77158255
predict(fit_n, n_ahead = 3, n_sim = 2000, seed = 4)
#>   horizon  location      mean mean_method    sigma    q_0.01    q_0.05
#> 1       1 -0.032246 -0.032246       exact 1.172140 -2.759051 -1.960245
#> 2       2        NA -0.032246       exact 1.159525 -2.697111 -1.968024
#> 3       3        NA -0.032246       exact 1.147911 -2.716736 -1.857086
```

## Relationship to other packages

This package builds on existing work. The ARMA-APARCH filter and model
structure follow fGarch (Wuertz, Chalabi and Luksan, 2006); maximum
likelihood optimization uses Rsolnp (Galanos and Ye, 2025); and the
stable density uses libstable4u (Royuela-del-Val, Simmross-Wattenberg
and Alberola-Lopez, 2017). Its own contribution is the pluggable
innovation-distribution interface and the closed-form APARCH moments for
the stable, GEV and GAt families, which extend GARCH stationarity and
forecasting (the standard fGarch-style recursion) to these heavy-tailed
innovations.

## References

Galanos, A. and Ye, Y. (2025). Rsolnp: General Non-Linear Optimization.
R package version 2.0.1. <https://CRAN.R-project.org/package=Rsolnp>

Mittnik, S., Paolella, M. S. and Rachev, S. T. (2002). Stationarity of
stable power-GARCH processes. *Journal of Econometrics*, 106, 97-107.

Paolella, M. S. (1997). Tail estimation and conditional modeling of
heteroskedastic time series. PhD thesis, University of Kiel.

Royuela-del-Val, J., Simmross-Wattenberg, F. and Alberola-Lopez, C.
(2017). libstable: Fast, Parallel, and High-Precision Computation of
alpha-Stable Distributions in R, C/C++, and MATLAB. *Journal of
Statistical Software*, 78(1), 1-25.

Wuertz, D., Chalabi, Y. and Luksan, L. (2006). Parameter estimation of
ARMA models with GARCH/APARCH errors: an R and SPlus software
implementation. *Journal of Statistical Software*.
