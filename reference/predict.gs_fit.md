# Forecast an ARMA-GARCH/APARCH model

Forecasts from a fitted model (`predict.gs_fit`) or from a model with
given parameters and a data history (`predict.gs_model`, useful in
rolling windows with fixed parameters).

## Usage

``` r
# S3 method for class 'gs_fit'
predict(
  object,
  n_ahead = 1,
  level = c(0.01, 0.05),
  n_sim = 10000,
  seed = NULL,
  ...
)

# S3 method for class 'gs_model'
predict(
  object,
  data,
  n_ahead = 1,
  level = c(0.01, 0.05),
  n_sim = 10000,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A `gs_fit` or `gs_model` object.

- n_ahead:

  Forecast horizon.

- level:

  Probabilities of the reported quantiles.

- n_sim:

  Number of simulated paths used when `n_ahead > 1`.

- seed:

  Optional seed for the simulated paths.

- ...:

  Unused.

- data:

  For `gs_model`, the observed history (numeric vector).

## Value

A data frame with columns `horizon`, `location` (one step only), `mean`,
`mean_method`, `sigma` and one column per `level` named `q_<level>`.

## Details

One step ahead everything is exact: the conditional location \\m\_{T+1}
= \mu + \sum a_i x\_{T+1-i} + \sum b_j e\_{T+1-j}\\, the conditional
scale \\\sigma\_{T+1}\\, and the quantiles \\m\_{T+1} + \sigma\_{T+1}
F^{-1}(p)\\ (the value at risk at level `p` is minus the `p` quantile).
For horizons \\k \ge 2\\:

- `sigma` is \\(E\\\sigma^\delta\_{T+k})^{1/\delta}\\ from the APARCH
  recursion with \\\kappa = E(\|z\| - \gamma z)^\delta\\ (as in fGarch);

- quantiles come from `n_sim` simulated paths, since the predictive
  distribution is a scale mixture and has no closed form;

- `mean` is exact when \\E\[z\] = 0\\ or \\\delta = 1\\, otherwise it is
  the average of the simulated paths (column `mean_method`).

With stable innovations there is no forecast variance; `sigma` is a
scale. `mean` is `NA` when the innovations have no mean.

## Examples

``` r
spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
model <- gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.8,
                  dist_par = c(stable_alpha = 1.8, stable_beta = 0))
x <- gs_sim(model, n = 500, seed = 1)$y
predict(model, data = x, n_ahead = 5, n_sim = 2000, seed = 1)
#>   horizon    location          mean mean_method     sigma    q_0.01    q_0.05
#> 1       1 -0.09342497 -9.342497e-02       exact 0.6321559 -2.797025 -1.676901
#> 2       2          NA -9.342497e-03       exact 0.6359273 -2.748175 -1.621558
#> 3       3          NA -9.342497e-04       exact 0.6394230 -3.059095 -1.624391
#> 4       4          NA -9.342497e-05       exact 0.6426629 -2.747105 -1.683188
#> 5       5          NA -9.342497e-06       exact 0.6456660 -3.091593 -1.647842
```
