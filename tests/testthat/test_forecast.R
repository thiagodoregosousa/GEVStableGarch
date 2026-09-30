test_that("mean and scale forecasts equal fGarch predict (normal GARCH and APARCH)", {
  skip_if_not_installed("fGarch")
  data("dem2gbp", package = "fGarch", envir = environment())
  x <- dem2gbp[, 1]
  for (ap in c(FALSE, TRUE)) {
    form <- if (ap) ~arma(1, 1) + aparch(1, 1) else ~arma(1, 0) + garch(1, 1)
    f <- suppressWarnings(fGarch::garchFit(form, data = x, trace = FALSE))
    spec <- gs_spec(arma = if (ap) c(1, 1) else c(1, 0), garch = c(1, 1), aparch = ap, dist = .gs_norm())
    m <- gs_model(spec, par = fGarch::coef(f)[.par_names(spec)])
    ours <- predict(m, data = x, n_ahead = 10, n_sim = 500, seed = 1)
    theirs <- fGarch::predict(f, n.ahead = 10)
    expect_equal(ours$mean, theirs$meanForecast, tolerance = 1e-8)
    expect_equal(ours$sigma, theirs$standardDeviation, tolerance = 1e-8)
  }
})

test_that("one step value at risk equals fGarch's", {
  skip_if_not_installed("fGarch")
  data("dem2gbp", package = "fGarch", envir = environment())
  x <- dem2gbp[, 1]
  f <- fGarch::garchFit(~arma(1, 0) + garch(1, 1), data = x, trace = FALSE)
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = .gs_norm())
  m <- gs_model(spec, par = fGarch::coef(f)[.par_names(spec)])
  # beyond one step fGarch approximates the quantile by mean + q * meanError,
  # while the true predictive law is a scale mixture, so only one step is compared
  var_fgarch <- fGarch::predict(f, n.ahead = 1, p_loss = 0.01)$VaR
  expect_equal(-predict(m, data = x, level = 0.01)$q_0.01, var_fgarch, tolerance = 1e-8)
})

stable_model <- function(beta = 0.3)
{
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
  gs_model(spec, mu = 0.02, ar = 0.1, omega = 0.05, alpha = 0.08, beta = 0.85,
           dist_par = c(stable_alpha = 1.8, stable_beta = beta))
}

test_that("one step forecast is the next value of the filter", {
  m <- stable_model()
  x <- gs_sim(m, n = 1500, seed = 2)$y
  N <- length(x)
  fc <- predict(m, data = x[-N], level = 0.05)
  f <- .filter_model(x, .unpack(m$par, m$spec))
  expect_equal(fc$location, x[N] - f$e[N], tolerance = 1e-6)
  expect_equal(fc$sigma, f$sigma[N], tolerance = 1e-6)
})

test_that("one step quantiles are exact stable quantiles", {
  m <- stable_model()
  x <- gs_sim(m, n = 500, seed = 3)$y
  fc <- predict(m, data = x, level = c(0.01, 0.05))
  dp <- m$par[c("stable_alpha", "stable_beta")]
  expect_equal(unlist(fc[1, c("q_0.01", "q_0.05")], use.names = FALSE),
               fc$location + fc$sigma * gs_stable()$quantile(c(0.01, 0.05), dp))
  expect_equal(fc$mean, fc$location + fc$sigma * gs_stable()$mean(dp))
})

test_that("multi step stable scale and mean recursions agree with simulated paths", {
  skip_on_cran()
  m <- stable_model(beta = 0.5)
  x <- gs_sim(m, n = 500, seed = 4)$y
  fc <- predict(m, data = x, n_ahead = 5, n_sim = 2e5, seed = 5)
  expect_true(all(fc$mean_method == "exact"))  # delta = 1 for stable GARCH
  # simulate the same paths directly and compare sigma and mean
  u <- .unpack(m$par, m$spec); L <- .max_lag(u)
  f <- .filter_model(x, u); last <- length(x) - L + seq_len(L)
  set.seed(6)
  z <- matrix(gs_stable()$random(5 * 2e5, u$dist_par), nrow = 5)
  p <- .sim_paths(u, z, list(x = x[last], e = f$e[last], h = f$h[last]))
  expect_equal(fc$sigma, rowMeans(p$sigma), tolerance = 2e-3)
  expect_equal(fc$mean, rowMeans(p$x), tolerance = 0.05 * mean(fc$sigma))
})

test_that("one step value at risk has correct coverage (Kupiec) for a stable AR(1)-GARCH(1,1)", {
  m <- stable_model()
  s <- gs_sim(m, n = 20000, seed = 7)
  f <- .filter_model(s$y, .unpack(m$par, m$spec))
  loc <- s$y - f$e
  dp <- m$par[c("stable_alpha", "stable_beta")]
  kupiec <- function(p) {
    hits <- s$y < loc + f$sigma * gs_stable()$quantile(p, dp)
    k <- sum(hits); n <- length(hits); ph <- k / n
    lr <- -2 * ((n - k) * log(1 - p) + k * log(p) - (n - k) * log(1 - ph) - k * log(ph))
    stats::pchisq(lr, 1, lower.tail = FALSE)
  }
  expect_gt(kupiec(0.01), 0.01)
  expect_gt(kupiec(0.05), 0.01)
})

test_that("non centered innovations with delta != 1 use a simulated mean beyond one step", {
  spec <- gs_spec(garch = c(1, 1), dist = gs_gev())
  m <- gs_model(spec, mu = 0, omega = 0.05, alpha = 0.05, beta = 0.8, dist_par = c(xi = 0.1))
  x <- gs_sim(m, n = 300, seed = 8)$y
  fc <- predict(m, data = x, n_ahead = 3, n_sim = 1000, seed = 9)
  expect_equal(fc$mean_method, c("exact", "simulated", "simulated"))
  expect_true(all(diff(fc$q_0.01) != 0))
})

test_that("multi step forecasts warn when persistence is 1 or more", {
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  m <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.3, beta = 0.75)
  x <- suppressWarnings(gs_sim(m, n = 200, seed = 1))$y
  expect_warning(predict(m, data = x, n_ahead = 2, n_sim = 100), "Persistence")
  expect_silent(predict(m, data = x, n_ahead = 1))
})

test_that("predict works on a fitted model", {
  skip_if_not_installed("fGarch")
  data("dem2gbp", package = "fGarch", envir = environment())
  g <- gs_fit(dem2gbp[, 1], gs_spec(arma = c(1, 0), garch = c(1, 1), dist = .gs_norm()), hessian = FALSE)
  fc <- predict(g, n_ahead = 3, n_sim = 500, seed = 1)
  expect_equal(nrow(fc), 3)
  expect_named(fc, c("horizon", "location", "mean", "mean_method", "sigma", "q_0.01", "q_0.05"))
})
