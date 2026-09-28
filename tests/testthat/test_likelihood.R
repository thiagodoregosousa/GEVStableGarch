# fGarch fits used as oracle. Starting the recursion from fGarch's own first h,
# the filter and likelihood must reproduce fGarch's log likelihood exactly
fgarch_fit <- function(formula)
{
  data("dem2gbp", package = "fGarch", envir = environment())
  x <- dem2gbp[, 1]
  list(x = x, fit = suppressWarnings(fGarch::garchFit(formula, data = x, trace = FALSE)))
}

test_that("GARCH(1,1) normal likelihood equals fGarch at fGarch's estimates", {
  skip_if_not_installed("fGarch")
  r <- fgarch_fit(~garch(1, 1))
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  par <- fGarch::coef(r$fit)[c("mu", "omega", "alpha1", "beta1")]
  expect_equal(.neg_loglik(par, r$x, spec, h_init = r$fit@h.t[1]), unname(r$fit@fit$llh), tolerance = 1e-8)
})

test_that("ARMA(1,1)-APARCH(1,1) normal likelihood and scale equal fGarch", {
  skip_if_not_installed("fGarch")
  r <- fgarch_fit(~arma(1, 1) + aparch(1, 1))
  spec <- gs_spec(arma = c(1, 1), garch = c(1, 1), aparch = TRUE, dist = .gs_norm())
  par <- fGarch::coef(r$fit)[.par_names(spec)]
  expect_equal(.neg_loglik(par, r$x, spec, h_init = r$fit@h.t[1]), unname(r$fit@fit$llh), tolerance = 1e-8)
  f <- .filter_model(r$x, .unpack(par, spec), h_init = r$fit@h.t[1])
  expect_equal(f$sigma, r$fit@sigma.t, tolerance = 1e-8)
})

test_that("maximizing our likelihood from fGarch's estimates stays there (GARCH(1,1))", {
  skip_if_not_installed("fGarch")
  r <- fgarch_fit(~garch(1, 1))
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  par <- fGarch::coef(r$fit)[c("mu", "omega", "alpha1", "beta1")]
  h0 <- r$fit@h.t[1]
  o <- optim(par, .neg_loglik, x = r$x, spec = spec, h_init = h0,
             control = list(maxit = 5000, reltol = 1e-12))
  expect_equal(o$par[-1], par[-1], tolerance = 5e-3)
  expect_equal(o$value, unname(r$fit@fit$llh), tolerance = 1e-6)
})

test_that("penalties report the violated condition", {
  x <- c(0.1, -0.2, 0.3, -0.1, 0.05)
  spec <- gs_spec(garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  diag <- new.env()
  par <- c(mu = 0, omega = 0.1, alpha1 = 0.1, gamma1 = 0, beta1 = 0.8, delta = 1.8,
           stable_alpha = 1.7, stable_beta = 0)
  expect_equal(.neg_loglik(par, x, spec, diag = diag), .PENALTY)
  expect_match(diag$reason, "max_power")

  spec_gev <- gs_spec(garch = c(1, 1), dist = gs_gev())
  par_gev <- c(mu = 0, omega = 0.001, alpha1 = 0.01, beta1 = 0.5, xi = 0.45)
  expect_equal(.neg_loglik(par_gev, c(x, -50), spec_gev, diag = diag), .PENALTY)
  expect_match(diag$reason, "support")
})

test_that("default initialization is scale equivariant", {
  set.seed(1)
  x <- rnorm(300)
  spec <- gs_spec(garch = c(1, 1), aparch = TRUE, dist = .gs_norm())
  par <- c(mu = 0.1, omega = 0.1, alpha1 = 0.1, gamma1 = 0.2, beta1 = 0.8, delta = 1.4)
  s <- 3
  par_s <- par; par_s["mu"] <- par["mu"] * s; par_s["omega"] <- par["omega"] * s^par["delta"]
  # rescaling data and parameters shifts the log likelihood by N log(s)
  expect_equal(.neg_loglik(par_s, s * x, spec), .neg_loglik(par, x, spec) + length(x) * log(s),
               tolerance = 1e-10)
})
