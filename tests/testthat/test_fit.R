dem2gbp_x <- function()
{
  data("dem2gbp", package = "fGarch", envir = environment())
  dem2gbp[, 1]
}

test_that("normal GARCH(1,1) estimates and standard errors match fGarch", {
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  f <- fGarch::garchFit(~garch(1, 1), data = x, trace = FALSE)
  g <- gs_fit(x, gs_spec(garch = c(1, 1), dist = .gs_norm()))
  cf <- fGarch::coef(f)
  expect_equal(coef(g)[c("omega", "alpha1", "beta1")], cf[c("omega", "alpha1", "beta1")], tolerance = 5e-3)
  expect_equal(coef(g)[["mu"]], cf[["mu"]], tolerance = 5e-3)
  expect_equal(sqrt(diag(vcov(g))), f@fit$se.coef[names(coef(g))], tolerance = 0.02)
})

test_that("normal ARMA(1,1)-APARCH(1,1) reaches at least fGarch's likelihood", {
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  spec <- gs_spec(arma = c(1, 1), garch = c(1, 1), aparch = TRUE, dist = .gs_norm())
  f <- suppressWarnings(fGarch::garchFit(~arma(1, 1) + aparch(1, 1), data = x, trace = FALSE))
  g <- gs_fit(x, spec)
  # fGarch's reported llh cannot be compared directly: it starts the recursion from
  # omega + persistence * mean(z^2) (squares even when delta != 2), while we start from
  # mean(|e|^delta). For delta != 2 the two starting values differ, so the two
  # likelihoods differ slightly at the same parameters (0.07 here). Evaluating
  # fGarch's estimates with OUR likelihood puts both on the same objective:
  # our optimizer must reach at least that value, otherwise it stopped early.
  # (test_likelihood.R shows the likelihoods are identical with the same start.)
  expect_gte(g$loglik, -.neg_loglik(fGarch::coef(f)[.par_names(spec)], x, spec) - 1e-6)
  expect_equal(coef(g)[c("alpha1", "beta1")], fGarch::coef(f)[c("alpha1", "beta1")], tolerance = 0.02)
  se <- sqrt(diag(vcov(g)))
  expect_true(all(is.finite(se)))
  expect_equal(se[["delta"]], f@fit$se.coef[["delta"]], tolerance = 0.15)
})

test_that("stable AR(1)-APARCH(1,1) fit has finite standard errors and a reproducible optimum", {
  skip_on_cran()
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  g <- gs_fit(x, spec)
  expect_true(all(is.finite(sqrt(diag(vcov(g))))))
  expect_equal(g$loglik, -985.567, tolerance = 1e-5)
  expect_lt(coef(g)[["delta"]], coef(g)[["stable_alpha"]])
  # a warm start from the estimates returns the same optimum
  g2 <- gs_fit(x, spec, start = coef(g), hessian = FALSE)
  expect_equal(g2$loglik, g$loglik, tolerance = 1e-8)
})

test_that("residual types, fitted values and scale are consistent", {
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  g <- gs_fit(x, gs_spec(garch = c(1, 1), dist = gs_gev()), hessian = FALSE)
  e <- residuals(g); z <- residuals(g, type = "standardized"); u <- residuals(g, type = "pit")
  expect_equal(z, e / sigma(g))
  expect_equal(fitted(g) + e, x)
  expect_true(all(u > 0 & u < 1))
  expect_equal(u, gs_gev()$cdf(z, coef(g)["xi"]))
  expect_equal(length(sigma(g)), length(x))
})

test_that("S3 methods and printing", {
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  g <- gs_fit(x, gs_spec(garch = c(1, 1), dist = .gs_norm()))
  expect_s3_class(logLik(g), "logLik")
  expect_equal(attr(logLik(g), "df"), 4)
  expect_equal(nobs(g), length(x))
  expect_equal(AIC(g), -2 * g$loglik + 8)
  expect_output(print(g), "Persistence")
  expect_equal(g$persistence, sum(coef(g)[c("alpha1", "beta1")]))
})

test_that("invalid start values and data are reported", {
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  expect_error(gs_fit(x, spec, start = c(foo = 1)), "Unknown names")
  expect_error(gs_fit(x[1:20], spec), "at least 50")
  expect_error(gs_fit(c(x, NA), spec), "finite")
})

test_that("a fit with persistence of 1 or more warns", {
  skip_on_cran()
  skip_if_not_installed("fGarch")
  x <- dem2gbp_x()
  expect_warning(g <- gs_fit(x, gs_spec(garch = c(1, 1), dist = gs_gat()), hessian = FALSE),
                 "Persistence")
  expect_gte(g$persistence, 1)
})

test_that("stationarity of a model uses the APARCH moment of the innovations", {
  spec <- gs_spec(garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  dp <- c(stable_alpha = 1.7, stable_beta = 0)
  m <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, gamma = 0.2, beta = 0.8, delta = 1.2, dist_par = dp)
  expect_equal(gs_stationarity(m), 0.8 + 0.1 * gs_aparch_moment(gs_stable(), dp, 0.2, 1.2))
})
