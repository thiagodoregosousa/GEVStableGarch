test_that("built in distributions pass the self check at several parameter values", {
  expect_silent(gs_check_dist(.gs_norm()))
  expect_silent(gs_check_dist(gs_stable()))
  expect_silent(gs_check_dist(gs_stable(), c(stable_alpha = 1.3, stable_beta = -0.6)))
  expect_silent(gs_check_dist(gs_gev()))
  expect_silent(gs_check_dist(gs_gev(), c(xi = -0.3)))
  expect_silent(gs_check_dist(gs_gev(), c(xi = 0.3)))
  expect_silent(gs_check_dist(gs_gat()))
  expect_silent(gs_check_dist(gs_gat(), c(nu = 1.5, d = 1.8, xi = 0.7)))
})


# A complete, valid user distribution to break one field at a time
norm_args <- function() list(
  name = "norm",
  log_density = function(z, par) dnorm(z, log = TRUE),
  random = function(n, par) rnorm(n),
  cdf = function(q, par) pnorm(q),
  quantile = function(p, par) qnorm(p),
  par_names = character(0), lower = numeric(0), upper = numeric(0), start = numeric(0),
  default_delta = 2,
  max_power = function(par) Inf,
  mean = function(par) 0)

test_that("a complete user distribution is accepted", {
  expect_s3_class(do.call(gs_dist, norm_args()), "gs_dist")
})

test_that("missing moment declarations are rejected, not defaulted", {
  a <- norm_args(); a$max_power <- NULL; a$mean <- NULL
  expect_error(do.call(gs_dist, a), "`max_power`, `mean`")
})

test_that("inconsistent user distributions fail naming the field", {
  a <- norm_args(); a$mean <- function(par) 0.5
  expect_error(do.call(gs_dist, a), "`mean`")
  a <- norm_args(); a$cdf <- function(q, par) pnorm(q, sd = 2)
  expect_error(do.call(gs_dist, a), "`quantile`|`cdf`")
  a <- norm_args(); a$random <- function(n, par) rnorm(n, mean = 1)
  expect_error(do.call(gs_dist, a), "`random`")
  a <- norm_args(); a$log_density <- function(z, par) dnorm(z, sd = 2, log = TRUE)
  expect_error(do.call(gs_dist, a), "`cdf`|`log_density`")
  a <- norm_args(); a$max_power <- function(par) 1.5
  expect_error(do.call(gs_dist, a), "`default_delta`")
  a <- norm_args(); a$max_power <- function(par) 0.8
  a$default_delta <- 0.5
  expect_error(do.call(gs_dist, a), "`mean`")
  a <- norm_args(); a$aparch_moment <- function(gamma, delta, par) 1
  expect_error(do.call(gs_dist, a), "`aparch_moment`")
})

test_that("distribution parameters cannot clash with model parameter names", {
  a <- norm_args()
  a$par_names <- "alpha1"; a$lower <- 0; a$upper <- 1; a$start <- 0.5
  expect_error(do.call(gs_dist, a), "clash")
})

test_that("stable density from libstable4u matches stabledist, also in the tails", {
  skip_if_not_installed("stabledist")
  x <- c(-1e4, -100, -10, -1, 0, 1, 10, 100, 1e4)
  for (a in c(1.05, 1.5, 1.95)) for (b in c(-0.5, 0, 0.5)) {
    par <- c(stable_alpha = a, stable_beta = b)
    ours <- exp(gs_stable()$log_density(x, par))
    ref <- stabledist::dstable(x, a, b, 1, 0, pm = 0)
    expect_equal(ours, ref, tolerance = 1e-3)
  }
})

test_that("GEV density is zero outside its support", {
  d <- gs_gev()
  expect_equal(d$log_density(-3, c(xi = 0.4)), -Inf)
  expect_equal(d$cdf(-3, c(xi = 0.4)), 0)
  expect_equal(d$log_density(3, c(xi = -0.4)), -Inf)
  expect_equal(d$cdf(3, c(xi = -0.4)), 1)
})

test_that("max_power follows the tail index of each family", {
  expect_equal(gs_stable()$max_power(c(stable_alpha = 1.7, stable_beta = 0)), 1.7)
  expect_equal(gs_gev()$max_power(c(xi = 0.25)), 4)
  expect_equal(gs_gev()$max_power(c(xi = -0.25)), Inf)
  expect_equal(gs_gat()$max_power(c(nu = 2, d = 1.5, xi = 1)), 3)
})

test_that("closed form APARCH moments agree with numerical integration", {
  expect_equal(gs_aparch_moment(.gs_norm(), numeric(0), 0.5, 1.3),
               .aparch_moment_numeric(.gs_norm(), numeric(0), 0.5, 1.3), tolerance = 1e-6)
  pg <- c(nu = 2, d = 1.5, xi = 1.3)
  expect_equal(gs_aparch_moment(gs_gat(), pg, -0.4, 1.7),
               .aparch_moment_numeric(gs_gat(), pg, -0.4, 1.7), tolerance = 1e-6)
  # at beta = 0 S0 and S1 coincide, so the S1 closed form must match S0 integration
  ps <- c(stable_alpha = 1.6, stable_beta = 0)
  for (delta in c(0.7, 1, 1.4))
    expect_equal(gs_aparch_moment(gs_stable(), ps, 0.3, delta),
                 .aparch_moment_numeric(gs_stable(), ps, 0.3, delta), tolerance = 1e-5)
})

test_that("APARCH moment is infinite at or above max_power", {
  expect_equal(gs_aparch_moment(gs_stable(), c(stable_alpha = 1.5, stable_beta = 0.2), 0, 1.5), Inf)
  expect_equal(gs_aparch_moment(gs_gat(), c(nu = 1, d = 1, xi = 1), 0, 1.2), Inf)
})

test_that("stable APARCH moment with beta != 0 uses integration and is finite", {
  k <- gs_aparch_moment(gs_stable(), c(stable_alpha = 1.7, stable_beta = 0.3), 0.2, 1)
  expect_true(is.finite(k) && k > 0)
})

test_that("pgat keeps full accuracy next to zero", {
  eps <- 1e-6
  num <- (pgat(eps, nu = 2, d = 4) - pgat(-eps, nu = 2, d = 4)) / (2 * eps)
  expect_equal(num, dgat(0, nu = 2, d = 4), tolerance = 1e-6)
})

test_that("stable density is correct next to zeta, where libstable4u 1.0.5 is not", {
  skip_if_not_installed("stabledist")
  for (a in c(1.3, 1.75, 1.95)) for (b in c(-0.9, -0.3, 0.3, 0.9)) {
    zeta <- -b * tan(pi * a / 2)
    z <- zeta + c(-1e-3, -1e-4, -3e-5, -1e-5, -1e-6, 0, 1e-6, 1e-5, 3e-5, 1e-4, 1e-3)
    ours <- exp(gs_stable()$log_density(z, c(stable_alpha = a, stable_beta = b)))
    # stabledist is itself flat within ~1e-4 of zeta, so compare at 1e-4
    expect_equal(ours, stabledist::dstable(z, a, b, 1, 0, pm = 0), tolerance = 1e-4)
  }
})
