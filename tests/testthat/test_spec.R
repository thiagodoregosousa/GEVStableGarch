test_that("invalid specifications are rejected", {
  expect_error(gs_spec(garch = c(0, 1)), "p, number of alpha")
  expect_error(gs_spec(garch = c(1, -1)), "nonnegative integers")
  expect_error(gs_spec(arma = c(-1, 0)), "nonnegative integers")
  expect_error(gs_spec(aparch = "FALSE"), "TRUE or FALSE")
  expect_error(gs_spec(include_mean = "TRUE"), "TRUE or FALSE")
  expect_error(gs_spec(dist = "stable"), "gs_dist")
  expect_s3_class(gs_spec(), "gs_spec")
})

test_that("parameter names follow the specification", {
  spec <- gs_spec(arma = c(1, 2), garch = c(2, 1), aparch = TRUE, dist = gs_gat())
  expect_equal(.par_names(spec),
               c("mu", "ar1", "ma1", "ma2", "omega", "alpha1", "alpha2", "gamma1", "gamma2",
                 "beta1", "delta", "nu", "d", "xi"))
  spec0 <- gs_spec(garch = c(1, 0), include_mean = FALSE, dist = gs_gev())
  expect_equal(.par_names(spec0), c("omega", "alpha1", "xi"))
})

test_that("GARCH mode fixes gamma at 0 and delta at the distribution default", {
  spec <- gs_spec(garch = c(1, 1), dist = gs_stable())
  m <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.8,
                dist_par = c(stable_alpha = 1.8, stable_beta = 0))
  u <- .unpack(m$par, spec)
  expect_equal(u$gamma, 0)
  expect_equal(u$delta, 1)
})

test_that("gs_model checks lengths and admissible values", {
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
  dp <- c(stable_alpha = 1.8, stable_beta = 0)
  ok <- list(spec = spec, mu = 0, ar = 0.2, omega = 0.1, alpha = 0.1, beta = 0.8, dist_par = dp)
  expect_s3_class(do.call(gs_model, ok), "gs_model")
  expect_error(do.call(gs_model, modifyList(ok, list(ar = c(0.1, 0.1)))), "`ar`")
  expect_error(do.call(gs_model, modifyList(ok, list(ma = 0.1))), "`ma` must not")
  expect_error(do.call(gs_model, modifyList(ok, list(omega = -1))), "omega must be positive")
  expect_error(do.call(gs_model, modifyList(ok, list(ar = 1.2))), "AR polynomial")
  expect_error(do.call(gs_model, modifyList(ok, list(delta = 1))), "`delta` must not")
  expect_error(do.call(gs_model, modifyList(ok, list(dist_par = c(stable_alpha = 2.5, stable_beta = 0)))),
               "outside their bounds")
  expect_error(do.call(gs_model, modifyList(ok, list(dist_par = c(alpha = 1.8)))), "dist_par")
})

test_that("APARCH power must be below the moment bound of the innovations", {
  spec <- gs_spec(garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  expect_error(gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, gamma = 0, beta = 0.8, delta = 1.9,
                        dist_par = c(stable_alpha = 1.7, stable_beta = 0)), "max_power")
})

test_that("gs_model accepts a named parameter vector in any order", {
  spec <- gs_spec(garch = c(1, 1), dist = gs_gev())
  m <- gs_model(spec, par = c(xi = 0.1, beta1 = 0.8, alpha1 = 0.1, omega = 0.1, mu = 0))
  expect_equal(names(m$par), .par_names(spec))
})
