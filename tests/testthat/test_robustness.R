# Phase 1 robustness: graded penalty, convergence diagnostics, large-sample
# recovery and the parametric bootstrap.

test_that("the penalty slopes back toward the feasible region", {
  spec <- gs_spec(garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  x <- as.numeric(gs_sim(gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, gamma = 0,
                                  beta = 0.8, delta = 1,
                                  dist_par = c(stable_alpha = 1.8, stable_beta = 0)),
                         n = 200, seed = 1)$y)
  nms <- .par_names(spec)
  base <- c(mu = 0, omega = 0.1, alpha1 = 0.1, gamma1 = 0, beta1 = 0.8, delta = 1,
            stable_alpha = 1.8, stable_beta = 0)[nms]

  # delta >= stable_alpha (1.8) is infeasible; a larger breach must penalize more
  near <- replace(base, "delta", 1.9)   # breach 0.1
  far  <- replace(base, "delta", 3.0)   # breach 1.2
  v_near <- .neg_loglik(near, x, spec)
  v_far  <- .neg_loglik(far, x, spec)
  expect_gt(v_near, .PENALTY)            # infeasible, graded above the base penalty
  expect_gt(v_far, v_near)               # monotone in the size of the violation
})

test_that("a fit records gradient and conditioning diagnostics", {
  skip_if_not_installed("fGarch")
  data("dem2gbp", package = "fGarch", envir = environment())
  g <- gs_fit(dem2gbp[, 1], gs_spec(garch = c(1, 1), dist = .gs_norm()))
  d <- g$diagnostics
  expect_true(is.list(d))
  expect_true(is.finite(d$grad_rel) && d$grad_rel < 1e-2)   # near a stationary point
  expect_true(isTRUE(d$hess_pd))
  expect_true(is.finite(d$hess_rcond) && d$hess_rcond > 0)
  expect_false(isTRUE(d$cap_hit))
  expect_output(print(g), "scaled gradient")
})

test_that("GEV GARCH(1,1) estimates recover the truth at a large sample", {
  skip_on_cran()
  spec <- gs_spec(garch = c(1, 1), dist = gs_gev())
  # GEV innovations are standardized to scale 1, not unit variance, so E(z^2) ~ 3.2
  # at xi = 0.15: alpha * 3.25 + beta must stay below 1 for stationarity.
  truth <- c(mu = 0, omega = 0.05, alpha1 = 0.09, beta1 = 0.6, xi = 0.15)
  model <- gs_model(spec, mu = truth[["mu"]], omega = truth[["omega"]],
                    alpha = truth[["alpha1"]], beta = truth[["beta1"]],
                    dist_par = c(xi = truth[["xi"]]))
  expect_lt(gs_stationarity(model), 1)
  x <- gs_sim(model, n = 8000, seed = 42)$y
  g <- gs_fit(x, spec)
  se <- sqrt(diag(vcov(g)))
  expect_true(all(is.finite(se)))
  # every estimate within four standard errors of the true value
  expect_true(all(abs(coef(g)[names(truth)] - truth) < 4 * se[names(truth)]))
})

test_that("stable GARCH(1,1) estimates recover the truth at a large sample", {
  skip_on_cran()
  spec <- gs_spec(garch = c(1, 1), dist = gs_stable())
  truth <- c(mu = 0, omega = 0.04, alpha1 = 0.1, beta1 = 0.85,
             stable_alpha = 1.7, stable_beta = 0)
  model <- gs_model(spec, mu = truth[["mu"]], omega = truth[["omega"]],
                    alpha = truth[["alpha1"]], beta = truth[["beta1"]],
                    dist_par = truth[c("stable_alpha", "stable_beta")])
  x <- gs_sim(model, n = 6000, seed = 7)$y
  g <- gs_fit(x, spec)
  expect_equal(coef(g)[["stable_alpha"]], truth[["stable_alpha"]], tolerance = 0.1)
  expect_equal(coef(g)[["beta1"]], truth[["beta1"]], tolerance = 0.1)
})

test_that("gs_bootstrap returns standard errors and intervals of the right shape", {
  skip_on_cran()
  spec <- gs_spec(garch = c(1, 1), dist = gs_gev())
  model <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.07, beta = 0.6, dist_par = c(xi = 0.2))
  x <- gs_sim(model, n = 600, seed = 1)$y
  fit <- gs_fit(x, spec, hessian = FALSE)
  # an occasional bootstrap refit may fail and be dropped (warned); that is the
  # documented behaviour and is asserted through n_fail below.
  b <- suppressWarnings(gs_bootstrap(fit, n_boot = 20, seed = 1))
  expect_named(b$se, names(coef(fit)))
  expect_true(all(b$se > 0))
  expect_equal(dim(b$ci), c(length(coef(fit)), 2))
  expect_true(all(b$ci[, "lower"] <= b$ci[, "upper"]))
  expect_lt(b$n_fail, 20)
})
