# The pluggable interface, exercised with a genuinely external distribution: the
# bimodal GEV from the CRAN package bgev (Otiniano et al. 2025). Standardized to
# location 0, scale 1, with two shape parameters (renamed to avoid clashing with
# the APARCH power `delta`). Moment bound from Remark 1: E|z|^k < Inf iff
# k < (delta + 1) / xi for xi > 0.

bgev_dist <- function() {
  bgev_mean <- function(xi, delta)
    tryCatch(stats::integrate(function(p) bgev::qbgev(p, 0, 1, xi, delta), 0, 1,
                              subdivisions = 400L, rel.tol = 1e-6)$value,
             error = function(e) NA_real_)
  gs_dist(
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
    default_delta = 2,
    max_power = function(par) (par[["bgev_delta"]] + 1) / par[["bgev_xi"]],
    mean = function(par) bgev_mean(par[["bgev_xi"]], par[["bgev_delta"]]))
}

test_that("a CRAN bgev distribution passes the gs_dist self check", {
  skip_if_not_installed("bgev")
  expect_silent(gs_check_dist(bgev_dist()))
})

test_that("bgev max_power matches the (delta+1)/xi moment bound", {
  skip_if_not_installed("bgev")
  d <- bgev_dist()
  expect_equal(d$max_power(c(bgev_xi = 0.2, bgev_delta = 2)), 3 / 0.2)
  # E|z|^k is finite just below the bound and infinite above it
  p <- c(bgev_xi = 0.3, bgev_delta = 1)      # bound = 2/0.3 = 6.67
  expect_true(is.finite(gs_aparch_moment(d, p, gamma = 0, delta = 6)))
  expect_equal(gs_aparch_moment(d, p, gamma = 0, delta = 7), Inf)
})

test_that("GARCH(1,1) with bgev innovations fits from the default start and recovers", {
  skip_on_cran()
  skip_if_not_installed("bgev")
  spec <- gs_spec(garch = c(1, 1), dist = bgev_dist())
  truth <- c(mu = 0, omega = 0.1, alpha1 = 0.1, beta1 = 0.6, bgev_xi = 0.2, bgev_delta = 2)
  model <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.6,
                    dist_par = c(bgev_xi = 0.2, bgev_delta = 2))
  x <- gs_sim(model, n = 3000, seed = 7)$y
  fit <- suppressWarnings(gs_fit(x, spec, hessian = FALSE))   # default start: rescued automatically
  expect_equal(fit$convergence, 0)
  expect_equal(coef(fit)[c("alpha1", "beta1")], truth[c("alpha1", "beta1")], tolerance = 0.15)
  expect_lt(abs(coef(fit)[["bgev_xi"]] - 0.2), 0.1)
  expect_lt(abs(coef(fit)[["bgev_delta"]] - 2), 0.7)
})
