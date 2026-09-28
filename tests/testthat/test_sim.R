test_that("simulation is reproducible with a seed", {
  spec <- gs_spec(arma = c(1, 1), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  m <- gs_model(spec, mu = 0.01, ar = 0.2, ma = -0.1, omega = 0.05, alpha = 0.08, gamma = 0.3,
                beta = 0.85, delta = 1.2, dist_par = c(stable_alpha = 1.8, stable_beta = 0.2))
  a <- gs_sim(m, n = 300, seed = 42); b <- gs_sim(m, n = 300, seed = 42)
  expect_identical(a, b)
  expect_equal(dim(a), c(300, 3))
  expect_true(all(is.finite(as.matrix(a))))
})

test_that("simulation and filter implement the same equations", {
  spec <- gs_spec(arma = c(2, 1), garch = c(2, 1), aparch = TRUE, dist = gs_gat())
  m <- gs_model(spec, mu = 0.1, ar = c(0.3, -0.1), ma = 0.2, omega = 0.05, alpha = c(0.05, 0.03),
                gamma = c(0.2, -0.1), beta = 0.85, delta = 1.5, dist_par = c(nu = 2, d = 3, xi = 1.1))
  s <- gs_sim(m, n = 2000, seed = 1)
  f <- .filter_model(s$y, .unpack(m$par, spec))
  # different initial conditions, then the filter locks onto the simulated scale
  tail_idx <- 1901:2000
  expect_equal(f$sigma[tail_idx], s$sigma[tail_idx], tolerance = 1e-6)
  expect_equal(f$e[tail_idx] / f$sigma[tail_idx], s$z[tail_idx], tolerance = 1e-6)
})

test_that("normal GARCH(1,1) simulation reaches the theoretical variance", {
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  m <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.8)
  s <- gs_sim(m, n = 1e5, seed = 3)
  expect_equal(var(s$y), 0.1 / (1 - 0.1 - 0.8), tolerance = 0.05)
})

test_that("non stationary models warn", {
  spec <- gs_spec(garch = c(1, 1), dist = .gs_norm())
  m <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.3, beta = 0.75)
  expect_warning(gs_sim(m, n = 100, seed = 1), "Persistence")
})

# Simulate and estimate: estimates within 4 standard errors of the truth
round_trip <- function(dist, dist_par, alpha = 0.1, n = 3000, seed = 11)
{
  spec <- gs_spec(garch = c(1, 1), dist = dist)
  m <- gs_model(spec, mu = 0.05, omega = 0.05, alpha = alpha, beta = 0.8, dist_par = dist_par)
  x <- gs_sim(m, n = n, seed = seed)$y
  fit <- gs_fit(x, spec)
  z <- (coef(fit) - m$par) / sqrt(diag(vcov(fit)))
  list(fit = fit, z = z)
}

test_that("round trip: stable GARCH(1,1)", {
  skip_on_cran()
  r <- round_trip(gs_stable(), c(stable_alpha = 1.7, stable_beta = 0.3))
  expect_true(all(abs(r$z) < 4), info = paste(round(r$z, 2), collapse = " "))
})

test_that("round trip: GEV GARCH(1,1)", {
  skip_on_cran()
  # GEV innovations are not centered (E z^2 near 2.8), so a smaller alpha keeps it stationary
  r <- round_trip(gs_gev(), c(xi = 0.1), alpha = 0.05)
  expect_true(all(abs(r$z) < 4), info = paste(round(r$z, 2), collapse = " "))
})

test_that("round trip: GAt GARCH(1,1)", {
  skip_on_cran()
  r <- round_trip(gs_gat(), c(nu = 2, d = 2.5, xi = 1.2))
  expect_true(all(abs(r$z) < 4), info = paste(round(r$z, 2), collapse = " "))
})

test_that("round trip: stable AR(1)-APARCH(1,1), the model of the VaR paper", {
  skip_on_cran()
  spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
  m <- gs_model(spec, mu = 0.0005, ar = 0.05, omega = 0.0002, alpha = 0.05, gamma = 0.3, beta = 0.85,
                delta = 1.2, dist_par = c(stable_alpha = 1.75, stable_beta = 0))
  x <- gs_sim(m, n = 2000, seed = 5)$y
  fit <- gs_fit(x, spec)
  z <- (coef(fit) - m$par) / sqrt(diag(vcov(fit)))
  expect_true(all(abs(z) < 4), info = paste(round(z, 2), collapse = " "))
})
