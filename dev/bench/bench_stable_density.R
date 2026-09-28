# Benchmark: stable density from libstable4u (used by the package) against
# stabledist, per density call and inside a full likelihood evaluation.
# Run from the package root: Rscript dev/bench/bench_stable_density.R

devtools::load_all(quiet = TRUE)

set.seed(1)
z <- libstable4u::stable_rnd(2000, c(1.75, 0.2, 1, 0), 0L)

per_call <- microbenchmark::microbenchmark(
  libstable4u = libstable4u::stable_pdf(z, c(1.75, 0.2, 1, 0), 0L),
  stabledist  = stabledist::dstable(z, 1.75, 0.2, 1, 0, pm = 0),
  times = 5)
print(per_call)

# One likelihood evaluation of a stable AR(1)-APARCH(1,1) on 2000 observations
spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
model <- gs_model(spec, mu = 0, ar = 0.05, omega = 0.05, alpha = 0.05, gamma = 0.3,
                  beta = 0.85, delta = 1.2, dist_par = c(stable_alpha = 1.75, stable_beta = 0.2))
x <- gs_sim(model, n = 2000, seed = 2)$y
llh <- microbenchmark::microbenchmark(
  gs_likelihood = .neg_loglik(model$par, x, spec),
  times = 10)
print(llh)

# Full fit, with and without standard errors, and with a warm start
print(system.time(fit <- gs_fit(x, spec)))
print(system.time(gs_fit(x, spec, hessian = FALSE)))
print(system.time(gs_fit(x, spec, start = coef(fit), hessian = FALSE)))
