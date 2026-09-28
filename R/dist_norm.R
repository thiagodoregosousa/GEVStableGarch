# Standard normal innovations written with gs_dist(). Not exported: it is the
# worked example of a user defined distribution and the oracle for tests
# against fGarch.
.gs_norm <- function()
{
  gs_dist(
    name = "norm",
    log_density = function(z, par) stats::dnorm(z, log = TRUE),
    random = function(n, par) stats::rnorm(n),
    cdf = function(q, par) stats::pnorm(q),
    quantile = function(p, par) stats::qnorm(p),
    par_names = character(0),
    lower = numeric(0), upper = numeric(0), start = numeric(0),
    default_delta = 2,
    max_power = function(par) Inf,
    mean = function(par) 0,
    aparch_moment = function(gamma, delta, par) .norm_aparch_moment(gamma, delta),
    check = FALSE)
}
