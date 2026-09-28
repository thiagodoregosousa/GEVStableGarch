#' Stable innovations (S0 parametrization)
#'
#' Standardized stable distribution \eqn{S(\alpha, \beta, 1, 0; 0)} in
#' Nolan's S0 parametrization, computed with `libstable4u`. Parameters:
#' `stable_alpha` (tail index, in (1, 2)) and `stable_beta` (skewness, in
#' (-1, 1)).
#'
#' The variance is infinite, so \eqn{\sigma_t} is a conditional scale, not a
#' standard deviation. In GARCH mode the power is \eqn{\delta = 1} (absolute
#' value GARCH), since \eqn{E|z|^\delta} is finite only for
#' \eqn{\delta < \alpha}. The mean is \eqn{-\beta \tan(\pi\alpha/2)}.
#'
#' @param lower,upper,start Optional named vectors replacing the default
#'   bounds and start value.
#' @return A `gs_dist` object.
#' @references Nolan, J. P. (2020). Univariate Stable Distributions. Springer.
#' @examples
#' d <- gs_stable()
#' d$log_density(c(-1, 0, 1), c(stable_alpha = 1.7, stable_beta = 0.2))
#' @export
gs_stable <- function(lower = c(stable_alpha = 1.01, stable_beta = -0.99),
                      upper = c(stable_alpha = 1.99, stable_beta = 0.99),
                      start = c(stable_alpha = 1.8, stable_beta = 0))
{
  pars <- function(par) c(par[["stable_alpha"]], par[["stable_beta"]], 1, 0)
  gs_dist(
    name = "stable",
    log_density = function(z, par) log(.stable_pdf(z, par[["stable_alpha"]], par[["stable_beta"]])),
    random = function(n, par) libstable4u::stable_rnd(n, pars(par), 0L),
    cdf = function(q, par) libstable4u::stable_cdf(q, pars(par), 0L),
    quantile = function(p, par) libstable4u::stable_q(p, pars(par), 0L),
    par_names = c("stable_alpha", "stable_beta"),
    lower = lower, upper = upper, start = start,
    default_delta = 1,
    max_power = function(par) par[["stable_alpha"]],
    mean = function(par) -par[["stable_beta"]] * tan(pi * par[["stable_alpha"]] / 2),
    # closed form only for the symmetric case, where S0 and S1 coincide
    aparch_moment = function(gamma, delta, par)
      if (par[["stable_beta"]] == 0) .stable_sym_aparch_moment(par[["stable_alpha"]], gamma, delta) else NA,
    # McCulloch quantile estimates, kept inside the bounds
    start_fun = function(z) {
      est <- libstable4u::stable_fit_init(z, 0L)[1:2]
      est <- pmin(pmax(est, lower + 0.02), upper - 0.02)
      c(stable_alpha = est[[1]], stable_beta = est[[2]])
    },
    check = FALSE)
}

# libstable4u (1.0.5) returns about half the density for points at distance ~1e-5
# from zeta = -beta tan(pi alpha / 2), the special point of the S0 integral
# representation (wider when alpha is close to 1). The resulting jumps in the log
# likelihood break optimization and the Hessian. Inside a small window around zeta
# the density is replaced by the quadratic through three points where it is exact.
.stable_pdf <- function(z, alpha, beta)
{
  pars <- c(alpha, beta, 1, 0)
  f <- libstable4u::stable_pdf(z, pars, 0L)
  zeta <- -beta * tan(pi * alpha / 2)
  w <- if (alpha < 1.2) 0.02 else 1e-3
  near <- abs(z - zeta) < w
  if (any(near)) {
    y <- libstable4u::stable_pdf(zeta + c(-w, 0, w), pars, 0L)
    t <- (z[near] - zeta) / w
    f[near] <- y[2] + t * (y[3] - y[1]) / 2 + t^2 * (y[3] - 2 * y[2] + y[1]) / 2
  }
  f
}
