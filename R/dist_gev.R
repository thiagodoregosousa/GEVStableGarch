#' GEV innovations
#'
#' Generalized extreme value distribution with location 0, scale 1 and shape
#' `xi`, \eqn{F(z) = \exp\{-(1 + \xi z)^{-1/\xi}\}} on \eqn{1 + \xi z > 0}
#' (Gumbel when \eqn{\xi = 0}). The shape is restricted to (-0.5, 0.5):
#' above -0.5 for regular likelihood properties, below 0.5 for a finite
#' variance, which the GARCH mode power \eqn{\delta = 2} requires.
#'
#' The innovations are not centered: \eqn{E[z] = (\Gamma(1-\xi) - 1)/\xi}, so
#' the model parameter `mu` is a location, not the conditional mean.
#'
#' @inheritParams gs_stable
#' @return A `gs_dist` object.
#' @examples
#' d <- gs_gev()
#' d$mean(c(xi = 0.1))
#' @export
gs_gev <- function(lower = c(xi = -0.49), upper = c(xi = 0.49), start = c(xi = 0.01))
{
  gumbel <- function(xi) abs(xi) < 1e-7
  gs_dist(
    name = "gev",
    log_density = function(z, par) {
      xi <- par[["xi"]]
      if (gumbel(xi)) return(-z - exp(-z))
      y <- 1 + xi * z
      out <- rep(-Inf, length(z))
      inside <- y > 0
      out[inside] <- -(1 + 1 / xi) * log(y[inside]) - y[inside]^(-1 / xi)
      out
    },
    random = function(n, par) .gev_quantile(stats::runif(n), par[["xi"]]),
    cdf = function(q, par) {
      xi <- par[["xi"]]
      if (gumbel(xi)) return(exp(-exp(-q)))
      y <- pmax(1 + xi * q, 0)
      exp(-y^(-1 / xi))
    },
    quantile = function(p, par) .gev_quantile(p, par[["xi"]]),
    par_names = "xi",
    lower = lower, upper = upper, start = start,
    default_delta = 2,
    max_power = function(par) if (par[["xi"]] > 0) 1 / par[["xi"]] else Inf,
    mean = function(par) {
      xi <- par[["xi"]]
      if (gumbel(xi)) return(-digamma(1))
      if (xi >= 1) return(NA_real_)
      (gamma(1 - xi) - 1) / xi
    },
    check = FALSE)
}

.gev_quantile <- function(p, xi)
{
  if (abs(xi) < 1e-7) return(-log(-log(p)))
  ((-log(p))^(-xi) - 1) / xi
}
