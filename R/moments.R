#' APARCH moment of the innovation distribution
#'
#' Computes \eqn{\kappa = E(|z| - \gamma z)^\delta} for standardized
#' innovations \eqn{z}, the quantity that enters the APARCH stationarity
#' condition and the multi step scale forecast. A closed form is used when the
#' distribution provides one, otherwise numerical integration.
#'
#' @param dist A `gs_dist` object.
#' @param par Distribution parameters, named by `dist$par_names`.
#' @param gamma Asymmetry parameter in (-1, 1).
#' @param delta Power, positive.
#' @return \eqn{\kappa}, or `Inf` when \eqn{\delta} is not below
#'   `dist$max_power(par)`.
#' @references Ding, Z., Granger, C. W. J. and Engle, R. F. (1993). A long
#'   memory property of stock market returns and a new model. Journal of
#'   Empirical Finance, 1, 83-106.
#'
#'   Mittnik, S., Paolella, M. S. and Rachev, S. T. (2002). Stationarity of
#'   stable power-GARCH processes. Journal of Econometrics, 106, 97-107.
#' @examples
#' gs_aparch_moment(gs_stable(), c(stable_alpha = 1.7, stable_beta = 0),
#'                  gamma = 0.3, delta = 1)
#' @export
gs_aparch_moment <- function(dist, par, gamma, delta)
{
  par <- .named(par, dist$par_names)
  if (abs(gamma) >= 1 || delta <= 0) stop("Need |gamma| < 1 and delta > 0.", call. = FALSE)
  if (delta >= dist$max_power(par)) return(Inf)
  if (!is.null(dist$aparch_moment)) {
    closed <- dist$aparch_moment(gamma, delta, par)
    if (!is.na(closed)) return(closed)
  }
  .aparch_moment_numeric(dist, par, gamma, delta)
}


# Numerical E(|z| - gamma z)^delta, split at zero where the integrand has a kink
.aparch_moment_numeric <- function(dist, par, gamma, delta)
{
  g <- function(x) (abs(x) - gamma * x)^delta
  .integrate_dist(g, dist, par)
}


# ------------------------------------------------------------------------------
# Closed forms

# Standard normal, Ding, Granger and Engle (1993), appendix B
.norm_aparch_moment <- function(gamma, delta)
{
  1 / sqrt(2 * pi) * ((1 + gamma)^delta + (1 - gamma)^delta) *
    2^((delta - 1) / 2) * gamma((delta + 1) / 2)
}

# Symmetric stable with unit scale (S0 and S1 coincide when beta = 0), Diongue (2008)
.stable_sym_aparch_moment <- function(alpha, gamma, delta)
{
  k <- if (delta == 1) pi / 2 else gamma(1 - delta) * cos(pi * delta / 2)
  gamma(1 - delta / alpha) / k * 0.5 * ((1 + gamma)^delta + (1 - gamma)^delta)
}

# GAt distribution, Mittnik and Paolella (2000)
.gat_aparch_moment <- function(nu, d, xi, gamma, delta)
{
  a <- (1 + gamma)^delta * xi^(-delta - 1) + (1 - gamma)^delta * xi^(delta + 1)
  b <- nu^(delta / d) * beta((delta + 1) / d, nu - delta / d)
  c <- (1 / xi + xi) * beta(1 / d, nu)
  a * b / c
}
