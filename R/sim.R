#' Simulate an ARMA-GARCH/APARCH model
#'
#' Simulates the model equations of [gs_spec()] with innovations drawn from
#' the model's distribution. The recursion starts at the unconditional
#' location and scale level and the first `burnin` values are discarded.
#'
#' @param model A `gs_model` object.
#' @param n Length of the returned series.
#' @param burnin Number of initial values discarded.
#' @param seed Optional seed passed to [set.seed()].
#' @return A data frame with columns `y` (series), `sigma` (conditional scale)
#'   and `z` (innovations).
#' @examples
#' spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
#' model <- gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.8,
#'                   dist_par = c(stable_alpha = 1.8, stable_beta = 0))
#' x <- gs_sim(model, n = 500, seed = 1)
#' plot(x$y, type = "l")
#' @export
gs_sim <- function(model, n = 1000, burnin = 1000, seed = NULL)
{
  if (!inherits(model, "gs_model")) stop("`model` must be a gs_model object.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  spec <- model$spec
  u <- .unpack(model$par, spec)
  persistence <- gs_stationarity(model)
  if (persistence >= 1)
    warning(sprintf("Persistence %.4f >= 1: the simulated scale has no stationary level.", persistence),
            call. = FALSE)

  z <- matrix(spec$dist$random(n + burnin, u$dist_par), ncol = 1)
  paths <- .sim_paths(u, z, .initial_state(u, persistence))
  keep <- burnin + seq_len(n)
  data.frame(y = paths$x[keep, 1], sigma = paths$sigma[keep, 1], z = z[keep, 1])
}


# ------------------------------------------------------------------------------
# Recursion engine shared by simulation and forecasting

# Number of past values the recursion needs
.max_lag <- function(u) max(length(u$ar), length(u$ma), length(u$alpha), length(u$beta), 1)

# Start at the unconditional location, zero residuals and the unconditional h level
.initial_state <- function(u, persistence)
{
  L <- .max_lag(u)
  level <- if (persistence < 1) u$omega / (1 - persistence) else u$omega
  list(x = rep(u$mu / (1 - sum(u$ar)), L), e = rep(0, L), h = rep(level, L))
}

# Run the model forward. `z` is a (horizon x paths) matrix of innovations and `state`
# holds the last .max_lag(u) values of x, e and h (oldest first), shared by all paths.
.sim_paths <- function(u, z, state)
{
  H <- nrow(z); P <- ncol(z); L <- .max_lag(u)
  x <- e <- h <- matrix(0, L + H, P)
  x[1:L, ] <- state$x; e[1:L, ] <- state$e; h[1:L, ] <- state$h
  for (t in L + seq_len(H)) {
    ht <- u$omega
    for (i in seq_along(u$alpha)) ht <- ht + u$alpha[i] * (abs(e[t - i, ]) - u$gamma[i] * e[t - i, ])^u$delta
    for (j in seq_along(u$beta)) ht <- ht + u$beta[j] * h[t - j, ]
    h[t, ] <- ht
    e[t, ] <- ht^(1 / u$delta) * z[t - L, ]
    xt <- u$mu + e[t, ]
    for (i in seq_along(u$ar)) xt <- xt + u$ar[i] * x[t - i, ]
    for (j in seq_along(u$ma)) xt <- xt + u$ma[j] * e[t - j, ]
    x[t, ] <- xt
  }
  rows <- L + seq_len(H)
  list(x = x[rows, , drop = FALSE], e = e[rows, , drop = FALSE],
       h = h[rows, , drop = FALSE], sigma = h[rows, , drop = FALSE]^(1 / u$delta))
}
