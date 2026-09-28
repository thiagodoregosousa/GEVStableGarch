#' Forecast an ARMA-GARCH/APARCH model
#'
#' Forecasts from a fitted model (`predict.gs_fit`) or from a model with
#' given parameters and a data history (`predict.gs_model`, useful in rolling
#' windows with fixed parameters).
#'
#' One step ahead everything is exact: the conditional location
#' \eqn{m_{T+1} = \mu + \sum a_i x_{T+1-i} + \sum b_j e_{T+1-j}}, the
#' conditional scale \eqn{\sigma_{T+1}}, and the quantiles
#' \eqn{m_{T+1} + \sigma_{T+1} F^{-1}(p)} (the value at risk at level `p` is
#' minus the `p` quantile). For horizons \eqn{k \ge 2}:
#' * `sigma` is \eqn{(E\,\sigma^\delta_{T+k})^{1/\delta}} from the APARCH
#'   recursion with \eqn{\kappa = E(|z| - \gamma z)^\delta} (as in fGarch);
#' * quantiles come from `n_sim` simulated paths, since the predictive
#'   distribution is a scale mixture and has no closed form;
#' * `mean` is exact when \eqn{E[z] = 0} or \eqn{\delta = 1}, otherwise it is
#'   the average of the simulated paths (column `mean_method`).
#'
#' With stable innovations there is no forecast variance; `sigma` is a scale.
#' `mean` is `NA` when the innovations have no mean.
#'
#' @param object A `gs_fit` or `gs_model` object.
#' @param data For `gs_model`, the observed history (numeric vector).
#' @param n_ahead Forecast horizon.
#' @param level Probabilities of the reported quantiles.
#' @param n_sim Number of simulated paths used when `n_ahead > 1`.
#' @param seed Optional seed for the simulated paths.
#' @param ... Unused.
#' @return A data frame with columns `horizon`, `location` (one step only),
#'   `mean`, `mean_method`, `sigma` and one column per `level` named `q_<level>`.
#' @examples
#' spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
#' model <- gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.8,
#'                   dist_par = c(stable_alpha = 1.8, stable_beta = 0))
#' x <- gs_sim(model, n = 500, seed = 1)$y
#' predict(model, data = x, n_ahead = 5, n_sim = 2000, seed = 1)
#' @export
predict.gs_fit <- function(object, n_ahead = 1, level = c(0.01, 0.05), n_sim = 10000,
                           seed = NULL, ...)
  .forecast(object$model, object$data, n_ahead, level, n_sim, seed)

#' @rdname predict.gs_fit
#' @export
predict.gs_model <- function(object, data, n_ahead = 1, level = c(0.01, 0.05), n_sim = 10000,
                             seed = NULL, ...)
  .forecast(object, as.numeric(data), n_ahead, level, n_sim, seed)


.forecast <- function(model, x, n_ahead, level, n_sim, seed)
{
  if (n_ahead < 1 || n_ahead %% 1 != 0) stop("`n_ahead` must be a positive integer.", call. = FALSE)
  if (any(level <= 0 | level >= 1)) stop("`level` must be in (0, 1).", call. = FALSE)
  spec <- model$spec; dist <- spec$dist
  u <- .unpack(model$par, spec)
  L <- .max_lag(u); N <- length(x)
  if (N < L + 1) stop("`data` is shorter than the model lags.", call. = FALSE)
  if (n_ahead > 1 && gs_stationarity(model) >= 1)
    warning("Persistence >= 1: multi step scale forecasts do not revert to a level.", call. = FALSE)

  # State at time T from the filtered history
  f <- .filter_model(x, u)
  last <- N - L + seq_len(L)
  state <- list(x = x[last], e = f$e[last], h = f$h[last])

  # One step ahead: zero innovations give the conditional location and the exact scale
  step1 <- .sim_paths(u, matrix(0, 1, 1), state)
  loc1 <- step1$x[1, 1]; sig1 <- step1$sigma[1, 1]

  # Expected h by the APARCH recursion; future lagged shocks replaced by kappa * E h
  kappa <- vapply(u$gamma, function(g) gs_aparch_moment(dist, u$dist_par, g, u$delta), numeric(1))
  arch_past <- sapply(seq_along(u$alpha), function(i) (abs(state$e) - u$gamma[i] * state$e)^u$delta)
  arch_past <- matrix(arch_past, nrow = L)
  Eh <- c(state$h, numeric(n_ahead))
  for (k in seq_len(n_ahead)) {
    t <- L + k; v <- u$omega
    for (i in seq_along(u$alpha))
      v <- v + u$alpha[i] * (if (t - i <= L) arch_past[t - i, i] else kappa[i] * Eh[t - i])
    for (j in seq_along(u$beta)) v <- v + u$beta[j] * Eh[t - j]
    Eh[t] <- v
  }
  sigma <- Eh[L + seq_len(n_ahead)]^(1 / u$delta)

  # Mean: exact via the ARMA recursion when E[sigma] E[z] is known exactly
  Ez <- dist$mean(u$dist_par)
  exact_mean <- is.na(Ez) || Ez == 0 || u$delta == 1
  Ex <- c(state$x, numeric(n_ahead)); Ee <- c(state$e, sigma * if (is.na(Ez)) 0 else Ez)
  Ee[L + 1] <- sig1 * if (is.na(Ez)) 0 else Ez
  for (k in seq_len(n_ahead)) {
    t <- L + k; v <- u$mu + Ee[t]
    for (i in seq_along(u$ar)) v <- v + u$ar[i] * Ex[t - i]
    for (j in seq_along(u$ma)) v <- v + u$ma[j] * Ee[t - j]
    Ex[t] <- v
  }
  mean <- if (is.na(Ez)) rep(NA_real_, n_ahead) else Ex[L + seq_len(n_ahead)]
  mean_method <- rep("exact", n_ahead)

  # Quantiles: exact at one step, simulated paths beyond
  q <- matrix(NA_real_, n_ahead, length(level))
  q[1, ] <- loc1 + sig1 * dist$quantile(level, u$dist_par)
  if (n_ahead > 1) {
    if (!is.null(seed)) set.seed(seed)
    z <- matrix(dist$random(n_ahead * n_sim, u$dist_par), nrow = n_ahead)
    paths <- .sim_paths(u, z, state)
    q[-1, ] <- t(apply(paths$x[-1, , drop = FALSE], 1, stats::quantile, probs = level, names = FALSE))
    if (!exact_mean) {
      mean[-1] <- rowMeans(paths$x[-1, , drop = FALSE])
      mean_method[-1] <- "simulated"
    }
  }

  out <- data.frame(horizon = seq_len(n_ahead), location = c(loc1, rep(NA_real_, n_ahead - 1)),
                    mean = mean, mean_method = mean_method, sigma = sigma)
  q <- as.data.frame(q); names(q) <- paste0("q_", level)
  cbind(out, q)
}
