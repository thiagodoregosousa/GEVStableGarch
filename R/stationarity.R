#' Persistence of the APARCH scale recursion
#'
#' Computes \eqn{\sum_j \beta_j + \sum_i \alpha_i \kappa_i}, with
#' \eqn{\kappa_i = E(|z| - \gamma_i z)^\delta} from [gs_aparch_moment()]. A
#' value below 1 means \eqn{E(\sigma_t^\delta)} is finite (Ding, Granger and
#' Engle, 1993; Mittnik, Paolella and Rachev, 2002). It is `Inf` when
#' \eqn{\delta} is not below the moment bound of the innovations.
#'
#' @param object A `gs_model` or `gs_fit` object.
#' @return The persistence, a single number.
#' @examples
#' spec <- gs_spec(garch = c(1, 1), dist = gs_stable())
#' m <- gs_model(spec, mu = 0, omega = 0.05, alpha = 0.1, beta = 0.8,
#'               dist_par = c(stable_alpha = 1.8, stable_beta = 0))
#' gs_stationarity(m)
#' @export
gs_stationarity <- function(object)
{
  model <- if (inherits(object, "gs_fit")) object$model else object
  if (!inherits(model, "gs_model")) stop("`object` must be a gs_model or gs_fit.", call. = FALSE)
  spec <- model$spec
  u <- .unpack(model$par, spec)
  kappa <- vapply(u$gamma, function(g) gs_aparch_moment(spec$dist, u$dist_par, g, u$delta), numeric(1))
  sum(u$beta) + sum(u$alpha * kappa)
}
