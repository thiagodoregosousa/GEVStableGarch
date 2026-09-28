#' Specify an ARMA-GARCH/APARCH model
#'
#' The model, in the intercept form of Wuertz et al. (2006) and fGarch, is
#' \deqn{x_t = \mu + \sum_{i=1}^m a_i x_{t-i} + \sum_{j=1}^n b_j e_{t-j} + e_t,
#'       \qquad e_t = \sigma_t z_t,}
#' \deqn{\sigma_t^\delta = \omega + \sum_{i=1}^p \alpha_i (|e_{t-i}| - \gamma_i e_{t-i})^\delta
#'       + \sum_{j=1}^q \beta_j \sigma_{t-j}^\delta,}
#' with \eqn{z_t} i.i.d. from `dist`. In GARCH mode (`aparch = FALSE`) the
#' asymmetry \eqn{\gamma_i} is 0 and \eqn{\delta} is fixed at
#' `dist$default_delta` (1 for stable, 2 otherwise). `mu` is an intercept:
#' the unconditional location is \eqn{\mu / (1 - \sum a_i)}. Since the
#' innovations are location scale (not standardized), \eqn{\sigma_t} is a
#' conditional scale.
#'
#' @param arma Integer vector `c(m, n)`, ARMA orders.
#' @param garch Integer vector `c(p, q)`, number of \eqn{\alpha} (p >= 1) and
#'   \eqn{\beta} (q >= 0) terms.
#' @param aparch Logical, estimate the asymmetry \eqn{\gamma} and the power
#'   \eqn{\delta}.
#' @param include_mean Logical, include the intercept \eqn{\mu}.
#' @param dist Innovation distribution, a `gs_dist` object.
#' @return An object of class `gs_spec`.
#' @references Wuertz, D., Chalabi, Y. and Luksan, L. (2006). Parameter
#'   estimation of ARMA models with GARCH/APARCH errors: an R and SPlus
#'   software implementation. Journal of Statistical Software.
#' @examples
#' gs_spec(arma = c(1, 0), garch = c(1, 1), aparch = TRUE, dist = gs_stable())
#' @export
gs_spec <- function(arma = c(0, 0), garch = c(1, 1), aparch = FALSE,
                    include_mean = TRUE, dist = gs_stable())
{
  .check_order(arma, "arma")
  .check_order(garch, "garch")
  if (garch[1] < 1) stop("`garch[1]` (p, number of alpha terms) must be >= 1.", call. = FALSE)
  .check_flag(aparch, "aparch")
  .check_flag(include_mean, "include_mean")
  .check_dist_structure(dist)

  spec <- list(order = c(m = as.integer(arma[1]), n = as.integer(arma[2]),
                         p = as.integer(garch[1]), q = as.integer(garch[2])),
               aparch = aparch, include_mean = include_mean, dist = dist)
  class(spec) <- "gs_spec"
  spec
}


#' Set the parameters of an ARMA-GARCH/APARCH model
#'
#' Attaches parameter values to a [gs_spec()], checking lengths, admissible
#' ranges and moment conditions. Parameters can be given one by one or as a
#' named vector `par` (for example `coef(fit)`).
#'
#' @param spec A `gs_spec` object.
#' @param mu Intercept (only when `spec$include_mean`).
#' @param ar,ma AR and MA coefficients, lengths m and n.
#' @param omega Positive constant of the scale recursion.
#' @param alpha,beta Nonnegative coefficients, lengths p and q.
#' @param gamma Asymmetry coefficients in (-1, 1), length p (APARCH only).
#' @param delta Positive power (APARCH only).
#' @param dist_par Distribution parameters, named by `spec$dist$par_names`.
#' @param par Alternatively, a named vector with all free parameters.
#' @return An object of class `gs_model`.
#' @examples
#' spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
#' gs_model(spec, mu = 0, ar = 0.1, omega = 0.05, alpha = 0.1, beta = 0.85,
#'          dist_par = c(stable_alpha = 1.8, stable_beta = 0))
#' @export
gs_model <- function(spec, mu = NULL, ar = NULL, ma = NULL, omega = NULL,
                     alpha = NULL, gamma = NULL, beta = NULL, delta = NULL,
                     dist_par = NULL, par = NULL)
{
  if (!inherits(spec, "gs_spec")) stop("`spec` must be a gs_spec object.", call. = FALSE)
  if (is.null(par)) {
    o <- spec$order
    .check_len(mu, if (spec$include_mean) 1 else 0, "mu", "include_mean")
    .check_len(ar, o[["m"]], "ar", "m")
    .check_len(ma, o[["n"]], "ma", "n")
    .check_len(omega, 1, "omega", "the model")
    .check_len(alpha, o[["p"]], "alpha", "p")
    .check_len(gamma, if (spec$aparch) o[["p"]] else 0, "gamma", "aparch")
    .check_len(beta, o[["q"]], "beta", "q")
    .check_len(delta, if (spec$aparch) 1 else 0, "delta", "aparch")
    dist_par <- unlist(dist_par)
    if (!setequal(names(dist_par), spec$dist$par_names))
      stop("`dist_par` must be named: ", paste(spec$dist$par_names, collapse = ", "), ".", call. = FALSE)
    par <- c(mu, ar, ma, omega, alpha, gamma, beta, delta, dist_par[spec$dist$par_names])
    names(par) <- .par_names(spec)
  }
  model <- list(spec = spec, par = .validate_par(par, spec))
  class(model) <- "gs_model"
  model
}


#' @export
print.gs_spec <- function(x, ...)
{
  cat(.model_label(x), "with", x$dist$name, "innovations\n")
  cat("Parameters:", paste(.par_names(x), collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.gs_model <- function(x, ...)
{
  cat(.model_label(x$spec), "with", x$spec$dist$name, "innovations\n")
  print(x$par)
  invisible(x)
}


# ------------------------------------------------------------------------------
# Internal helpers

.model_label <- function(spec)
{
  o <- spec$order
  mean_part <- if (o[["m"]] + o[["n"]] > 0) sprintf("ARMA(%d,%d)-", o[["m"]], o[["n"]]) else ""
  sprintf("%s%s(%d,%d)", mean_part, if (spec$aparch) "APARCH" else "GARCH", o[["p"]], o[["q"]])
}

.check_order <- function(x, what)
{
  if (!is.numeric(x) || length(x) != 2 || any(is.na(x)) || any(x < 0) || any(x %% 1 != 0))
    stop(sprintf("`%s` must be two nonnegative integers.", what), call. = FALSE)
}

.check_flag <- function(x, what)
{
  if (!is.logical(x) || length(x) != 1 || is.na(x))
    stop(sprintf("`%s` must be TRUE or FALSE.", what), call. = FALSE)
}

.check_len <- function(x, k, what, because)
{
  if (k == 0 && !is.null(x))
    stop(sprintf("`%s` must not be supplied (see `%s`).", what, because), call. = FALSE)
  if (k > 0 && (is.null(x) || !is.numeric(x) || length(x) != k || any(!is.finite(x))))
    stop(sprintf("`%s` must be numeric of length %d.", what, k), call. = FALSE)
}

# Every free parameter named and admissible; returns the vector in canonical order
.validate_par <- function(par, spec)
{
  nms <- .par_names(spec)
  if (is.null(names(par)) || !setequal(names(par), nms))
    stop("Parameters must be named: ", paste(nms, collapse = ", "), ".", call. = FALSE)
  par <- par[nms]
  why <- .inadmissible(.unpack(par, spec), spec)
  if (!is.null(why)) stop("Inadmissible parameters: ", why, ".", call. = FALSE)
  par
}
