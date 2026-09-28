#' Fit an ARMA-GARCH/APARCH model by maximum likelihood
#'
#' Estimates the model of [gs_spec()] by conditional maximum likelihood. The
#' data is rescaled internally to unit standard deviation (as in fGarch) and
#' the estimates are transformed back. Parameters move in plain intervals;
#' stationarity is not imposed but reported after fitting (see
#' [gs_stationarity()]), with a warning when the persistence is 1 or more.
#'
#' Standard errors come from a numerical Hessian (central differences
#' with steps adapted to the curvature). When a parameter sits on a bound, or a
#' differentiation step leaves the parameter space, the affected standard
#' errors are `NA` and a warning explains why.
#'
#' @param data Numeric vector, the time series (for example log returns).
#' @param spec A `gs_spec` object.
#' @param algorithm `"sqp"` ([Rsolnp::solnp()]) or `"nlminb"`.
#' @param start Optional named vector of start values on the original scale,
#'   for example `coef(previous_fit)` in rolling windows. Missing names are
#'   filled with the default start values.
#' @param hessian Logical, compute standard errors. Set `FALSE` to save time
#'   when only point estimates are needed.
#' @param control List passed to the optimizer.
#' @return An object of class `gs_fit` with methods [coef()], [vcov()],
#'   [logLik()], [residuals()], [sigma()], [fitted()], [predict()] and
#'   `print()`.
#' @examples
#' \donttest{
#' if (requireNamespace("fGarch", quietly = TRUE)) {
#'   data("dem2gbp", package = "fGarch")
#'   spec <- gs_spec(arma = c(1, 0), garch = c(1, 1), dist = gs_stable())
#'   fit <- gs_fit(dem2gbp[, 1], spec)
#'   fit
#' }
#' }
#' @export
gs_fit <- function(data, spec, algorithm = c("sqp", "nlminb"), start = NULL,
                   hessian = TRUE, control = list())
{
  algorithm <- match.arg(algorithm)
  if (!inherits(spec, "gs_spec")) stop("`spec` must be a gs_spec object.", call. = FALSE)
  x <- as.numeric(data)
  if (length(x) < 50 || any(!is.finite(x)))
    stop("`data` must be a finite numeric vector with at least 50 observations.", call. = FALSE)
  N <- length(x); s <- stats::sd(x); xs <- x / s
  bounds <- .opt_bounds(spec, xs)

  # Start values: defaults, overwritten by the user's (original scale) values
  par0 <- .start_values(spec, xs, bounds)
  if (!is.null(start)) {
    unknown <- setdiff(names(start), names(par0))
    if (length(unknown)) stop("Unknown names in `start`: ", paste(unknown, collapse = ", "), call. = FALSE)
    full <- .rescale_par(par0, spec, s, to_unit = FALSE)
    full[names(start)] <- start
    par0 <- .clamp(.rescale_par(full, spec, s, to_unit = TRUE), bounds)
  }
  diag <- new.env()
  if (.neg_loglik(par0, xs, spec, diag = diag) >= .PENALTY)
    stop("Start values are inadmissible: ", diag$reason, call. = FALSE)

  # Optimize on the unit scale
  obj <- function(p) .neg_loglik(p, xs, spec)
  opt <- if (algorithm == "sqp") {
    r <- Rsolnp::solnp(par0, obj, LB = bounds$lower, UB = bounds$upper,
                       control = utils::modifyList(list(trace = 0), control))
    list(par = r$pars, value = utils::tail(r$values, 1), convergence = r$convergence)
  } else {
    r <- stats::nlminb(par0, obj, lower = bounds$lower, upper = bounds$upper, control = control)
    list(par = r$par, value = r$objective, convergence = r$convergence)
  }
  par_unit <- stats::setNames(opt$par, names(par0))
  if (.neg_loglik(par_unit, xs, spec, diag = diag) >= .PENALTY)
    stop("The optimizer ended outside the parameter space: ", diag$reason, call. = FALSE)

  # Estimates and covariance on the original scale
  par <- .rescale_par(par_unit, spec, s, to_unit = FALSE)
  vc <- matrix(NA_real_, length(par), length(par), dimnames = list(names(par), names(par)))
  at_bound <- .at_bound(par_unit, bounds)
  if (hessian) {
    vc_unit <- .vcov_unit(par_unit, xs, spec, at_bound)
    J <- .rescale_jacobian(par_unit, spec, s)
    vc <- J %*% vc_unit %*% t(J)
  }

  model <- gs_model(spec, par = par)
  f <- .filter_model(x, .unpack(par, spec))
  fit <- list(call = match.call(), spec = spec, model = model, coef = par, vcov = vc,
              loglik = -(opt$value + N * log(s)), nobs = N, data = x,
              residuals = f$e, sigma = f$sigma, h = f$h,
              persistence = gs_stationarity(model),
              convergence = opt$convergence, algorithm = algorithm, at_bound = at_bound)
  class(fit) <- "gs_fit"

  if (opt$convergence != 0) warning("The optimizer did not report convergence (code ",
                                    opt$convergence, ").", call. = FALSE)
  if (fit$persistence >= 1)
    warning(sprintf("Persistence %.4f >= 1: E(sigma_t^delta) is not finite and long horizon scale forecasts diverge.",
                    fit$persistence), call. = FALSE)
  fit
}


# Parameters closer to a bound than a small fraction of the interval width
.at_bound <- function(par, bounds, tol = 1e-4)
{
  width <- bounds$upper - bounds$lower
  near <- par - bounds$lower < tol * width | bounds$upper - par < tol * width
  names(par)[near]
}

# Covariance of the unit scale estimates from the numerical Hessian. Steps that
# hit the penalty make the result meaningless, so they are detected and reported.
.vcov_unit <- function(par, xs, spec, at_bound)
{
  k <- length(par); nms <- names(par)
  na <- matrix(NA_real_, k, k, dimnames = list(nms, nms))
  hits <- new.env(); hits$reason <- NULL
  obj <- function(p) {
    v <- .neg_loglik(p, xs, spec, diag = hits)
    v
  }
  free <- setdiff(nms, at_bound)
  if (length(at_bound))
    warning("Parameter(s) at a bound, standard errors set to NA: ",
            paste(at_bound, collapse = ", "), ".", call. = FALSE)
  if (!length(free)) return(na)

  # Differentiate only in the free parameters, keeping bound ones fixed
  obj_free <- function(pf) { p <- par; p[free] <- pf; obj(p) }
  H <- tryCatch(.hessian(obj_free, par[free]), error = function(e) NULL)
  if (!is.null(hits$reason)) {
    warning("A differentiation step left the parameter space (", hits$reason,
            "); standard errors set to NA.", call. = FALSE)
    return(na)
  }
  V <- if (is.null(H)) NULL else tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V)) || any(diag(V) <= 0)) {
    warning("The Hessian is not invertible or not positive definite; standard errors set to NA.",
            call. = FALSE)
    return(na)
  }
  na[free, free] <- V
  na
}


# Central difference Hessian with steps adapted to the curvature: each step moves the
# log likelihood by about `target`. That keeps steps far above the small numerical noise
# of densities like the stable one, and inside the region where the surface is quadratic.
.hessian <- function(f, par, target = 0.05)
{
  k <- length(par); f0 <- f(par)
  second <- function(i, h) {
    e <- replace(numeric(k), i, h)
    (f(par + e) - 2 * f0 + f(par - e)) / h^2
  }
  # Pilot curvature per coordinate, then the adapted step
  h <- vapply(seq_len(k), function(i) {
    h0 <- 1e-3 * max(abs(par[i]), 1e-2)
    for (try in 1:6) {
      c0 <- second(i, h0)
      if (is.finite(c0) && c0 > 0) return(sqrt(2 * target / c0))
      h0 <- h0 * 4
    }
    h0
  }, numeric(1))
  H <- matrix(0, k, k)
  shift <- function(i, j, si, sj) { p <- par; p[i] <- p[i] + si * h[i]; p[j] <- p[j] + sj * h[j]; p }
  for (i in seq_len(k)) {
    H[i, i] <- second(i, h[i])
    for (j in seq_len(i - 1))
      H[i, j] <- H[j, i] <- (f(shift(i, j, 1, 1)) - f(shift(i, j, 1, -1)) -
                             f(shift(i, j, -1, 1)) + f(shift(i, j, -1, -1))) / (4 * h[i] * h[j])
  }
  H
}
