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
#' These Hessian standard errors are asymptotic and rely on the usual
#' regularity conditions. They become unreliable *near* a boundary of the
#' parameter space, even when not exactly on it: the stable tail index
#' approaching 2 (where the skewness is no longer identified), the GEV shape
#' approaching \eqn{\pm 0.5}, or an APARCH asymmetry approaching \eqn{\pm 1}.
#' The GEV support depends on its shape parameter, so for GEV innovations the
#' usual maximum likelihood confidence intervals are not guaranteed valid;
#' prefer the parametric bootstrap of [gs_bootstrap()] there and in the
#' near-boundary cases above. The fit records convergence diagnostics in
#' `$diagnostics` (scaled gradient, Hessian reciprocal condition number and
#' definiteness, optimizer iterations), which flag these situations.
#'
#' @param data Numeric vector, the time series (for example log returns).
#' @param spec A `gs_spec` object.
#' @param algorithm `"sqp"` (sequential quadratic programming via
#'   [Rsolnp::solnp()], Galanos and Ye 2025) or `"nlminb"`.
#' @param start Optional named vector of start values on the original scale,
#'   for example `coef(previous_fit)` in rolling windows. Missing names are
#'   filled with the default start values.
#' @param hessian Logical, compute standard errors. Set `FALSE` to save time
#'   when only point estimates are needed.
#' @param restarts Integer, number of additional optimizer runs from perturbed
#'   start values; the fit with the highest log likelihood is kept. A few
#'   restarts (for example 3 to 5) make the estimate robust to the local optima
#'   of multimodal surfaces, such as ARMA-GARCH models with heavy-tailed
#'   innovations, at a proportional increase in time. The default 0 keeps the
#'   single run from the data-driven start.
#' @param control List passed to the optimizer.
#' @return An object of class `gs_fit` with methods [coef()], [vcov()],
#'   [logLik()], [residuals()], [sigma()], [fitted()], [predict()] and
#'   `print()`.
#' @references Galanos, A. and Ye, Y. (2025). Rsolnp: General Non-Linear
#'   Optimization. R package version 2.0.1.
#'   \url{https://CRAN.R-project.org/package=Rsolnp}
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
                   hessian = TRUE, restarts = 0L, control = list())
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
  if (.neg_loglik(par0, xs, spec, diag = diag) >= .PENALTY) {
    # A bounded-support innovation (for example GEV with a large shape, or a
    # user distribution) can leave the data outside the support when the
    # starting scale is too small. Enlarge omega toward its bound before failing.
    cap <- bounds$upper[["omega"]] * (1 - 1e-4)
    while (par0[["omega"]] < cap && .neg_loglik(par0, xs, spec, diag = diag) >= .PENALTY)
      par0["omega"] <- min(par0[["omega"]] * 3, cap)
    if (.neg_loglik(par0, xs, spec, diag = diag) >= .PENALTY)
      stop("Start values are inadmissible: ", diag$reason, call. = FALSE)
  }

  # Optimize on the unit scale, with explicit caps so unattended runs terminate.
  # With `restarts > 0`, re-optimize from perturbed start values and keep the
  # best log likelihood, to escape the local optima of multimodal surfaces
  # (for example ARMA-GARCH with heavy-tailed innovations).
  obj <- function(p) .neg_loglik(p, xs, spec)
  opt <- .optimize_once(par0, obj, bounds, algorithm, control)
  for (r in seq_len(restarts)) {
    start_r <- .perturb_start(par0, bounds)
    if (.neg_loglik(start_r, xs, spec) >= .PENALTY) next
    opt_r <- tryCatch(.optimize_once(start_r, obj, bounds, algorithm, control),
                      error = function(e) NULL)
    if (!is.null(opt_r) && is.finite(opt_r$value) && opt_r$value < opt$value) opt <- opt_r
  }
  par_unit <- stats::setNames(opt$par, names(par0))
  if (.neg_loglik(par_unit, xs, spec, diag = diag) >= .PENALTY)
    stop("The optimizer ended outside the parameter space: ", diag$reason, call. = FALSE)

  # Estimates, gradient/conditioning verdict and covariance on the original scale
  par <- .rescale_par(par_unit, spec, s, to_unit = FALSE)
  at_bound <- .at_bound(par_unit, bounds)
  inf <- .inference(par_unit, xs, spec, at_bound, hessian)
  vc <- if (hessian) {
    J <- .rescale_jacobian(par_unit, spec, s); J %*% inf$vcov %*% t(J)
  } else inf$vcov
  cap_hit <- !is.null(opt$iterations) && !is.null(opt$cap) && opt$iterations >= opt$cap

  model <- gs_model(spec, par = par)
  f <- .filter_model(x, .unpack(par, spec))
  fit <- list(call = match.call(), spec = spec, model = model, coef = par, vcov = vc,
              loglik = -(opt$value + N * log(s)), nobs = N, data = x,
              residuals = f$e, sigma = f$sigma, h = f$h,
              persistence = gs_stationarity(model),
              convergence = opt$convergence, algorithm = algorithm, at_bound = at_bound,
              diagnostics = list(grad_max = inf$grad_max, grad_rel = inf$grad_rel,
                                 hess_rcond = inf$hess_rcond, hess_pd = inf$hess_pd,
                                 iterations = opt$iterations, cap_hit = cap_hit))
  class(fit) <- "gs_fit"

  if (opt$convergence != 0) warning("The optimizer did not report convergence (code ",
                                    opt$convergence, ").", call. = FALSE)
  if (cap_hit) warning("The optimizer hit its iteration cap (", opt$iterations,
                       "); the solution may not be converged.", call. = FALSE)
  if (isTRUE(inf$grad_rel > 1e-2))
    warning(sprintf("The gradient is not close to zero at the solution (scaled max %.1e); the optimizer may not have converged.",
                    inf$grad_rel), call. = FALSE)
  if (fit$persistence >= 1)
    warning(sprintf("Persistence %.4f >= 1: E(sigma_t^delta) is not finite and long horizon scale forecasts diverge.",
                    fit$persistence), call. = FALSE)
  fit
}


# One optimizer run from a given (unit scale) start, returning a common result shape
.optimize_once <- function(par_start, obj, bounds, algorithm, control)
{
  if (algorithm == "sqp") {
    ctrl <- utils::modifyList(list(trace = 0, outer.iter = 400, inner.iter = 800), control)
    r <- Rsolnp::solnp(par_start, obj, LB = bounds$lower, UB = bounds$upper, control = ctrl)
    list(par = r$pars, value = utils::tail(r$values, 1), convergence = r$convergence,
         iterations = r$outer.iter, cap = ctrl$outer.iter)
  } else {
    ctrl <- utils::modifyList(list(iter.max = 500, eval.max = 1000), control)
    r <- stats::nlminb(par_start, obj, lower = bounds$lower, upper = bounds$upper, control = ctrl)
    list(par = r$par, value = r$objective, convergence = r$convergence,
         iterations = r$iterations, cap = ctrl$iter.max)
  }
}

# A jittered start for a restart: each coordinate perturbed by about `sd` of its
# magnitude (or of the box width when near zero), kept strictly inside the box.
.perturb_start <- function(par0, bounds, sd = 0.3)
{
  width <- bounds$upper - bounds$lower
  step <- stats::rnorm(length(par0), 0, sd) * pmax(abs(par0), 0.05 * width)
  .clamp(par0 + step, bounds)
}

# Parameters closer to a bound than a small fraction of the interval width
.at_bound <- function(par, bounds, tol = 1e-4)
{
  width <- bounds$upper - bounds$lower
  near <- par - bounds$lower < tol * width | bounds$upper - par < tol * width
  names(par)[near]
}

# Gradient, Hessian-based covariance and conditioning of the unit scale estimates.
# Returns the covariance (NA where a parameter is on a bound or the Hessian is
# unusable), the scaled gradient at the solution (a first-order convergence
# check), and the Hessian's reciprocal condition number and definiteness. Steps
# that hit the penalty make the Hessian meaningless, so they are detected.
.inference <- function(par, xs, spec, at_bound, hessian)
{
  k <- length(par); nms <- names(par)
  na <- matrix(NA_real_, k, k, dimnames = list(nms, nms))
  free <- setdiff(nms, at_bound)
  obj_nodiag <- function(p) .neg_loglik(p, xs, spec)
  f0 <- obj_nodiag(par)

  # Scaled gradient over the free parameters: |g_i| |par_i| / |nll|, dimensionless
  g <- if (length(free))
    .gradient(function(pf) { p <- par; p[free] <- pf; obj_nodiag(p) }, par[free]) else numeric(0)
  grad_max <- if (length(g)) max(abs(g)) else 0
  grad_rel <- if (length(g)) max(abs(g) * pmax(abs(par[free]), 1)) / max(abs(f0), 1) else 0
  out <- list(vcov = na, grad_max = grad_max, grad_rel = grad_rel,
              hess_rcond = NA_real_, hess_pd = NA)

  if (length(at_bound))
    warning("Parameter(s) at a bound, standard errors set to NA: ",
            paste(at_bound, collapse = ", "), ".", call. = FALSE)
  if (!hessian || !length(free)) return(out)

  # Differentiate only in the free parameters, keeping bound ones fixed
  hits <- new.env(); hits$reason <- NULL
  obj_free <- function(pf) { p <- par; p[free] <- pf; .neg_loglik(p, xs, spec, diag = hits) }
  H <- tryCatch(.hessian(obj_free, par[free]), error = function(e) NULL)
  if (!is.null(hits$reason)) {
    warning("A differentiation step left the parameter space (", hits$reason,
            "); standard errors set to NA.", call. = FALSE)
    return(out)
  }
  if (!is.null(H)) {
    ev <- tryCatch(eigen(H, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
    if (!is.null(ev)) {
      out$hess_pd <- min(ev) > 0
      out$hess_rcond <- min(abs(ev)) / max(abs(ev))
    }
  }
  V <- if (is.null(H)) NULL else tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(V)) || !isTRUE(out$hess_pd)) {
    warning("The Hessian is not invertible or not positive definite; standard errors set to NA.",
            call. = FALSE)
    return(out)
  }
  if (isTRUE(out$hess_rcond < 1e-8))
    warning(sprintf("The Hessian is ill conditioned (reciprocal condition number %.1e); standard errors may be unreliable.",
                    out$hess_rcond), call. = FALSE)
  na[free, free] <- V
  out$vcov <- na
  out
}

# Central difference gradient with a small relative step
.gradient <- function(f, par)
{
  k <- length(par); h <- 1e-4 * pmax(abs(par), 1e-2)
  vapply(seq_len(k), function(i) {
    e <- replace(numeric(k), i, h[i])
    (f(par + e) - f(par - e)) / (2 * h[i])
  }, numeric(1))
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
