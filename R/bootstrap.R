#' Parametric bootstrap standard errors and confidence intervals
#'
#' Simulates `n_boot` series of the same length from the fitted model, refits
#' each one, and summarizes the resulting sampling distribution of the
#' estimates. This is the recommended alternative to the Hessian standard
#' errors of [gs_fit()] when the asymptotic normal approximation is unreliable:
#' near a boundary (stable tail index near 2, GEV shape near \eqn{\pm 0.5},
#' APARCH asymmetry near \eqn{\pm 1}) or, for GEV innovations, because the
#' support depends on the shape parameter and the usual regularity conditions
#' do not hold (Efron, 1979).
#'
#' Each replication is refit with `hessian = FALSE`; replications whose fit
#' fails are dropped and counted in `n_fail`.
#'
#' @param fit A `gs_fit` object.
#' @param n_boot Number of bootstrap replications.
#' @param level Confidence level of the percentile intervals.
#' @param seed Optional seed for reproducibility.
#' @param ... Passed to [gs_fit()] (for example `algorithm`).
#' @return A list with `se` (named vector of bootstrap standard errors), `ci`
#'   (a matrix with `lower` and `upper` percentile limits), `replicates` (the
#'   `n_boot` by number-of-parameters matrix of estimates, failed fits as `NA`
#'   rows), `n_fail` and `level`.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'   jackknife. The Annals of Statistics, 7(1), 1-26.
#' @seealso [gs_fit()]
#' @examples
#' \donttest{
#' spec <- gs_spec(garch = c(1, 1), dist = gs_gev())
#' model <- gs_model(spec, mu = 0, omega = 0.1, alpha = 0.1, beta = 0.8,
#'                   dist_par = c(xi = 0.2))
#' x <- gs_sim(model, n = 800, seed = 1)$y
#' fit <- gs_fit(x, spec)
#' gs_bootstrap(fit, n_boot = 99, seed = 1)
#' }
#' @export
gs_bootstrap <- function(fit, n_boot = 199, level = 0.95, seed = NULL, ...)
{
  if (!inherits(fit, "gs_fit")) stop("`fit` must be a gs_fit object.", call. = FALSE)
  if (n_boot < 2) stop("`n_boot` must be at least 2.", call. = FALSE)
  if (level <= 0 || level >= 1) stop("`level` must be in (0, 1).", call. = FALSE)
  spec <- fit$spec; model <- fit$model; n <- fit$nobs; nms <- names(fit$coef)
  if (!is.null(seed)) set.seed(seed)

  reps <- matrix(NA_real_, n_boot, length(nms), dimnames = list(NULL, nms))
  for (b in seq_len(n_boot)) {
    y <- suppressWarnings(gs_sim(model, n = n)$y)
    est <- tryCatch(suppressWarnings(stats::coef(gs_fit(y, spec, hessian = FALSE, ...))),
                    error = function(e) NULL)
    if (!is.null(est)) reps[b, ] <- est[nms]
  }

  ok <- stats::complete.cases(reps)
  if (sum(ok) < 2)
    stop("Too few bootstrap fits succeeded (", sum(ok), "); cannot summarize.", call. = FALSE)
  if (any(!ok))
    warning(sum(!ok), " of ", n_boot, " bootstrap fits failed and were dropped.", call. = FALSE)

  good <- reps[ok, , drop = FALSE]
  a <- (1 - level) / 2
  ci <- t(apply(good, 2, stats::quantile, probs = c(a, 1 - a), names = FALSE))
  dimnames(ci) <- list(nms, c("lower", "upper"))
  list(se = apply(good, 2, stats::sd), ci = ci, replicates = reps,
       n_fail = sum(!ok), level = level)
}
