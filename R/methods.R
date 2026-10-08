#' Methods for fitted models
#'
#' Accessors for `gs_fit` objects. `residuals()` returns, by `type`:
#' `"raw"` the residuals \eqn{e_t}, `"standardized"` the innovations
#' \eqn{e_t / \sigma_t}, `"pit"` their probability integral transform
#' \eqn{F(e_t / \sigma_t)} under the fitted innovation distribution (uniform
#' under a correct model; the input of copula models and of KS diagnostics).
#' `sigma()` returns the conditional scale \eqn{\sigma_t} (not a standard
#' deviation when the innovation variance is infinite), `fitted()` the
#' conditional location \eqn{x_t - e_t}.
#'
#' @param object,x A `gs_fit` object.
#' @param type Type of residuals.
#' @param ... Unused.
#' @name gs_fit-methods
NULL

#' @rdname gs_fit-methods
#' @export
coef.gs_fit <- function(object, ...) object$coef

#' @rdname gs_fit-methods
#' @export
vcov.gs_fit <- function(object, ...) object$vcov

#' @rdname gs_fit-methods
#' @export
logLik.gs_fit <- function(object, ...)
  structure(object$loglik, df = length(object$coef), nobs = object$nobs, class = "logLik")

#' @rdname gs_fit-methods
#' @export
nobs.gs_fit <- function(object, ...) object$nobs

#' @rdname gs_fit-methods
#' @export
residuals.gs_fit <- function(object, type = c("raw", "standardized", "pit"), ...)
{
  type <- match.arg(type)
  if (type == "raw") return(object$residuals)
  z <- object$residuals / object$sigma
  if (type == "standardized") return(z)
  object$spec$dist$cdf(z, .unpack(object$coef, object$spec)$dist_par)
}

#' @rdname gs_fit-methods
#' @export
sigma.gs_fit <- function(object, ...) object$sigma

#' @rdname gs_fit-methods
#' @export
fitted.gs_fit <- function(object, ...) object$data - object$residuals

#' @rdname gs_fit-methods
#' @export
print.gs_fit <- function(x, ...)
{
  cat(.model_label(x$spec), "with", x$spec$dist$name, "innovations\n")
  se <- sqrt(diag(x$vcov))
  tval <- x$coef / se
  tab <- cbind(Estimate = x$coef, `Std. Error` = se, `t value` = tval,
               `Pr(>|t|)` = 2 * stats::pnorm(-abs(tval)))
  stats::printCoefmat(tab, na.print = "NA")
  cat(sprintf("\nLog likelihood: %.4f   AIC: %.4f   BIC: %.4f   n = %d\n",
              x$loglik, stats::AIC(x), stats::BIC(x), x$nobs))
  cat(sprintf("Persistence: %.4f%s\n", x$persistence,
              if (x$persistence >= 1) "  (not stationary in the sense of finite E sigma^delta)" else ""))
  if (length(x$at_bound)) cat("At a bound:", paste(x$at_bound, collapse = ", "), "\n")
  d <- x$diagnostics
  if (!is.null(d))
    cat(sprintf("Convergence: scaled gradient %.1e, Hessian rcond %s%s\n",
                d$grad_rel, if (is.na(d$hess_rcond)) "NA" else sprintf("%.1e", d$hess_rcond),
                if (isTRUE(d$cap_hit)) ", iteration cap hit" else ""))
  invisible(x)
}

#' @rdname gs_fit-methods
#' @export
summary.gs_fit <- function(object, ...) print.gs_fit(object, ...)
