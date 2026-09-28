#' Innovation distribution for ARMA-GARCH/APARCH models
#'
#' Builds the object that describes the innovation distribution \eqn{z_t} of
#' the model. The built in families [gs_stable()], [gs_gev()] and [gs_gat()]
#' are created with this constructor, and a user defined distribution gets
#' exactly the same fitting, simulation, stationarity and forecasting
#' machinery.
#'
#' All functions work on the standardized innovation (location 0, scale 1).
#' `z` is a numeric vector and `par` a numeric vector named by `par_names`.
#'
#' Every field is required except `valid`, `aparch_moment` and `start_fun`,
#' which have well defined fallbacks. In particular `max_power` and `mean`
#' have no defaults: the user must state when moments exist, so that fitting
#' never fails silently. Unless `check = FALSE`, [gs_check_dist()] runs a
#' self check at `start` and stops with an error naming the failing field.
#'
#' @param name Character, name of the distribution.
#' @param log_density `function(z, par)` returning \eqn{\log f(z)}, vectorized
#'   in `z`, and `-Inf` outside the support.
#' @param random `function(n, par)` returning `n` draws.
#' @param cdf `function(q, par)`, distribution function.
#' @param quantile `function(p, par)`, quantile function.
#' @param par_names Character vector with the parameter names. They must not
#'   clash with model parameter names (`mu`, `ar1`, `omega`, `alpha1`, ...).
#' @param lower,upper,start Numeric vectors named by `par_names`: box bounds
#'   used in estimation and the default start value (strictly inside).
#' @param default_delta Power \eqn{\delta} used when the model is a GARCH
#'   (not APARCH) model, e.g. 2 for a variance recursion, 1 for a scale
#'   recursion.
#' @param max_power `function(par)` returning the supremum of the powers
#'   \eqn{\delta} with \eqn{E|z|^\delta < \infty} (`Inf` if all moments exist).
#' @param mean `function(par)` returning \eqn{E[z]}; it must return `NA`
#'   explicitly where the mean does not exist.
#' @param valid Optional `function(par)` returning `TRUE`/`FALSE` for
#'   constraints that are not box bounds.
#' @param aparch_moment Optional `function(gamma, delta, par)` returning
#'   \eqn{E(|z| - \gamma z)^\delta} in closed form. It may return `NA` where no
#'   closed form is available; numerical integration is then used.
#' @param start_fun Optional `function(z)` returning data driven start values
#'   from standardized data.
#' @param check Logical, run [gs_check_dist()] at `start`.
#'
#' @return An object of class `gs_dist`.
#' @seealso [gs_check_dist()], [gs_stable()], [gs_gev()], [gs_gat()]
#' @examples
#' # Normal innovations written as a user defined distribution
#' norm <- gs_dist(
#'   name = "norm",
#'   log_density = function(z, par) dnorm(z, log = TRUE),
#'   random = function(n, par) rnorm(n),
#'   cdf = function(q, par) pnorm(q),
#'   quantile = function(p, par) qnorm(p),
#'   par_names = character(0),
#'   lower = numeric(0), upper = numeric(0), start = numeric(0),
#'   default_delta = 2,
#'   max_power = function(par) Inf,
#'   mean = function(par) 0)
#' norm
#' @export
gs_dist <- function(name, log_density, random, cdf, quantile,
                    par_names, lower, upper, start,
                    default_delta, max_power, mean,
                    valid = NULL, aparch_moment = NULL, start_fun = NULL,
                    check = TRUE)
{
  # Required fields have no defaults: report every missing one at once
  required <- c("name", "log_density", "random", "cdf", "quantile", "par_names",
                "lower", "upper", "start", "default_delta", "max_power", "mean")
  frame <- environment()
  missing_fields <- required[vapply(required, function(f) eval(call("missing", as.name(f)), frame), logical(1))]
  if (length(missing_fields))
    stop("gs_dist: missing required field(s): ", paste0("`", missing_fields, "`", collapse = ", "),
         call. = FALSE)

  dist <- list(name = name, log_density = log_density, random = random, cdf = cdf,
               quantile = quantile, par_names = as.character(par_names),
               lower = .named(lower, par_names), upper = .named(upper, par_names),
               start = .named(start, par_names), default_delta = default_delta,
               max_power = max_power, mean = mean, valid = valid,
               aparch_moment = aparch_moment, start_fun = start_fun)
  class(dist) <- "gs_dist"

  .check_dist_structure(dist)
  if (check) gs_check_dist(dist)
  dist
}


#' Self check of an innovation distribution
#'
#' Verifies that the pieces of a [gs_dist()] object are consistent at a
#' parameter value: names and bounds, finite and vectorized log density,
#' density integrating to one, `quantile` inverting `cdf`, `cdf` matching the
#' density, `random` matching `cdf` (Kolmogorov-Smirnov), moment conditions
#' (`max_power` above `default_delta`, `mean` and `aparch_moment` against
#' numerical integration). It stops with an error naming the failing field.
#'
#' @param dist A `gs_dist` object.
#' @param par Parameter value where the checks are done (default `dist$start`).
#' @param tol Tolerance for the numerical comparisons.
#' @return `dist`, invisibly.
#' @export
gs_check_dist <- function(dist, par = dist$start, tol = 1e-3)
{
  .check_dist_structure(dist)
  fail <- function(field, msg)
    stop(sprintf("gs_dist '%s': `%s` %s", dist$name, field, msg), call. = FALSE)
  par <- .named(par, dist$par_names)
  if (any(par <= dist$lower | par >= dist$upper))
    fail("start", "must lie strictly inside (lower, upper)")
  if (!is.null(dist$valid) && !isTRUE(dist$valid(par)))
    fail("valid", "returns FALSE at the check parameters")

  # Quantiles, density and distribution function on a central grid
  p <- c(0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999)
  q <- dist$quantile(p, par)
  if (length(q) != length(p) || any(!is.finite(q)) || is.unsorted(q))
    fail("quantile", "must return finite increasing values, vectorized in p")
  ld <- dist$log_density(q, par)
  if (length(ld) != length(q) || any(!is.finite(ld)))
    fail("log_density", "must return finite values, vectorized in z, inside the support")
  if (max(abs(dist$cdf(q, par) - p)) > tol)
    fail("quantile", "is not the inverse of `cdf`")
  h <- 1e-4 * pmax(1, abs(q))
  num_density <- (dist$cdf(q + h, par) - dist$cdf(q - h, par)) / (2 * h)
  if (max(abs(num_density / exp(ld) - 1)) > 10 * tol)
    fail("cdf", "is not consistent with `log_density`")

  # Density integrates to one
  total <- tryCatch(.integrate_dist(function(x) 1, dist, par), error = function(e) NA)
  if (!is.finite(total) || abs(total - 1) > tol)
    fail("log_density", sprintf("does not integrate to 1 (got %s)", format(total)))

  # Random generator against the distribution function, without touching the user seed
  draws <- .with_seed(20260928, dist$random(2000, par))
  if (length(draws) != 2000 || any(!is.finite(draws)))
    fail("random", "must return n finite draws")
  if (stats::ks.test(draws, function(x) dist$cdf(x, par))$p.value < 1e-4)
    fail("random", "draws do not follow `cdf` (Kolmogorov-Smirnov test)")

  # Moment declarations
  mp <- dist$max_power(par)
  if (!is.numeric(mp) || length(mp) != 1 || is.na(mp) || mp <= 0)
    fail("max_power", "must return a single positive number (Inf allowed)")
  if (mp <= dist$default_delta)
    fail("default_delta", sprintf("(%s) must be below `max_power` (%s)", dist$default_delta, mp))
  m <- dist$mean(par)
  if (length(m) != 1 || !(is.numeric(m) || is.na(m)))
    fail("mean", "must return a single number, or NA when the mean does not exist")
  if (mp <= 1 && !is.na(m))
    fail("mean", "must return NA when `max_power` <= 1 (the mean does not exist)")
  if (!is.na(m)) {
    num_mean <- tryCatch(.integrate_dist(identity, dist, par), error = function(e) NA)
    if (!is.finite(num_mean) || abs(num_mean - m) > 10 * tol * max(1, abs(m)))
      fail("mean", sprintf("(%s) disagrees with numerical integration (%s)", format(m), format(num_mean)))
  }
  if (!is.null(dist$aparch_moment)) {
    closed <- dist$aparch_moment(0.2, dist$default_delta, par)
    if (!is.na(closed)) {
      num <- .aparch_moment_numeric(dist, par, 0.2, dist$default_delta)
      if (abs(closed / num - 1) > tol)
        fail("aparch_moment", sprintf("(%s) disagrees with numerical integration (%s)", format(closed), format(num)))
    }
  }
  invisible(dist)
}


#' @export
print.gs_dist <- function(x, ...)
{
  cat("Innovation distribution:", x$name, "\n")
  if (length(x$par_names)) {
    tab <- rbind(start = x$start, lower = x$lower, upper = x$upper)
    print(tab)
  } else cat("No parameters\n")
  cat("GARCH mode delta:", x$default_delta, "\n")
  invisible(x)
}


# ------------------------------------------------------------------------------
# Internal helpers

# Names that the model itself uses; distribution parameters must avoid them
.reserved_name <- function(x) grepl("^(mu|omega|delta)$|^(ar|ma|alpha|gamma|beta)[0-9]+$", x)

.named <- function(x, nms)
{
  x <- as.numeric(x)
  if (length(x) == length(nms)) names(x) <- nms
  x
}

# Structural checks that are cheap and always run
.check_dist_structure <- function(dist)
{
  if (!inherits(dist, "gs_dist")) stop("`dist` must be a gs_dist object.", call. = FALSE)
  fail <- function(field, msg)
    stop(sprintf("gs_dist '%s': `%s` %s", dist$name, field, msg), call. = FALSE)
  for (f in c("log_density", "random", "cdf", "quantile", "max_power", "mean"))
    if (!is.function(dist[[f]])) fail(f, "must be a function")
  for (f in c("valid", "aparch_moment", "start_fun"))
    if (!is.null(dist[[f]]) && !is.function(dist[[f]])) fail(f, "must be a function or NULL")
  k <- length(dist$par_names)
  if (anyDuplicated(dist$par_names)) fail("par_names", "must be unique")
  if (any(.reserved_name(dist$par_names)))
    fail("par_names", "clash with model parameter names (mu, omega, delta, ar1, alpha1, ...)")
  for (f in c("lower", "upper", "start"))
    if (length(dist[[f]]) != k || any(is.na(dist[[f]])))
      fail(f, "must be a numeric vector with one value per `par_names`")
  if (k && any(dist$lower >= dist$upper)) fail("lower", "must be below `upper`")
  if (k && any(dist$start <= dist$lower | dist$start >= dist$upper))
    fail("start", "must lie strictly inside (lower, upper)")
  if (!is.numeric(dist$default_delta) || length(dist$default_delta) != 1 || dist$default_delta <= 0)
    fail("default_delta", "must be a single positive number")
  invisible(TRUE)
}

# Integral of g(x) f(x) over the support, split at central quantiles for accuracy in the tails
.integrate_dist <- function(g, dist, par)
{
  f <- function(x) g(x) * exp(dist$log_density(x, par))
  cuts <- sort(unique(c(dist$quantile(c(0.001, 0.5, 0.999), par), 0)))
  bounds <- c(-Inf, cuts, Inf)
  # Heavy tails can trigger roundoff detection at tight tolerance; relax once before giving up
  piece <- function(a, b)
    tryCatch(stats::integrate(f, a, b, subdivisions = 2000L, rel.tol = 1e-9)$value,
             error = function(e) stats::integrate(f, a, b, subdivisions = 5000L, rel.tol = 1e-6,
                                                  stop.on.error = FALSE)$value)
  sum(vapply(seq_len(length(bounds) - 1), function(i) piece(bounds[i], bounds[i + 1]), numeric(1)))
}

# Evaluate an expression with a fixed seed and restore the caller's RNG state
.with_seed <- function(seed, expr)
{
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old, envir = globalenv())
  })
  set.seed(seed)
  expr
}

# Parameter value inside bounds and valid
.dist_par_ok <- function(dist, par)
{
  all(par > dist$lower & par < dist$upper) && (is.null(dist$valid) || isTRUE(dist$valid(par)))
}
