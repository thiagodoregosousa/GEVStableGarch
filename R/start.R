# Start values and optimization bounds. Both live on the internal scale where
# the data has unit standard deviation (see .rescale_par).

.opt_bounds <- function(spec, x)
{
  o <- spec$order; nms <- .par_names(spec)
  lower <- upper <- stats::setNames(numeric(length(nms)), nms)
  set <- function(pattern, lo, up) {
    i <- grepl(pattern, nms)
    lower[i] <<- lo; upper[i] <<- up
  }
  set("^mu$", mean(x) - 10, mean(x) + 10)
  set("^(ar|ma)[0-9]+$", -3, 3)
  set("^omega$", 1e-6, 10)
  set("^alpha[0-9]+$", 1e-8, 1)
  set("^gamma[0-9]+$", -0.999, 0.999)
  set("^beta[0-9]+$", 0, 0.9999)
  set("^delta$", 0.1, 4)
  lower[spec$dist$par_names] <- spec$dist$lower
  upper[spec$dist$par_names] <- spec$dist$upper
  list(lower = lower, upper = upper)
}

.start_values <- function(spec, x, bounds)
{
  o <- spec$order; dist <- spec$dist
  m <- o[["m"]]; n <- o[["n"]]; p <- o[["p"]]; q <- o[["q"]]

  # ARMA part from a quick CSS fit; arima reports the mean, the model uses the intercept
  ar <- rep(0, m); ma <- rep(0, n); mu <- mean(x)
  if (m + n > 0) {
    a <- tryCatch(stats::arima(x, order = c(m, 0, n), method = "CSS",
                               include.mean = spec$include_mean)$coef, error = function(e) NULL)
    if (!is.null(a)) {
      ar <- unname(a[seq_len(m)]); ma <- unname(a[m + seq_len(n)])
      if (spec$include_mean) mu <- a[["intercept"]] * (1 - sum(ar))
    }
  }

  # Distribution start from the data when the family provides it
  z <- (x - stats::median(x)) / stats::mad(x)
  dist_par <- dist$start
  if (!is.null(dist$start_fun)) {
    s <- tryCatch(dist$start_fun(z), error = function(e) NULL)
    if (!is.null(s) && .dist_par_ok(dist, s)) dist_par <- s[dist$par_names]
  }

  # Power below the moment bound, then omega matching the level of |e|^delta at persistence ~0.9
  delta <- dist$default_delta
  if (spec$aparch && delta >= 0.95 * dist$max_power(dist_par)) delta <- 0.5 * dist$max_power(dist_par)
  alpha <- rep(if (q > 0) 0.1 / p else 0.3 / p, p)
  beta <- rep(0.8 / max(q, 1), q)
  kappa <- gs_aparch_moment(dist, dist_par, 0, delta)
  persistence <- sum(alpha) * kappa + sum(beta)
  if (persistence >= 0.95) {
    alpha <- alpha * 0.5 / (sum(alpha) * kappa)
    persistence <- 0.5 + sum(beta)
  }
  level <- mean(abs(x - mean(x))^delta) / kappa
  omega <- max(level * (1 - min(persistence, 0.95)), 1e-4)

  par <- c(if (spec$include_mean) mu, ar, ma, omega, alpha,
           if (spec$aparch) rep(0, p), beta, if (spec$aparch) delta, dist_par)
  names(par) <- .par_names(spec)
  .clamp(par, bounds)
}

# Keep a parameter vector strictly inside the box
.clamp <- function(par, bounds, margin = 1e-4)
{
  width <- bounds$upper - bounds$lower
  pmin(pmax(par, bounds$lower + margin * width), bounds$upper - margin * width)
}

# Map parameters between the original and the unit variance scale:
# x' = x / s implies mu' = mu / s and omega' = omega / s^delta
.rescale_par <- function(par, spec, s, to_unit = TRUE)
{
  delta <- .unpack(par, spec)$delta
  f <- if (to_unit) 1 / s else s
  if (spec$include_mean) par["mu"] <- par["mu"] * f
  par["omega"] <- par["omega"] * f^delta
  par
}

# Jacobian of the map from unit scale parameters to original parameters
.rescale_jacobian <- function(par_unit, spec, s)
{
  k <- length(par_unit); J <- diag(k); nms <- names(par_unit)
  dimnames(J) <- list(nms, nms)
  delta <- .unpack(par_unit, spec)$delta
  if (spec$include_mean) J["mu", "mu"] <- s
  J["omega", "omega"] <- s^delta
  if (spec$aparch) J["omega", "delta"] <- par_unit[["omega"]] * s^delta * log(s)
  J
}
