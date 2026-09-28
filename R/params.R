# Single place that knows the layout of the parameter vector. Likelihood,
# bounds, start values, fitting output and simulation all go through here.

.par_names <- function(spec)
{
  o <- spec$order
  seq_named <- function(prefix, k) if (k > 0) paste0(prefix, seq_len(k)) else character(0)
  c(if (spec$include_mean) "mu",
    seq_named("ar", o[["m"]]), seq_named("ma", o[["n"]]),
    "omega", seq_named("alpha", o[["p"]]),
    if (spec$aparch) seq_named("gamma", o[["p"]]),
    seq_named("beta", o[["q"]]),
    if (spec$aparch) "delta",
    spec$dist$par_names)
}

# Named vector to a structured list, filling the fixed values of GARCH mode
.unpack <- function(par, spec)
{
  o <- spec$order
  get <- function(prefix, k) if (k > 0) unname(par[paste0(prefix, seq_len(k))]) else numeric(0)
  list(mu = if (spec$include_mean) unname(par[["mu"]]) else 0,
       ar = get("ar", o[["m"]]), ma = get("ma", o[["n"]]),
       omega = unname(par[["omega"]]), alpha = get("alpha", o[["p"]]),
       gamma = if (spec$aparch) get("gamma", o[["p"]]) else rep(0, o[["p"]]),
       beta = get("beta", o[["q"]]),
       delta = if (spec$aparch) unname(par[["delta"]]) else spec$dist$default_delta,
       dist_par = par[spec$dist$par_names])
}

# Reason why unpacked parameters are outside the model's parameter space, NULL if admissible
.inadmissible <- function(u, spec)
{
  dist <- spec$dist
  if (any(!is.finite(unlist(u[c("mu", "ar", "ma", "omega", "alpha", "gamma", "beta", "delta")]))))
    return("non finite parameter")
  if (u$omega <= 0) return("omega must be positive")
  if (any(u$alpha < 0) || any(u$beta < 0)) return("alpha and beta must be nonnegative")
  if (any(abs(u$gamma) >= 1)) return("gamma must lie in (-1, 1)")
  if (u$delta <= 0) return("delta must be positive")
  if (length(u$ar) && !.ar_stationary(u$ar)) return("AR polynomial not stationary")
  if (!.dist_par_ok(dist, u$dist_par))
    return(sprintf("%s parameters outside their bounds or not valid", dist$name))
  mp <- dist$max_power(u$dist_par)
  if (u$delta >= mp)
    return(sprintf("delta (%s) not below max_power (%s): E|z|^delta is infinite",
                   format(u$delta), format(mp)))
  NULL
}

.ar_stationary <- function(ar) all(Mod(polyroot(c(1, -ar))) > 1)
