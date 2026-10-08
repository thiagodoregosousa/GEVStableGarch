# Base value returned outside the parameter space. Large but finite so numerical
# optimizers and derivatives do not break on Inf. A breach adds `.PENALTY *
# violation` on top, so the penalty slopes back toward the feasible region
# rather than forming a flat plateau the optimizer can stall on.
.PENALTY <- 1e10

# Negative log likelihood of a named parameter vector. When `diag` is an
# environment, the reason of a penalty is recorded in diag$reason, so the fit
# can report which condition was hit instead of failing silently.
.neg_loglik <- function(par, x, spec, h_init = NULL, diag = NULL)
{
  penalty <- function(why, violation = 0) {
    if (!is.null(diag)) diag$reason <- why
    .PENALTY * (1 + violation)
  }
  names(par) <- .par_names(spec)
  u <- .unpack(par, spec)
  why <- .inadmissible(u, spec)
  if (!is.null(why)) return(penalty(why, .violation(u, spec)))

  f <- .filter_model(x, u, h_init)
  if (any(!is.finite(f$h)) || any(f$h <= 0)) return(penalty("non positive or non finite scale"))
  ld <- spec$dist$log_density(f$e / f$sigma, u$dist_par)
  if (any(!is.finite(ld))) return(penalty("log density not finite (observation outside the support?)"))
  -sum(ld - log(f$sigma))
}
