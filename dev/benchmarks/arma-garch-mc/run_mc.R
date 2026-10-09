# Replication driver for the Monte Carlo study. Parallel over replications with
# the base `parallel` package, capped at 5 cores so the machine stays usable.
# Each replication is seeded independently (base_seed + r), so results are
# reproducible regardless of how the forks are scheduled. See MC_PLAN.md.
#
# Requires GEVStableGarch loaded and dgp.R sourced.

# One replication: simulate a series from the cell's model, fit it, and return a
# tidy data frame with one row per parameter plus the Phase 1 diagnostics.
mc_fit_once <- function(cell, n, seed, hessian = FALSE) {
  nms <- names(cell$true)
  fail <- function(msg, elapsed = NA_real_)
    data.frame(cell = cell$cell, table = cell$table, param_set = cell$param_set,
               param = nms, true = unname(cell$true), estimate = NA_real_, seed = seed,
               converged = FALSE, grad_rel = NA_real_, hess_rcond = NA_real_,
               iterations = NA_real_, n_at_bound = NA_integer_, persistence = NA_real_,
               elapsed = elapsed, error = msg, stringsAsFactors = FALSE)

  res <- tryCatch({
    t <- system.time({
      y <- suppressWarnings(gs_sim(cell$model, n = n, seed = seed)$y)
      f <- suppressWarnings(gs_fit(y, cell$spec, hessian = hessian))
    })
    list(f = f, elapsed = unname(t[3]))
  }, error = function(e) e)

  if (inherits(res, "error")) return(fail(conditionMessage(res)))
  f <- res$f; dg <- f$diagnostics
  ok <- (f$convergence == 0) && !isTRUE(dg$cap_hit) && isTRUE(dg$grad_rel < 1e-2)
  data.frame(cell = cell$cell, table = cell$table, param_set = cell$param_set,
             param = nms, true = unname(cell$true), estimate = unname(coef(f)[nms]),
             seed = seed, converged = ok, grad_rel = dg$grad_rel,
             hess_rcond = if (is.null(dg$hess_rcond)) NA_real_ else dg$hess_rcond,
             iterations = if (is.null(dg$iterations)) NA_real_ else dg$iterations,
             n_at_bound = length(f$at_bound), persistence = f$persistence,
             elapsed = res$elapsed, error = NA_character_, stringsAsFactors = FALSE)
}

# All R replications of one cell, in parallel.
mc_run_cell <- function(cell, n = 2500, R = 100, base_seed = 1000,
                        hessian = FALSE, cores = 5) {
  cores <- max(1, min(cores, parallel::detectCores()))
  reps <- parallel::mclapply(seq_len(R),
    function(r) mc_fit_once(cell, n, seed = base_seed + r, hessian = hessian),
    mc.cores = cores)
  bad <- vapply(reps, function(x) inherits(x, "try-error") || !is.data.frame(x), logical(1))
  if (any(bad)) reps[bad] <- lapply(which(bad), function(i)
    mc_fit_once(cell, n, seed = base_seed + i, hessian = hessian))  # serial retry
  do.call(rbind, reps)
}

# Every cell, with a one-line progress report per cell.
mc_run_all <- function(dgps, n = 2500, R = 100, base_seed = 1000,
                       hessian = FALSE, cores = 5, verbose = TRUE) {
  out <- vector("list", length(dgps))
  for (i in seq_along(dgps)) {
    if (verbose) cat(sprintf("[%2d/%2d] %-14s ", i, length(dgps), dgps[[i]]$cell))
    tt <- system.time(out[[i]] <- mc_run_cell(dgps[[i]], n, R, base_seed, hessian, cores))
    if (verbose) {
      df <- out[[i]]; cr <- mean(df$converged[df$param == df$param[1]])
      cat(sprintf("%5.0fs  conv %3.0f%%\n", tt[3], 100 * cr))
    }
  }
  do.call(rbind, out)
}
