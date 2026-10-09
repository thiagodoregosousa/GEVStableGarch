# Aggregate the raw replications to per-parameter mean/bias/RMSE and compare to
# the 2012 thesis baseline. See MC_PLAN.md for the pass criterion.
#
# Requires dgp.R sourced (for mc_load_baseline).

# Per (cell, parameter): mean estimate, bias, RMSE over the converged replications,
# plus the convergence rate.
mc_summarise <- function(results) {
  parts <- split(results, list(results$cell, results$param), drop = TRUE)
  rows <- lapply(parts, function(df) {
    conv <- df[df$converged & !is.na(df$estimate), ]
    est <- conv$estimate; true <- df$true[1]
    data.frame(cell = df$cell[1], table = df$table[1], param_set = df$param_set[1],
               param = df$param[1], true = true, n_conv = nrow(conv),
               conv_rate = nrow(conv) / nrow(df),
               mean_est = if (nrow(conv)) mean(est) else NA_real_,
               bias = if (nrow(conv)) mean(est) - true else NA_real_,
               rmse = if (nrow(conv)) sqrt(mean((est - true)^2)) else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$table, out$param_set, out$param), ]
}

# Per (cell, parameter, n): RMSE, bias, mean Hessian SE, empirical SD and the
# 95% Wald coverage (fraction of converged fits whose 95% interval covers the
# truth). Requires an `n` column and the `se` column from a hessian = TRUE run.
mc_coverage <- function(results) {
  parts <- split(results, list(results$cell, results$param, results$n), drop = TRUE)
  rows <- lapply(parts, function(df) {
    conv <- df[df$converged & !is.na(df$estimate), ]
    est <- conv$estimate; se <- conv$se; true <- df$true[1]; has <- is.finite(se)
    data.frame(cell = df$cell[1], table = df$table[1], param_set = df$param_set[1],
               param = df$param[1], n = df$n[1], true = true, n_conv = nrow(conv),
               rmse = if (nrow(conv)) sqrt(mean((est - true)^2)) else NA_real_,
               bias = if (nrow(conv)) mean(est) - true else NA_real_,
               mean_se = if (any(has)) mean(se[has]) else NA_real_,
               emp_sd = if (nrow(conv) > 1) stats::sd(est) else NA_real_,
               se_rate = mean(has),
               coverage95 = if (any(has))
                 mean(abs(est[has] - true) <= stats::qnorm(0.975) * se[has]) else NA_real_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[order(out$cell, out$param, out$n), ]
}

# Join the 2012 RMSE and flag non-degradation. Stable location parameters
# (mu, omega) are marked because they are only comparable if the thesis used the
# S0 parametrization (see MC_PLAN.md caveat 2); they do not count toward the gate.
mc_compare <- function(summary, baseline, tol = 0.10) {
  key <- c("table", "param_set", "param")
  m <- merge(summary, baseline[c(key, "rmse_2012", "est_2012")], by = key, all.x = TRUE)
  stable_tables <- c("9.7", "9.11", "9.12", "9.13", "9.14")
  m$ratio <- m$rmse / m$rmse_2012
  m$comparable <- !(m$param %in% c("mu", "omega") & m$table %in% stable_tables)
  m$pass <- !m$comparable | (m$ratio <= (1 + tol))
  m[order(m$table, m$param_set, match(m$param, unique(m$param))), ]
}

# One-line-per-cell and overall verdict.
mc_report <- function(comparison, tol = 0.10) {
  cmp <- comparison[comparison$comparable, ]
  cat(sprintf("Non-degradation gate (rmse_new <= rmse_2012 * %.2f), comparable params only:\n",
              1 + tol))
  by_cell <- split(cmp, cmp$cell)
  for (nm in names(by_cell)) {
    d <- by_cell[[nm]]
    worst <- d$param[which.max(d$ratio)]
    cat(sprintf("  %-14s  pass %2d/%2d   median ratio %.2f   worst %s=%.2f\n",
                nm, sum(d$pass), nrow(d), stats::median(d$ratio, na.rm = TRUE),
                worst, max(d$ratio, na.rm = TRUE)))
  }
  cat(sprintf("\nOVERALL: %d/%d comparable parameters pass (%.0f%%).\n",
              sum(cmp$pass), nrow(cmp), 100 * mean(cmp$pass)))
  flagged <- comparison[!comparison$comparable, ]
  if (nrow(flagged))
    cat(sprintf("(%d stable mu/omega rows flagged not-comparable, excluded from the gate.)\n",
                nrow(flagged)))
  invisible(comparison)
}
