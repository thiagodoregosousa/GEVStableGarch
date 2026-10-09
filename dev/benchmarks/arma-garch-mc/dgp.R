# Data generating processes for the Monte Carlo study, reconstructed from the
# 2012 thesis baseline (baseline_2012_thesis.csv). Each cell is one (table,
# param_set): a gs_spec, a gs_model built from the true values, and the true
# parameter vector. See MC_PLAN.md.
#
# Reproducible: source this file after library(GEVStableGarch).

# Model code -> gs_spec orders. GARCH mode has no gamma/delta; APARCH estimates
# both. Stable GARCH fixes delta = 1, GEV GARCH fixes delta = 2 (package default).
.mc_spec <- function(model_code, dist_name) {
  dist <- switch(dist_name, gev = gs_gev(), stable = gs_stable(),
                 stop("unknown dist ", dist_name))
  cfg <- switch(model_code,
    ar1_garch11     = list(arma = c(1, 0), garch = c(1, 1), aparch = FALSE),
    arma11_garch11  = list(arma = c(1, 1), garch = c(1, 1), aparch = FALSE),
    arma22_garch22  = list(arma = c(2, 2), garch = c(2, 2), aparch = FALSE),
    arma22_aparch22 = list(arma = c(2, 2), garch = c(2, 2), aparch = TRUE),
    stop("unknown model_code ", model_code))
  gs_spec(arma = cfg$arma, garch = cfg$garch, aparch = cfg$aparch, dist = dist)
}

# Read the baseline and return the list of DGP cells.
mc_load_baseline <- function(path = "baseline_2012_thesis.csv") {
  b <- utils::read.csv(path, comment.char = "#", stringsAsFactors = FALSE,
                       colClasses = c(table = "character"))
  b$cell <- paste(b$table, b$param_set, sep = ":")
  b
}

mc_dgps <- function(baseline) {
  cells <- unique(baseline[c("cell", "table", "model", "dist", "param_set")])
  lapply(seq_len(nrow(cells)), function(i) {
    row <- cells[i, ]
    sub <- baseline[baseline$cell == row$cell, ]
    spec <- .mc_spec(row$model, row$dist)
    true <- stats::setNames(sub$true, sub$param)
    # gs_model validates the names against the spec, reorders to canonical order,
    # and checks admissibility; take the validated vector back from the model.
    model <- gs_model(spec, par = true)
    list(cell = row$cell, table = row$table, model_code = row$model,
         dist = row$dist, param_set = row$param_set, spec = spec, model = model,
         true = model$par, persistence = gs_stationarity(model))
  })
}

# Quick sanity print: every cell admissible and (for these stationary DGPs)
# with persistence below 1.
mc_check_dgps <- function(dgps) {
  data.frame(
    cell        = vapply(dgps, `[[`, "", "cell"),
    model       = vapply(dgps, `[[`, "", "model_code"),
    dist        = vapply(dgps, `[[`, "", "dist"),
    n_par       = vapply(dgps, function(d) length(d$true), 0L),
    persistence = round(vapply(dgps, `[[`, 0, "persistence"), 4),
    stringsAsFactors = FALSE)
}
