#!/usr/bin/env Rscript
# Part B: consistency grid + Wald standard-error calibration on a representative
# subset of the DGPs. One pass with hessian = TRUE gives both the RMSE/bias
# trend across sample size and the 95% Wald coverage. Parametric-bootstrap
# coverage is left to a small separate spot-check (gs_bootstrap). See MC_PLAN.md.
#
# Usage (from this directory):
#   Rscript run_partB.R <R> <cores> "<n1,n2,...>"
# Defaults: R = 300, cores = 5, n = 500,1000,2500,5000.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
args  <- commandArgs(trailingOnly = TRUE)
R_    <- if (length(args) >= 1) as.integer(args[1]) else 300L
cores <- if (length(args) >= 2) as.integer(args[2]) else 5L
ns    <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(500L, 1000L, 2500L, 5000L)

# Representative subset: a well-identified and a harder cell per family, plus one
# APARCH cell each (stable and GEV).
subset_cells <- c("9.2:theta2", "9.2:theta3", "9.7:phi1", "9.7:phi3",
                  "9.11:set1", "9.15:set1")

if (requireNamespace("GEVStableGarch", quietly = TRUE)) library(GEVStableGarch) else
  devtools::load_all("../../..", quiet = TRUE)
source("dgp.R"); source("run_mc.R"); source("summarise.R")

baseline <- mc_load_baseline("baseline_2012_thesis.csv")
dgps <- Filter(function(d) d$cell %in% subset_cells, mc_dgps(baseline))
cat(sprintf("Part B: cells {%s}, n = {%s}, R = %d, cores = %d, hessian = TRUE\n\n",
            paste(vapply(dgps, `[[`, "", "cell"), collapse = ", "),
            paste(ns, collapse = ", "), R_, cores))

out_dir <- "../../mc/results"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
all_res <- list()
t0 <- Sys.time()
for (n in ns) {
  cat(sprintf("=== n = %d ===\n", n))
  res <- mc_run_all(dgps, n = n, R = R_, base_seed = 2000, hessian = TRUE, cores = cores)
  res$n <- n
  all_res[[as.character(n)]] <- res
  saveRDS(res, file.path(out_dir, sprintf("partB_n%d_raw.rds", n)))
  cat("\n")
}
cat(sprintf("Done in %.1f min.\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

results <- do.call(rbind, all_res)
cov <- mc_coverage(results)
write.csv(cov, file.path(out_dir, "partB_coverage.csv"), row.names = FALSE)

# Consistency: RMSE by n per parameter.
rmse_wide <- reshape(cov[c("cell", "param", "n", "rmse")],
                     idvar = c("cell", "param"), timevar = "n", direction = "wide")
cat("=== Consistency: RMSE by sample size ===\n")
print(rmse_wide, row.names = FALSE, digits = 3)

cat("\n=== 95% Wald coverage by sample size (target 0.95) ===\n")
cov_wide <- reshape(cov[c("cell", "param", "n", "coverage95")],
                    idvar = c("cell", "param"), timevar = "n", direction = "wide")
print(cov_wide, row.names = FALSE, digits = 2)

big <- cov[cov$n == max(ns), ]
cat(sprintf("\nAt n=%d: median coverage %.2f; params within [0.90,0.98]: %d/%d; SE available %.0f%%.\n",
            max(ns), stats::median(big$coverage95, na.rm = TRUE),
            sum(big$coverage95 >= 0.90 & big$coverage95 <= 0.98, na.rm = TRUE),
            sum(!is.na(big$coverage95)), 100 * mean(big$se_rate)))
