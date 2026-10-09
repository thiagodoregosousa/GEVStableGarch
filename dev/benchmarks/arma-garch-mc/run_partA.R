#!/usr/bin/env Rscript
# Part A of the Monte Carlo study: non-degradation against the 2012 thesis.
# Replicates the thesis MLE tables (9.2, 9.7, 9.11-9.17) at n = 2500, R = 100,
# and checks that the new RMSE does not exceed the 2012 RMSE by more than 10%.
#
# Usage (from this directory):
#   Rscript run_partA.R [R] [cores] [n]
# Defaults: R = 100, cores = 5, n = 2500. Raw replications and the comparison
# table are written to ../../mc/results/ (results are not kept in the benchmark
# folder; the tables go into the pkgdown vignette later).

# Keep each forked worker single-threaded so `cores` parallel fits do not
# oversubscribe the machine (libstable4u/BLAS can otherwise spawn extra threads).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

args  <- commandArgs(trailingOnly = TRUE)
R_    <- if (length(args) >= 1) as.integer(args[1]) else 100L
cores <- if (length(args) >= 2) as.integer(args[2]) else 5L
n_    <- if (length(args) >= 3) as.integer(args[3]) else 2500L
only  <- if (length(args) >= 4) args[4] else "all"   # "all", "gev" or "stable"
restarts <- if (length(args) >= 5) as.integer(args[5]) else 0L

if (requireNamespace("GEVStableGarch", quietly = TRUE)) {
  library(GEVStableGarch)
} else {
  message("GEVStableGarch not installed; loading the source tree with devtools.")
  devtools::load_all("../../..", quiet = TRUE)
}

source("dgp.R"); source("run_mc.R"); source("summarise.R")

baseline <- mc_load_baseline("baseline_2012_thesis.csv")
dgps     <- mc_dgps(baseline)
if (only != "all") dgps <- Filter(function(d) d$dist == only, dgps)

cat(sprintf("Part A (%s): %d cells, R = %d, n = %d, cores = %d\n\n",
            only, length(dgps), R_, n_, cores))
t0 <- Sys.time()
results <- mc_run_all(dgps, n = n_, R = R_, base_seed = 1000, hessian = FALSE,
                      cores = cores, restarts = restarts)
cat(sprintf("\nDone in %.1f min.\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

out_dir <- "../../mc/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("partA_R%d_n%d%s%s", R_, n_, if (only == "all") "" else paste0("_", only),
               if (restarts > 0) paste0("_rs", restarts) else "")
saveRDS(results, file.path(out_dir, paste0(tag, "_raw.rds")))

summary_tab <- mc_summarise(results)
comparison  <- mc_compare(summary_tab, baseline)
write.csv(summary_tab, file.path(out_dir, paste0(tag, "_summary.csv")), row.names = FALSE)
write.csv(comparison,  file.path(out_dir, paste0(tag, "_comparison.csv")), row.names = FALSE)

cat("Session: ", R.version.string, " | GEVStableGarch ",
    as.character(utils::packageVersion("GEVStableGarch")), "\n\n", sep = "")
mc_report(comparison)
