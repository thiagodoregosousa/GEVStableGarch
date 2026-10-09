#!/usr/bin/env Rscript
# Consistency check: run a subset of cells across several sample sizes and show
# that per-parameter RMSE falls as n grows. Used first to validate that the GEV
# ARMA(2,2) soft spot at n = 2500 is small-sample variance (the estimator is
# median-unbiased there), and reused for Part B's consistency grid.
#
# Usage (from this directory):
#   Rscript run_consistency.R "<cell regex>" <R> <cores> "<n1,n2,...>"
# Defaults: cells "9.15|9.16", R = 100, cores = 5, n = 2500,5000,10000.

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
args    <- commandArgs(trailingOnly = TRUE)
pattern <- if (length(args) >= 1) args[1] else "9.15|9.16"
R_      <- if (length(args) >= 2) as.integer(args[2]) else 100L
cores   <- if (length(args) >= 3) as.integer(args[3]) else 5L
ns      <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else c(2500L, 5000L, 10000L)

if (requireNamespace("GEVStableGarch", quietly = TRUE)) library(GEVStableGarch) else
  devtools::load_all("../../..", quiet = TRUE)
source("dgp.R"); source("run_mc.R"); source("summarise.R")

baseline <- mc_load_baseline("baseline_2012_thesis.csv")
dgps     <- Filter(function(d) grepl(pattern, d$cell), mc_dgps(baseline))
cat(sprintf("Consistency: cells {%s}, n = {%s}, R = %d, cores = %d\n\n",
            paste(vapply(dgps, `[[`, "", "cell"), collapse = ", "),
            paste(ns, collapse = ", "), R_, cores))

out_dir <- "../../mc/results"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
summ_by_n <- list()
for (n in ns) {
  cat(sprintf("=== n = %d ===\n", n))
  res <- mc_run_all(dgps, n = n, R = R_, base_seed = 1000, hessian = FALSE, cores = cores)
  s <- mc_summarise(res); s$n <- n
  summ_by_n[[as.character(n)]] <- s
  saveRDS(res, file.path(out_dir, sprintf("consistency_%s_n%d_raw.rds", gsub("[^0-9]", "", pattern), n)))
  cat("\n")
}

summ <- do.call(rbind, summ_by_n)
write.csv(summ, file.path(out_dir, "consistency_summary_long.csv"), row.names = FALSE)

# Trend report: for each (cell, param), RMSE at each n and the 2012 n=2500 value.
base_key <- baseline; base_key$cell <- paste(base_key$table, base_key$param_set, sep = ":")
ref <- base_key[c("cell", "param", "rmse_2012")]
tab <- reshape(summ[c("cell", "param", "n", "rmse")],
               idvar = c("cell", "param"), timevar = "n", direction = "wide")
tab <- merge(tab, ref, by = c("cell", "param"), all.x = TRUE)
names(tab) <- sub("^rmse\\.", "rmse_n", names(tab))
tab <- tab[order(tab$cell, tab$param), ]
write.csv(tab, file.path(out_dir, "consistency_trend.csv"), row.names = FALSE)

cat("\n=== RMSE by sample size (consistency) ===\n")
print(tab, row.names = FALSE, digits = 3)
largest <- paste0("rmse_n", max(ns))
cat(sprintf("\nParams where RMSE at n=%d <= 2012 (n=2500): %d of %d\n",
            max(ns), sum(tab[[largest]] <= tab$rmse_2012, na.rm = TRUE), nrow(tab)))
cat(sprintf("Params where RMSE falls monotonically with n: %d of %d\n",
            sum(apply(tab[grep("^rmse_n", names(tab))], 1,
                      function(r) all(diff(as.numeric(r)) <= 1e-9)), na.rm = TRUE), nrow(tab)))
