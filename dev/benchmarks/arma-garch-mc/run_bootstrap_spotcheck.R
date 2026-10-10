#!/usr/bin/env Rscript
# Bootstrap spot-check: for the stable AR(1)-GARCH(1,1) cell 9.7:phi1, compare
# 95% Wald coverage (Hessian SE) with parametric-bootstrap percentile coverage.
# Part B found that the mean-equation parameters mu and ar1 undercover badly with
# Wald under stable (infinite-variance) innovations; this shows the bootstrap
# restores calibration. Small by design (expensive: R x B refits).
#
# Usage (from this directory): Rscript run_bootstrap_spotcheck.R <R> <B> <n> <cores>

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
args  <- commandArgs(trailingOnly = TRUE)
R_    <- if (length(args) >= 1) as.integer(args[1]) else 80L
B_    <- if (length(args) >= 2) as.integer(args[2]) else 99L
n_    <- if (length(args) >= 3) as.integer(args[3]) else 2500L
cores <- if (length(args) >= 4) as.integer(args[4]) else 5L

if (requireNamespace("GEVStableGarch", quietly = TRUE)) library(GEVStableGarch) else
  devtools::load_all("../../..", quiet = TRUE)
source("dgp.R")

baseline <- mc_load_baseline("baseline_2012_thesis.csv")
cell <- Filter(function(d) d$cell == "9.7:phi1", mc_dgps(baseline))[[1]]
nms <- names(cell$true); truth <- cell$true
cat(sprintf("Bootstrap spot-check %s: R=%d, B=%d, n=%d, cores=%d\n\n", cell$cell, R_, B_, n_, cores))

one <- function(r) {
  y <- suppressWarnings(gs_sim(cell$model, n = n_, seed = 3000 + r)$y)
  fit <- tryCatch(suppressWarnings(gs_fit(y, cell$spec, hessian = TRUE)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  est <- coef(fit)[nms]; se <- sqrt(diag(vcov(fit)))[nms]
  wald_lo <- est - stats::qnorm(0.975) * se; wald_hi <- est + stats::qnorm(0.975) * se
  bt <- tryCatch(suppressWarnings(gs_bootstrap(fit, n_boot = B_, seed = 7000 + r)),
                 error = function(e) NULL)
  bo_lo <- bo_hi <- stats::setNames(rep(NA_real_, length(nms)), nms)
  if (!is.null(bt)) { bo_lo <- bt$ci[nms, "lower"]; bo_hi <- bt$ci[nms, "upper"] }
  data.frame(param = nms, true = unname(truth),
             wald = unname(truth >= wald_lo & truth <= wald_hi),
             boot = unname(truth >= bo_lo & truth <= bo_hi),
             stringsAsFactors = FALSE)
}

t0 <- Sys.time()
reps <- parallel::mclapply(seq_len(R_), one, mc.cores = max(1, min(cores, parallel::detectCores())))
res <- do.call(rbind, Filter(is.data.frame, reps))
cat(sprintf("Done in %.1f min (%d/%d reps usable).\n\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), length(unique(res)) , R_))

agg <- aggregate(cbind(wald, boot) ~ param, data = res, FUN = function(z) mean(z, na.rm = TRUE))
agg <- agg[match(nms, agg$param), ]
out_dir <- "../../mc/results"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(agg, file.path(out_dir, sprintf("bootstrap_spotcheck_phi1_n%d.csv", n_)), row.names = FALSE)

cat("=== 95% coverage: Wald (Hessian) vs parametric bootstrap (percentile) ===\n")
print(data.frame(param = agg$param, wald = round(agg$wald, 2), bootstrap = round(agg$boot, 2)),
      row.names = FALSE)
