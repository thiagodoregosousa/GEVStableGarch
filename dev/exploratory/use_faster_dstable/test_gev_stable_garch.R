###############################################################################
#  test_gev_stable_garch.R
#
#  Demonstrates the full pipeline with stableS0 conditional distribution:
#    1. Source the self-contained gsStableGarch.R
#    2. Fit an ARMA(1,1)-GARCH(1,1) model to real data (dem2gbp)
#    3. Simulate a new series using the estimated parameters
#    4. Re-estimate the model on the simulated data
#    5. Compare the true (step 2) and re-estimated (step 4) parameters
###############################################################################

# Step 1: Load everything (auto-installs missing packages)
source("gsStableGarch.R")

# Step 2: Fit model to real data
data(dem2gbp, package = "fGarch")
x = dem2gbp[, 1]
cat("\n=================== Step 2: Fit to real data ===================\n")
fit_real = gsFit(data = x, formula = ~arma(1,1)+garch(1,1), cond.dist = "stableS0")

# Step 3: Simulate from the fitted parameters
cat("\n=============== Step 3: Simulate from fitted params ===============\n")
est = fit_real@fit$par
spec = gsSpec(
  model = list(
    mu     = est["mu"],
    ar     = est["ar1"],
    ma     = est["ma1"],
    omega  = est["omega"],
    alpha  = est["alpha1"],
    beta   = est["beta1"],
    delta  = 1,
    skew   = est["skew"],
    shape  = est["shape1"]
  ),
  cond.dist = "stableS0",
  rseed = 9731
)
n_sim = 5000
sim_data = gsSim(spec, n = n_sim, n.start = 500)
x_sim = as.numeric(sim_data[, "Series"])
cat("Simulated", n_sim, "observations from the fitted model\n")

# Step 4: Re-estimate on the simulated data
cat("\n============ Step 4: Re-estimate on simulated data ============\n")
fit_sim = gsFit(data = x_sim, formula = ~arma(1,1)+garch(1,1), cond.dist = "stableS0")

# Step 5: Compare parameters
cat("\n============ Step 5: Parameter comparison ============\n")
true_params = est
est_params = fit_sim@fit$par

comparison = data.frame(
  Parameter = names(true_params),
  True      = round(as.numeric(true_params), 6),
  Estimated = round(as.numeric(est_params), 6),
  Abs_Diff  = round(abs(as.numeric(true_params) - as.numeric(est_params)), 6)
)
print(comparison, row.names = FALSE)
