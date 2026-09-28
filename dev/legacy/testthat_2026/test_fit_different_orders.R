test_that("GARCH(1,0) with GAT simulation and estimation", {
  
  skip_on_cran()
  
  set.seed(122)
  
  true <- c(mu = 0, omega = 0.1, alpha1 = 0.1,
            skew = 0, shape1 = 1, shape2 = 2)
  
  
  spec.gat = GEVStableGarch::gsSpec(model = list(mu = true["mu"], omega = true["omega"],
                                                 alpha = true["alpha1"], 
                                                 shape = c(true["shape1"], true["shape2"])), 
                                    cond.dist = "gat")
  sim.gat = as.numeric(GEVStableGarch::gsSim(spec = spec.gat, n = 10000)[,1])
  fit <- GEVStableGarch::gsFit(
    data = sim.gat, formula = ~garch(1,0),
    cond.dist = "gat", include.mean = TRUE,
    algorithm = "sqp")@fit$par
  cbind(true,fit)
  expect_equal(fit["mu"],     true["mu"],     tolerance = 0.01)
  expect_equal(fit["omega"],  true["omega"],  tolerance = 0.01)
  expect_equal(fit["alpha1"], true["alpha1"], tolerance = 0.01)
  expect_equal(fit["shape1"], true["shape1"], tolerance = 0.01)
  
})
