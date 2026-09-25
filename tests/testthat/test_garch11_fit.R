test_that("GARCH(1,1) with stableS0 simulation and estimation", {

    skip_on_cran()
    skip_if_not_installed("fGarch")
    skip_if_not_installed("stabledist")
    
    set.seed(123)
    
    true <- c(mu = 0, omega = 0.1, alpha1 = 0.1, beta1 = 0.8,
              skew = 0.2, shape1 = 1.95)
    
    
    spec.stable = GEVStableGarch::gsSpec(model = list(mu = true["mu"], omega = true["omega"],
                                      alpha = true["alpha1"], beta = true["beta1"], 
                                      skew = true["skew"], shape = true["shape1"]), 
                         cond.dist = "stableS0")
    sim.stable = as.numeric(GEVStableGarch::gsSim(spec = spec.stable, n = 10000)[,1])
    fit <- GEVStableGarch::gsFit(
      data = sim.stable, formula = ~garch(1,1),
      cond.dist = "stableS0", include.mean = TRUE,
      algorithm = "sqp")@fit$par
    
    expect_equal(fit["mu"],     true["mu"],     tolerance = 0.02)
    expect_equal(fit["omega"],  true["omega"],  tolerance = 0.01)
    expect_equal(fit["alpha1"], true["alpha1"], tolerance = 0.01)
    expect_equal(fit["beta1"],  true["beta1"],  tolerance = 0.01)
    expect_equal(fit["skew"],   true["skew"],   tolerance = 0.1)
    expect_equal(fit["shape1"], true["shape1"], tolerance = 0.05)
    
})


test_that("GARCH(1,1) with GEV simulation and estimation", {
  
  skip_on_cran()
  
  set.seed(123)
  
  true <- c(mu = 0, omega = 0.1, alpha1 = 0.1, beta1 = 0.4,
            shape1 = 0.15)
  
  
  spec.gev = GEVStableGarch::gsSpec(model = list(mu = true["mu"], omega = true["omega"],
                                                    alpha = true["alpha1"], beta = true["beta1"], 
                                                    shape = true["shape1"]), 
                                       cond.dist = "gev")
  sim.gev = as.numeric(GEVStableGarch::gsSim(spec = spec.gev, n = 10000)[,1])
  fit <- GEVStableGarch::gsFit(
    data = sim.gev, formula = ~garch(1,1),
    cond.dist = "gev", include.mean = TRUE,
    algorithm = "sqp")@fit$par
  expect_equal(fit["mu"],     true["mu"],     tolerance = 0.01)
  expect_equal(fit["omega"],  true["omega"],  tolerance = 0.01)
  expect_equal(fit["alpha1"], true["alpha1"], tolerance = 0.01)
  expect_equal(fit["beta1"],  true["beta1"],  tolerance = 0.01)
  expect_equal(fit["shape1"], true["shape1"], tolerance = 0.01)
  
})


test_that("GARCH(1,1) with GAT simulation and estimation", {
  
  skip_on_cran()
  
  set.seed(123)
  
  true <- c(mu = 0, omega = 0.4, alpha1 = 0.2, beta1 = 0.2,
            skew = 0.5, shape1 = 2.3, shape2 = 2)
  
  
  spec.gat = GEVStableGarch::gsSpec(model = list(mu = true["mu"], omega = true["omega"],
                                                 alpha = true["alpha1"], beta = true["beta1"], 
                                                 skew = true["skew"],
                                                 shape = c(true["shape1"], true["shape2"])), 
                                    cond.dist = "gat")
  sim.gat = as.numeric(GEVStableGarch::gsSim(spec = spec.gat, n = 10000)[,1])
  fit <- GEVStableGarch::gsFit(
    data = sim.gat, formula = ~garch(1,1),
    cond.dist = "gat", include.mean = TRUE,
    algorithm = "nlminb")@fit$par
  cbind(fit,true)
  expect_equal(fit["mu"],     true["mu"],     tolerance = 0.03)
  expect_equal(fit["omega"],  true["omega"],  tolerance = 0.03)
  expect_equal(fit["alpha1"], true["alpha1"], tolerance = 0.01)
  expect_equal(fit["beta1"],  true["beta1"],  tolerance = 0.03)
  expect_equal(fit["skew"], true["skew"], tolerance = 0.01)
  expect_equal(fit["shape1"], true["shape1"], tolerance = 0.03)
  expect_equal(fit["shape2"], true["shape2"], tolerance = 0.02)
  
})
  


test_that("GARCH(1,1) normal intercept matches fGarch estimation values", {
  
  data("dem2gbp", package = "fGarch")
  x <- dem2gbp[, 1]
  
  fgarch_garch11_norm <- fGarch::garchFit(data = x, formula = ~garch(1,1),
                                          cond.dist = "norm", include.mean = TRUE,
                                          algorithm = "nlminb+nm", trace = FALSE)@fit$params$params
  
  gevstablegarch_garch11_norm <- GEVStableGarch::gsFit(data = x, formula = ~garch(1,1),
                                       cond.dist = "norm", include.mean = TRUE,
                                       algorithm = "sqp", control = list(trace = FALSE))@fit$par
  
  expect_equal(fgarch_garch11_norm["mu"],     gevstablegarch_garch11_norm["mu"],     tolerance = 0.01)
  expect_equal(fgarch_garch11_norm["omega"],  gevstablegarch_garch11_norm["omega"],  tolerance = 0.01)
  expect_equal(fgarch_garch11_norm["alpha1"], gevstablegarch_garch11_norm["alpha1"], tolerance = 0.01)
  expect_equal(fgarch_garch11_norm["beta1"],  gevstablegarch_garch11_norm["beta1"],  tolerance = 0.01)
  
})
