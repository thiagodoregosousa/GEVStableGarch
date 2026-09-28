test_that("Enter invalid models gives error", {
  
  expect_error(GEVStableGarch_model_spec(p = 0))
  expect_error(GEVStableGarch_model_spec(q = -1))
  expect_error(GEVStableGarch_model_spec(m = -1))
  expect_error(GEVStableGarch_model_spec(cond.dist = "ged"))
  expect_error(GEVStableGarch_model_spec(aparch = FALSE, include_mean = "TRUE"))
  expect_error(GEVStableGarch_model_spec(aparch = "FALSE"))
  expect_equal(class(GEVStableGarch_model_spec()), "GEVStableGarch_spec")
  
})


test_that("Enter invalid parameters for fixed model spec", {
  
  spec_garch11_stable = GEVStableGarch_model_spec(q = 1)
  
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, ar = c(2,2)))
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, ma = c(2,2)))
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, alpha = c(1), beta = c(2)))
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, omega = 0.01, alpha = 3))
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, omega = 0.01, alpha = 3, beta = 0.2, aparch = 2))
  expect_error(GEVStableGarch_set_params(spec = spec_garch11_stable, omega = 0.01, alpha = 3, beta = 0.2, stable_alpha = 3))
  
})


