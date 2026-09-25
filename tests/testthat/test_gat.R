test_that("Passing invalid input parameters for dgat", {
  # sd cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1))
  
  # d cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0.1, d = -1, xi = 1))
  
  # xi cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 1, d = 0.1, xi = 0))
  
  # nu cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0, d = 0.1, xi = 1))
})


test_that("Passing invalid input parameters for dgat", {
  # sd cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = -1, nu = 2, d = 3, xi = 1))
  
  # d cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0.1, d = -1, xi = 1))
  
  # xi cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 1, d = 0.1, xi = 0))
  
  # nu cannot be <= 0
  expect_error(dgat(rnorm(1000), mean = 0, sd = 3, nu = 0, d = 0.1, xi = 1))
})


test_that("Integration of dgat is approximately 1", {
  expect_equal(integrate(dgat , lower = -Inf, upper = Inf, 
                         mean = 0.3, sd = 1, nu = 1,
                         d = 2, xi = 1)$value, 1, tolerance = 1e-5)
})


test_that("Distribution function behaviour close to -Inf and Inf", {
  expect_equal(pgat(-1000, mean = 0.3, sd = 1, nu = 1,
                    d = 2, xi = 1), 0, tolerance = 1e-5)
  
  expect_equal(pgat(1000, mean = 0.3, sd = 1, nu = 1,
                    d = 2, xi = 1), 1, tolerance = 1e-5)
})


test_that("Distribution and Density property: integral( pgat, a, b) = pgat(b) - pgat(a)", {
  a = -3
  b = 5
  prob_with_integral = integrate(dgat , lower = a, upper = b, 
                                 mean = 0.3, sd = 2, nu = 1,
                                 d = 2, xi = 1)$value
  pgat_at_b = pgat(b, mean = 0.3, sd = 2, nu = 1,
                   d = 2, xi = 1)
  pgat_at_a = pgat(a, mean = 0.3, sd = 2, nu = 1,
                   d = 2, xi = 1)
  
  expect_equal(prob_with_integral, pgat_at_b - pgat_at_a, tolerance = 1e-5)
})


test_that("Distribution and Density property: integral( pgat, a, b) = pgat(b) - pgat(a)", {
  
  prob_value = 0.344
  quantile_value = qgat(prob_value, mean = 0.3, sd = 2, nu = 1,
       d = 2, xi = 1)
  prob_value_with_distribution = pgat(quantile_value, mean = 0.3, sd = 2, nu = 1,
       d = 2, xi = 1)
  expect_equal(prob_value, prob_value_with_distribution, tolerance = 1e-5)

})



test_that("Using random function to compute expectation and comparing with integral", {
  set.seed(123)
  
  a <- -Inf; b <- Inf
  mean <- 0.3; sd <- 2; nu <- 1; d <- 2; xi <- 1
  
  integral_mean <- integrate(function(x) x * dgat(x, mean = mean, sd = sd, nu = nu, d = d, xi = xi),
                             lower = a, upper = b)$value
  
  n <- 1e5
  mc_mean <- mean(rgat(n, mean = mean, sd = sd, nu = nu, d = d, xi = xi))
  
  expect_equal(mc_mean, integral_mean, tolerance = 0.05)
})


test_that("Function to detect invalid parameters", {
  
  expect_true(gat.valid.pars(0, 1, 1, 1, 1))
  expect_true(gat.valid.pars(-100, 100, 1000, 10000, 10))
  expect_false(gat.valid.pars(0, -1, 1, 1, 1))
  expect_false(gat.valid.pars(0, 1, 1, -1, 1))
  expect_false(gat.valid.pars(0, 1, 1, 1, -0.001))
})



test_that("fitting parameters with MLE gives reasonable bias and RMSE for big sample size", {

  Design <- SimDesign::createDesign(N = c(1000,5000), 
                         mean.gat = 2,
                         sd.gat = 2,
                         nu.gat = 2,
                         d.gat = 2,
                         xi.gat = 2)
  Generate <- function(condition, fixed_objects) {
    dat <- with(condition, rgat(n = N, mean = mean.gat, sd = sd.gat, nu = nu.gat, d = d.gat, xi = xi.gat)  ) 
    dat
  }
  Analyse <- function(condition, dat, fixed_objects) {
    ret = gat.fit(dat, control = list(trace = FALSE))$pars
    names(ret) = c("mean.gat", "sd.gat", "nu.gat", "d.gat", "xi.gat")
    ret
  }
  Summarise <- function(condition, results, fixed_objects) {
    # assuming your Design object columns match these names)
    true_mean_gat <- condition$mean.gat
    true_sd_gat <- condition$sd.gat
    true_nu_gat <- condition$nu.gat
    true_d_gat <- condition$d.gat
    true_xi_gat <- condition$xi.gat
    
    # Calculate bias for each parameter using the built-in SimDesign functions
    bias_mean <- bias(results[, "mean.gat"], parameter = true_mean_gat)
    bias_sd <- bias(results[, "sd.gat"], parameter = true_sd_gat)
    bias_nu <- bias(results[, "nu.gat"], parameter = true_nu_gat)
    bias_d <- bias(results[, "d.gat"], parameter = true_d_gat)
    bias_xi <- bias(results[, "xi.gat"], parameter = true_xi_gat)
    
    # Calculate RMSE for each parameter using the built-in SimDesign functions
    RMSE_mean <- RMSE(results[, "mean.gat"], parameter = true_mean_gat)
    RMSE_sd <- RMSE(results[, "sd.gat"], parameter = true_sd_gat)
    RMSE_nu <- RMSE(results[, "nu.gat"], parameter = true_nu_gat)
    RMSE_d <- RMSE(results[, "d.gat"], parameter = true_d_gat)
    RMSE_xi <- RMSE(results[, "xi.gat"], parameter = true_xi_gat)
    
    # Return a named vector of the summary statistics
    ret <- c(
      Bias_mean_gat = bias_mean,
      Bias_sd_gat = bias_sd,
      Bias_nu_gat = bias_nu,
      Bias_d_gat = bias_d,
      Bias_xi_gat = bias_xi,
      RMSE_mean_gat = RMSE_mean,
      RMSE_sd_gat = RMSE_sd,
      RMSE_nu_gat = RMSE_nu,
      RMSE_d_gat = RMSE_d,
      RMSE_xi_gat = RMSE_xi
    )
    
    return(ret)
  }
  
  Final <- SimDesign::runSimulation(design=Design, replications=1000,
                         generate=Generate, analyse=Analyse, summarise=Summarise,
                         seed = 1:nrow(Design), progress = FALSE, verbose = FALSE)
  
  # expect RMSE to decay
  Final.RMSE.values = Final %>% dplyr::select(contains("RMSE")) %>% as.matrix()
  expect_true(sum ((Final.RMSE.values[2,] - Final.RMSE.values[1,]) < 0) == 5)
  
  # expect small RMSE, except for nu which is around 0.4 for N = 5000 
  expect_equal(Final$RMSE_mean_gat[2], 0, 0.2)
  expect_equal(Final$RMSE_sd_gat[2], 0, 0.2)
  expect_equal(Final$RMSE_nu_gat[2], 0, 0.4) # currently 'nu' is more difficult to estimate
  expect_equal(Final$RMSE_d_gat[2], 0, 0.2)
  expect_equal(Final$RMSE_xi_gat[2], 0, 0.2)
  
})


