
library(ggplot2)
#-------------------------------------------------------------------------------
# Example 1: Sampling distribution of mean

# This example demonstrate some of the simpler uses of SimDesign,
# particularly for classroom settings. The only factor varied in this simulation
# is sample size.

# skeleton functions to be saved and edited
SimFunctions()

#### Step 1 --- Define your conditions under study and create design data.frame

Design <- createDesign(N = c(100, 300, 500), 
                       mean.gat = c(-3,1),
                       sd.gat = c(0.5, 3),
                       nu.gat = c(1,3),
                       d.gat = c(2,4),
                       xi.gat = c(0.5, 3))


#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 2 --- Define generate, analyse, and summarise functions

# help(Generate)
Generate <- function(condition, fixed_objects) {
  dat <- with(condition, rgat(n = N, mean = mean.gat, sd = sd.gat, nu = nu.gat, d = d.gat, xi = xi.gat)  ) 
  dat
}

# help(Analyse)
Analyse <- function(condition, dat, fixed_objects) {
  ret = gat.fit(dat, control = list(trace = FALSE))$pars
  names(ret) = c("mean.gat", "sd.gat", "nu.gat", "d.gat", "xi.gat")
  ret
}

# help(Summarise)
Summarise <- function(condition, results, fixed_objects) {
  # 'results' is a matrix of all the estimates from Analyse()
  # 'condition' contains the true population parameters defined in your Design object
  
  # Define the true population parameters (psi) based on your conditions
  # (Assuming your Design object columns match these names)
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



#~~~~~~~~~~~~~~~~~~~~~~~~
#### Step 3 --- Collect results by looping over the rows in design

# run with more standard number of replications
Final <- runSimulation(design=Design, replications=1000,
                       generate=Generate, analyse=Analyse, summarise=Summarise,
)
Final
print(Final, n  = 96)
plot(Final)

library(dplyr)
library(reshape2)
library(tidyr)
plot.var = "bias" # "bias" or "RMESE"
Final.plot <- Final %>% select(N, mean.gat, sd.gat, nu.gat, d.gat, xi.gat, contains(plot.var))
Final.longer <- Final.plot %>%
  pivot_longer(
    contains(plot.var),
    names_to = "type",
    values_to = "value"
  )

Final.longer <- Final.longer %>% mutate(
  param_combo = paste0(
    "mean=", mean.gat,
    ", sd=", sd.gat,
    ", nu=", nu.gat,
    ", d=", d.gat,
    ", xi=", xi.gat
  ))

pdf(paste0("../../vignettes/benchmarks/dist-gat/",plot.var,"_vs_N_gat.pdf"), width = 16, height = 16)
ggplot(data = Final.longer) +
  aes(x = N, y =  value, colour = type) +
  geom_line() +  geom_point(size = 3, stroke = 1.1) +
  xlab("Sample size") +
  ylab("RMSE") + 
  facet_wrap(~param_combo) + 
  ggtitle(plot.var) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.margin = margin(20, 20, 20, 20)
  )
dev.off()

