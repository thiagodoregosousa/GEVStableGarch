
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA



################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  Several                 This is a miscelaneous script that must be run and
#                          analysed after building the package in both 
#                          mac and window sistems
################################################################################


# Clean the Environment
rm(list=ls())



# ------------------------------------------------------------------------------
# gsFit 
# ------------------------------------------------------------------------------



# 1: Estimate garch(1,1) with all conditional distributions and the
# "sqp" algorithm
# Expect: Must achieve convergence for all models
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]
# x = 100*sp500dge[, 1]
mylist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", 
           "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
  print(mylist[i])
  model = gsFit(data = x, formula = ~garch(1,1), cond.dist = mylist[i]) 
}



# 2: Estimate garch(1,1) with all conditional distributions and the
# "nlminb" algorithm
# Expect: Must achieve convergence for all models
mylist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", 
           "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
  print(mylist[i])
  model = gsFit(data = x, formula = ~garch(1,1), cond.dist = mylist[i],
                algorithm = "nlminb") 
}



# 3: Estimate garch(1,1) with all conditional distributions and the
# "nlminb+nm" algorithm
# Expect: Must achieve convergence for all models
mylist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", 
           "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
  print(mylist[i])
  model = gsFit(data = x, formula = ~garch(1,1), cond.dist = mylist[i],
                algorithm = "nlminb+nm") 
}



# 4: Estimate garch(1,1) with all conditional distributions and the
# "sqp.restriction" algorithm
# Expect: Must achieve convergence for all models
mylist = c("stableS1", "gev", "gat", "norm", "std", 
           "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
  print(mylist[i])
  model = gsFit(data = x, formula = ~garch(1,1), cond.dist = mylist[i],
                algorithm = "sqp.restriction") 
}




# ------------------------------------------------------------------------------
# gsSelect
# ------------------------------------------------------------------------------



# 1: Select the best arma-aparch model for the dem2gbp dataset 
# with conditional stableS1 distribution and sqp algorithm
# Expect: Best Model is arma(0, 0) + aparch(1, 2)
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]
model.best = gsSelect(data = x, order.max = c(2,2,2,2), 
             is.aparch = TRUE, cond.dist = "stableS1")



# 2: Select the best arch model for the dem2gbp dataset 
# with conditional norm distribution
# Expect: Best Model is garch(10, 0)
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]
model.best = gsSelect(data = x, order.max = c(0,0,10,0),cond.dist = "norm")



# 3: Select the best arma-garch model for the dem2gbp dataset 
# with conditional gat distribution
# Expect: Best Model is garch(1, 1)
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]
model.best = gsSelect(data = x, order.max = c(1,1,1,1),cond.dist = "gat")



# ------------------------------------------------------------------------------
# gsMomentAparch
# ------------------------------------------------------------------------------



# Computation of the Moment E( |Z| - gamma Z) ^ delta for all the conditional 
# distributions

gsMomentAparch(cond.dist = "stableS1", shape = 1.1, skew = 0, delta = 1.01, gm = 0.99999)

gsMomentAparch(cond.dist = "gev", shape = -4, skew = 0, delta = 1.4, gm = 0)

gsMomentAparch(cond.dist = "gat", shape = c(1.9,2.3), skew = 0.5, delta = 0.4, gm = 0)

gsMomentAparch(cond.dist = "norm", shape = c(1.9,2.3), skew =1, delta = 11.4, gm = -0.999)

gsMomentAparch(cond.dist = "std", shape = 2.001, skew = -0.5, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "sstd", shape = 2.001, skew = 0.11, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "skstd", shape = 5.001, skew = 0.11, delta = 3, gm = -0.5)

gsMomentAparch(cond.dist = "ged", shape = 6, skew = 0.11, delta = 5.11, gm = -0.5)



# ------------------------------------------------------------------------------
# gsSelect and gsSim
# ------------------------------------------------------------------------------



# 1: Simulate ar(1)-garch(1,1)-gev and estimate the model
# Expect: estimated parameters next to the simulated model
gev.spec = gsSpec(model = list(ar = 0.4, alpha = 0.1, omega = 0.11, 
                               beta = 0.3, shape = 0.3, delta = 2),
                  cond.dist = "gev", rseed = 1001)

gev.sample = gsSim( spec = gev.spec, n = 2000)
model  <- gsFit(data = as.vector(gev.sample[,1]), formula = ~arma(1,0)+garch(1,1),
              cond.dist = "gev", algorithm = "sqp")


# 2: Simulate arma(1,2)-aparch(2,1)-stableS1 and estimate the model
# Expect: estimated parameters next to the simulated model
stable.spec = gsSpec(model = list(mu = -0.7, ar = c(-0.5), ma = c(0.7,0.3), 
		   omega = 0.4, alpha = 0.05,  
               beta = 0.7, shape = 1.3, delta = 1, skew = -0.7),
                  cond.dist = "stableS1", rseed = 1001)

stable.sample = gsSim( spec = stable.spec, n = 2000)
model  <- gsFit(data = as.vector(stable.sample[,1]), 
	formula = ~arma(1,2)+garch(1,1),
      cond.dist = "stableS1", algorithm = "sqp")



# 4: Simulate arma(1,1)-garch(1,1)-norm non-stationary and estimate the model
#with both "sqp" and "sqp.restriction algorithms" and conditional stableS1 
# distribution.
# Expect: a non-stationary model with the sqp algorithm and a stationary 
# model with the sqp.restriction algorithm.

norm.spec = gsSpec(model = list(mu = -0.7, ar = c(-0.5), ma = c(0.3), 
		   omega = 0.4, alpha = 0.4,  
               beta = 0.7, shape = 1.3, delta = 1, skew = -0.7),
                  cond.dist = "norm", rseed = 1001)
norm.sample = gsSim( spec = norm.spec, n = 2000)

model.sqp  <- gsFit(data = as.vector(norm.sample[,1]), 
	formula = ~arma(1,1)+garch(1,1),
      cond.dist = "stableS1", algorithm = "sqp")
stationarity.condition = .stationarityAparch(model = list(alpha = model.sqp@fit$par[5], beta = model.sqp@fit$par[6],
                                 skew = model.sqp@fit$par[7], shape = model.sqp@fit$par[8], gm = 0,
                                 delta = 1), formula = .getFormula(~arma(1,1)+garch(1,1)),
                    cond.dist = "stableS1")
stationarity.condition < 1

model.sqp  <- gsFit(data = as.vector(norm.sample[,1]), 
	formula = ~arma(1,1)+garch(1,1),
      cond.dist = "stableS1", algorithm = "sqp.restriction")
stationarity.condition = .stationarityAparch(model = list(alpha = model.sqp@fit$par[5], beta = model.sqp@fit$par[6],
                                 skew = model.sqp@fit$par[7], shape = model.sqp@fit$par[8], gm = 0,
                                 delta = 1), formula = .getFormula(~arma(1,1)+garch(1,1)),
                    cond.dist = "stableS1")
stationarity.condition < 1




