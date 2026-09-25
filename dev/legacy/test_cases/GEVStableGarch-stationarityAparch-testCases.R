
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
# TEST CASES FOF FUNCTIONS: 
#                
#  .stationarityAparch
#  gsMomentAparch
#  .normMomentAparch 
#  .stdMomentAparch
#  .skstdMomentAparch              
#  .gedMomentAparch                
#  .gatMomentAparch 
#  .gevMomentAparch
#  .stableS1MomentAparch            
#  .stableS1SymmetricMomentGarch    
#  .stableS1SymmetricMomentAparch   
#  .stableS1MomentPowerGarch
################################################################################



# ------------------------------------------------------------------------------
# Test Cases for functions .normMomentAparch 
# ------------------------------------------------------------------------------

# E(z^2) = Var(z) = 1, where z ~ N(0,1) - OK
.normMomentAparch (delta = 2, gm = 0)

# E(z^0) = 1, where z ~ N(0,1) - OK
.normMomentAparch (delta = 0.0000001, gm = 0)

# Comparison with simulated data - OK
n = 100
nSimulations = 1000
deltaValues = runif(n, 0, 10)
gmValues = runif(n, -1, 1)
error = rep(NA,n)
for(i in 1:n)
{
  z = rnorm(nSimulations)
  error[i] = abs((mean((abs(z)-gmValues[i]*z)^deltaValues[i]) - 
                    .normMomentAparch (delta = deltaValues[i], gm = gmValues[i]))/
                   .normMomentAparch (delta = deltaValues[i], gm = gmValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for functions .gevMomentAparch
# ------------------------------------------------------------------------------

# E(z^0) = 1, where z ~ GEV(0,1) - OK
.gevMomentAparch(delta = 0.0000001, gm = 0, shape = 5)




# ------------------------------------------------------------------------------
# Test Cases for functions .stdMomentAparch (standard t-Student)
# ------------------------------------------------------------------------------


# Tests with the numerical integration computation - OK
n = 1000
shapeValues = runif(n, 2+0.01, 5)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,0,shapeValues[i]-0.09)
gmValues = runif(n,-1,1)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz(fun = "dstd",gm = gmValues[i], 
                                                    delta = deltaValues[i], nu = shapeValues[i])[1])
  functionValues[i] = .stdMomentAparch(shape = shapeValues[i], delta = deltaValues[i],
                                        gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for functions .gedMomentAparch (standard GED distribution)
# ------------------------------------------------------------------------------

# Tests with the numerical integration computation - OK
n = 100
shapeValues = runif(n, 0.1, 5)
deltaValues = runif(n,0,10)
gmValues = runif(n,-1,1)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz(fun = "dged",gm = gmValues[i], 
                                                    delta = deltaValues[i], nu = shapeValues[i])[1])
  functionValues[i] = .gedMomentAparch(shape = shapeValues[i], delta = deltaValues[i],
                                        gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)


# ------------------------------------------------------------------------------
# Test Cases for functions .skstdMomentAparch (skewed t-Student from Fernandez and Steel)
# ------------------------------------------------------------------------------
# We get numerical problems when delta/skew approaches 1.


# Tests with simulated data
# works well with errors around 0.01% of the true value.
# Sometimes the function .trueAparchMomentsWurtz fail because of the integration
# routine.


n = 1e4
shapeValues = runif(n, 2, 5)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
skewValues = rep(NA,n)
for(i in 1:n)
{
  deltaValues[i] = runif(n=1,0.05,shapeValues[i]-0.2) 
  if(i %% 2 == 0)
    skewValues[i] = runif(n=1,0.015,deltaValues[i]-0.02)  
  else
    skewValues[i] = runif(n=1,deltaValues[i]+0.02,max(deltaValues[i]+2,2))
}
trueValues = rep(NA,n)
functionValues = rep(NA,n)
error = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz(fun = "dskstd", gm = gmValues[i], delta = deltaValues[i], 
                                                    nu = shapeValues[i], xi = skewValues[i])[1])
  functionValues[i] = .skstdMomentAparch(shape = shapeValues[i], skew = skewValues[i], 
                                         delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(gmValues,skewValues,deltaValues,shapeValues,trueValues,functionValues)
summary(error)

# Compare the calculation of the moment with the G@RCH software
# For the dem2gbp-APARCH(1,1) the G@RCH reports: persistency = 0.9901654 but the true
# persistency is 0.9846524 as calculated with numerical integration on 
# the dsstd distribution

.skstdMomentAparch(shape = 4.211667, skew = exp(-0.100622), 
                   delta = 1.183602, gm = 0.149934)*0.135215 + 0.885756


as.numeric(.trueAparchMomentsWurtz(fun = "dsstd", gm = 0.149934, delta = 1.183602, 
                       nu = 4.211667, xi = exp(-0.100622))[1])*0.135215 + 0.885756

# For the dem2gbp-ARMA(1,1)-APARCH the G@RCH reports: persistency = 0.989947
# but the true persistency is 0.9845671.
.skstdMomentAparch(shape = 4.221416, skew = exp(-0.095899), 
                   delta = 1.202501, gm = 0.151121)*0.136845 + 0.884074

as.numeric(.trueAparchMomentsWurtz(fun = "dsstd", gm = 0.151121, delta = 1.202501, 
                 nu = 4.221416, xi = exp(-0.095899))[1])*0.136845 + 0.884074

# Conclusion: The G@RCH software is using the expression to calculate the 
# moments of skewed t-Student definied by Lambert and Laurent (2000, 2001),
# but the density defined in the garch model was reparameterized to be a 
# zero mean and unit variance.
# Another interesting conclusion is that if we work with the density 
# defined by Lambert and Laurent (without reparameterizing it to have 
# a zero mean and unit variance) we can also get the same garch estimated 
# parameters.


# ------------------------------------------------------------------------------
# Test Cases for functions .gatMomentAparch (standard GAt distribution)
# ------------------------------------------------------------------------------

# Tests with the numerical integration computation
# 20150413 - Test OK with error around 0.05% of the true value. 
# Sometimes the integrate function fails because of the integration interval.
# Note that the GAt density becames very picky when d approaches zero and thus,
# the integrate function will fail for these values. 
n = 2000
shapeValues = cbind(runif(n, 0.1, 5),runif(n, 0.5, 5))
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,0,shapeValues[i,1]*shapeValues[i,2]-0.09)
gmValues = runif(n,-1,1)
skewValues = runif(n,0,3)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz(fun = "dgat",gm = gmValues[i], 
                  delta = deltaValues[i], nu = shapeValues[i,1], d = shapeValues[i,2], 
                  xi = skewValues[i], lower = -Inf, upper = Inf)[1])
  functionValues[i] = .gatMomentAparch(shape = shapeValues[i,], delta = deltaValues[i],
                                        gm = gmValues[i], skew = skewValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)[which(error > 1),]
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for function .stableS1SymmetricMomentGarch
# ------------------------------------------------------------------------------

# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002) - OK
.stableS1SymmetricMomentGarch(shape = 2)
sqrt(2)*.normMomentAparch (delta = 1, gm = 0)

# Comparison with simulated data - OK
library(stabledist)
n = 100
nSimulations = 10000000
shapeValues = runif(n, 1.5, 2) # the mean of simulated values is more well behaved.
error = rep(NA,n)

for(i in 1:n)
{
  z = stabledist::rstable(n=nSimulations,alpha=shapeValues[i],beta = 0,
                          gm=1,delta=0,pm=1)
  error[i] = abs((mean(abs(z)) - 
                    .stableS1SymmetricMomentGarch(shapeValues[i]))/
                   .stableS1SymmetricMomentGarch(shapeValues[i]))
}
plot(error, ylab = "error In percentage")
summary(error)


# Comparison with the numerical integration routine
library(stabledist)
n = 20
shapeValues = runif(n, 1, 2)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
error = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz("dstable",gm = 0, delta = 1,
                                                    alpha = shapeValues[i], beta = 0,pm = 1)[1])
  functionValues[i] = .stableS1SymmetricMomentGarch(shapeValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(shapeValues,trueValues,functionValues)
summary(error)







# ------------------------------------------------------------------------------
# Test Cases for functions .stableS1MomentPowerGarch (Mittnik et al. (2002))
# ------------------------------------------------------------------------------

# E(|z|^1) = E(|x|^1) where x ~ N(0,sd = sqrt(2)) see Mittnik et al. (2002) - OK
.stableS1MomentPowerGarch(shape = 2, skew = 0, delta = 1)
.stableS1SymmetricMomentGarch(shape = 2)
sqrt(2)*.normMomentAparch (delta = 1, gm = 0)

# Comparison with symmetric stable garch model - OK
n = 10000
shapeValues = runif(n, 1, 2)
error = rep(NA,n)
for(i in 1:n)
{
  error[i] = abs((.stableS1MomentPowerGarch(shape=shapeValues[i],skew=0,delta=1) - 
                    .stableS1SymmetricMomentGarch(shape=shapeValues[i]) )/
                   .stableS1SymmetricMomentGarch(shape=shapeValues[i]))
}
plot(error, ylab = "error In percentage", ylim = c(0,1e-15))
summary(error)

# Comparison with simulated data - OK
# Works really well for all values with errors around 0.002%!!!
# The trueValues are the estimated moments using numerical integration 
# from function .trueAparchMomentsWurtz addapted from function .truePersistence
# on the garch-Stats.R file.
n = 100
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz("dstable",gm = 0, delta = deltaValues[i],
                                                    alpha = shapeValues[i], beta = skewValues[i],pm = 1)[1])
  functionValues[i] = .stableS1MomentPowerGarch(shape = shapeValues[i], skew = skewValues[i],
                                                delta = deltaValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(skewValues,deltaValues,shapeValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for functions .stableS1MomentAparch (GEVStableGarch papper)
# ------------------------------------------------------------------------------

# E(z^1) = E(x^1), where x ~ N(0,2) see Mittnik et al. (2002) - OK
z = stabledist::rstable(n=10^7,alpha=1.999,beta = 0,
                        gm=1,delta=0,pm=0)
lambdaSim = mean((abs(z)-gm*z)^delta)
lambdaSim
.stableS1MomentAparch(shape = 1.999999,skew = skew,delta = delta,gm = gm)

# Comparison with .stableS1MomentPowerGarch function - OK
# from Mittnik et al. (2002)
n = 10000
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = .stableS1MomentPowerGarch(shape = shapeValues[i], skew = skewValues[i],
                                            delta = deltaValues[i])
  functionValues[i] = .stableS1MomentAparch(shape = shapeValues[i], skew = skewValues[i],
                                           delta = deltaValues[i], gm = 0)
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% GEVStableGarch eq. VS Mittnik)", type = "l", col = "red")
cbind(deltaValues,shapeValues,skewValues,trueValues,functionValues)
summary(error)

# Comparison with the numerical integration - OK
# Works really well for all values with errors around 0.005%!!!
n = 100
shapeValues = runif(n, 1+0.1, 2)
skewValues = runif(n,-1,1)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.05)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = as.numeric(.trueAparchMomentsWurtz("dstable",gm = gmValues[i], 
                                                    delta = deltaValues[i],alpha = shapeValues[i], beta = skewValues[i],
                                                    pm = 1)[1])
  functionValues[i] = .stableS1MomentAparch(shape = shapeValues[i], skew = skewValues[i],
                                           delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% TrueValues-numerical integration)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,skewValues,trueValues,functionValues)
summary(error)



# ------------------------------------------------------------------------------
# Test Cases for functions .stableS0MomentAparch (GEVStableGarch papper)
# ------------------------------------------------------------------------------
.stableS0MomentAparch(shape = 1.5, skew = 0.5, delta = 1.3, gm = 0)

# Comparison with simulated data

n = 100
nSimulations = 100000
shapeValues = runif(n, 1.8, 2) # the mean of simulated values is more well behaved.
skewValues = rep(0.5,n)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.05)
error = rep(NA,n)
integrationValue = rep(NA,n)
integrationValueS1 = rep(NA,n)
simulatedValue= rep(NA,n)
for(i in 1:n)
{
  z = stable::rstable(n=nSimulations,alpha=shapeValues[i],
		beta =skewValues[i], param = 0)
  simulatedValue[i] = mean( ( abs(z) - gmValues[i] * z ) ^ deltaValues[i] )
  integrationValue[i] = .stableS0MomentAparch(shape = shapeValues[i], 
         skew = skewValues[i], delta =  deltaValues[i], gm = gmValues[i])
  integrationValueS1[i] = .stableS1MomentAparch(shape = shapeValues[i], 
         skew = skewValues[i], delta =  deltaValues[i], gm = gmValues[i])
 
}
error1 = 100 * abs( simulatedValue - integrationValue ) / simulatedValue 
error2 = 100 * abs( simulatedValue - integrationValueS1 ) / simulatedValue 
plot(error, type = 'l')
cbind(simulatedValue, integrationValue, integrationValueS1)
plot(simulatedValue, type = 'l')
lines(integrationValue, col = 2)
lines(integrationValueS1, col = 3)


sum(error1)
sum(error2)



# ------------------------------------------------------------------------------
# Test Cases for functions .stableS2MomentAparch (GEVStableGarch papper)
# ------------------------------------------------------------------------------

# Comparison with simulated data

n = 100
nSimulations = 100000
shapeValues = runif(n, 1.05, 1.5) # the mean of simulated values is more well behaved.
skewValues = rep(-0.5,n)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.05)
error = rep(NA,n)
integrationValue = rep(NA,n)
integrationValueS1 = rep(NA,n)
simulatedValue= rep(NA,n)
for(i in 1:n)
{
  z = stable::rstable(n=nSimulations,alpha=shapeValues[i],
		beta =skewValues[i], param = 2)
  simulatedValue[i] = mean( ( abs(z) - gmValues[i] * z ) ^ deltaValues[i] )
  integrationValue[i] = .stableS0MomentAparch(shape = shapeValues[i], 
         skew = skewValues[i], delta =  deltaValues[i], gm = gmValues[i])
  integrationValueS1[i] = .stableS1MomentAparch(shape = shapeValues[i], 
         skew = skewValues[i], delta =  deltaValues[i], gm = gmValues[i])
 
}
error1 = 100 * abs( simulatedValue - integrationValue ) / simulatedValue 
error2 = 100 * abs( simulatedValue - integrationValueS1 ) / simulatedValue 
plot(error, type = 'l')
cbind(simulatedValue, integrationValue, integrationValueS1)
plot(simulatedValue, type = 'l')
lines(integrationValue, col = 2)
lines(integrationValueS1, col = 3)


sum(error1)
sum(error2)









# ------------------------------------------------------------------------------
# Test Cases for function .stableS1SymmetricMomentAparch (Diongue papper)
# ------------------------------------------------------------------------------

# Comparison with stableS1.moment.garch function 
# from GEVStableGarch papper on JSS - OK
# Error around 1e-9%!!!
n = 100000
shapeValues = runif(n, 1+0.1, 2)
gmValues = runif(n,-1,1)
deltaValues = rep(NA,n)
for(i in 1:n) 
  deltaValues[i] = runif(n=1,1,shapeValues[i]-0.09)
trueValues = rep(NA,n)
functionValues = rep(NA,n)
for(i in 1:n)
{
  trueValues[i] = .stableS1MomentAparch(shape = shapeValues[i], skew = 0,
                                       delta = deltaValues[i], gm = gmValues[i])
  functionValues[i] = .stableS1SymmetricMomentAparch(shape = shapeValues[i],
                                                     delta = deltaValues[i], gm = gmValues[i])
}
error = 100*abs((trueValues - functionValues)/trueValues)
plot(error, main = "Error (% GEVStableGarch eq. VS Mittnik)", type = "l", col = "red")
cbind(deltaValues,shapeValues,gmValues,trueValues,functionValues)
summary(error)




# ------------------------------------------------------------------------------
# Test Cases for function moment.Aparch
# ------------------------------------------------------------------------------

# GARCH(1,0)
# GARCH(1,1)
# GARCH(5,1)
# APARCH(1,0)
# APARCH(5,1)

# GARCH(1,0) Test Cases
# GARCH(1,0)-gev
spec <- gsSpec(model = list(alpha = 1.3), 
                    presample = NULL,cond.dist = c("gev"), rseed = 3)
# GARCH(1,0)-stableS1
spec <- gsSpec(model = list(alpha = 1.3, delta = 1, shape = 1.5, skew = 0), 
                    presample = NULL,cond.dist = c("stableS1"), rseed = 3)
# GARCH(1,0)-gat
spec <- gsSpec(model = list(alpha = 1.3, delta = 2, shape = c(3,3)), 
                    presample = NULL,cond.dist = c("gat"), rseed = 3)
# GARCH(1,0)-norm
spec <- gsSpec(model = list(alpha = 1.3, delta = 2), 
                    presample = NULL,cond.dist = c("norm"), rseed = 3)
# GARCH(1,0)-std
spec <- gsSpec(model = list(alpha = 1.3, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("std"), rseed = 3)
# GARCH(1,0)-sstd
spec <- gsSpec(model = list(alpha = 1.3, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("sstd"), rseed = 3)
# GARCH(1,0)-ged
spec <- gsSpec(model = list(alpha = 1.3, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("ged"), rseed = 3)


# GARCH(1,1) Test Cases
# GARCH(1,1)-gev
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4), 
                    presample = NULL,cond.dist = c("gev"), rseed = 3)
# GARCH(1,1)-stableS1
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 1, shape = 1.5), 
                    presample = NULL,cond.dist = c("stableS1"), rseed = 3)
# GARCH(1,1)-gat
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 2, shape = c(3,3)), 
                    presample = NULL,cond.dist = c("gat"), rseed = 3)
# GARCH(1,1)-norm
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 2), 
                    presample = NULL,cond.dist = c("norm"), rseed = 3)
# GARCH(1,1)-std
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("std"), rseed = 3)
# GARCH(1,1)-sstd
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("sstd"), rseed = 3)
# GARCH(1,1)-ged
spec <- gsSpec(model = list(alpha = 1.3, beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("ged"), rseed = 3)


# GARCH(5,1) Test Cases
# GARCH(5,1)-gev
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4), 
                    presample = NULL,cond.dist = c("gev"), rseed = 3)
# GARCH(5,1)-stableS1
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 1, shape = 4), 
                    presample = NULL,cond.dist = c("stableS1"), rseed = 3)
# GARCH(5,1)-gat
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 2, shape = c(3,3)), 
                    presample = NULL,cond.dist = c("gat"), rseed = 3)
# GARCH(5,1)-norm
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 2), 
                    presample = NULL,cond.dist = c("norm"), rseed = 3)
# GARCH(5,1)-std
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("std"), rseed = 3)
# GARCH(5,1)-sstd
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("sstd"), rseed = 3)
# GARCH(5,1)-ged
spec <- gsSpec(model = list(alpha = runif(5,0,1), beta = 1.4, delta = 2, shape = 4), 
                    presample = NULL,cond.dist = c("ged"), rseed = 3)


# APARCH(1,0) Test Cases
# APARCH(1,0)-gev
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1)), 
                    presample = NULL,cond.dist = c("gev"), rseed = 3)
# APARCH(1,0)-stableS1
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1),
                                 shape = runif(1,0,5), skew = runif(1,0,5) ), 
                    presample = NULL,cond.dist = c("stableS1"), rseed = 3)
# APARCH(1,0)-gat
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1),
                                 shape = runif(2,0,5), skew = runif(1,0,5)), 
                    presample = NULL,cond.dist = c("gat"), rseed = 3)
# APARCH(1,0)-norm
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1)), 
                    presample = NULL,cond.dist = c("norm"), rseed = 3)
# APARCH(1,0)-std
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1),
                                 shape = runif(1,0,100)), 
                    presample = NULL,cond.dist = c("std"), rseed = 3)
# APARCH(1,0)-sstd
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(5,-1,1), 
                                 shape = runif(1,0,5), skew = runif(1,0,5)), 
                    presample = NULL,cond.dist = c("sstd"), rseed = 3)
# APARCH(1,0)-ged
spec <- gsSpec(model = list(alpha = runif(5,1,3), delta = runif(1,0,5), gm = runif(1,-1,1), shape = 4), 
                    presample = NULL,cond.dist = c("ged"), rseed = 3)


# APARCH(5,1) Test Cases
# APARCH(5,1)-gev
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1)), 
                    presample = NULL,cond.dist = c("gev"), rseed = 3)
# APARCH(5,1)-stableS1
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1),
                                 shape = runif(1,0,5), skew = runif(1,0,5) ), 
                    presample = NULL,cond.dist = c("stableS1"), rseed = 3)
# APARCH(5,1)-gat
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1),
                                 shape = runif(2,0,5), skew = runif(1,0,5)), 
                    presample = NULL,cond.dist = c("gat"), rseed = 3)
# APARCH(5,1)-norm
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1)), 
                    presample = NULL,cond.dist = c("norm"), rseed = 3)
# APARCH(5,1)-std
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1),
                                 shape = runif(1,0,100)), 
                    presample = NULL,cond.dist = c("std"), rseed = 3)
# APARCH(5,1)-sstd
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1), 
                                 shape = 3, skew = 2), 
                    presample = NULL,cond.dist = c("sstd"), rseed = 3)
# APARCH(5,1)-ged
spec <- gsSpec(model = list(alpha = runif(5,0,3), delta = runif(1,0,5), gm = runif(5,-1,1), shape = 4), 
                    presample = NULL,cond.dist = c("ged"), rseed = 3)


# General testing
.stationarityAparch(model = list(alpha = spec@model$alpha, beta = spec@model$beta, gm = spec@model$gm, 
                                 delta = spec@model$delta, skew = spec@model$skew, shape = spec@model$shape), 
                    formula = .getFormula(spec@formula), cond.dist = spec@distribution)



# ------------------------------------------------------------------------------
# Test Cases for function gsMomentAparch
# ------------------------------------------------------------------------------

gsMomentAparch(cond.dist = "stableS1", shape = 1.1, skew = 0, delta = 1.01, gm = 0.99999)

gsMomentAparch(cond.dist = "gev", shape = -4, skew = 0, delta = 1.4, gm = 0)

gsMomentAparch(cond.dist = "gat", shape = c(1.9,2.3), skew = 0.5, delta = 0.4, gm = 0)

gsMomentAparch(cond.dist = "norm", shape = c(1.9,2.3), skew =1, delta = 11.4, gm = -0.999)

gsMomentAparch(cond.dist = "std", shape = 2.001, skew = -0.5, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "sstd", shape = 2.001, skew = 0.11, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "skstd", shape = 5.001, skew = 0.11, delta = 3, gm = -0.5)

gsMomentAparch(cond.dist = "ged", shape = 6, skew = 0.11, delta = 5.11, gm = -0.5)






################################################################################

