
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
# TEST CASES FOF FUNCTIONS:               SPECIFICATION:
#  .filterAparch                    Test how our filter.Aparh function
#  .filterGarch11Fit               build with matrices operation is related  
# .filterAparchForLoop   				  to the other filter functions described
#                                   in Wuertz et al. (2006)
################################################################################

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1][1:100]

###############
# Comparison between our function and the filter function from Wurtz.
# garch(1,1)
filter1Result <- .filterAparch(data = x, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3), beta = c(0.5), gamma = c(0),delta = 2)

filter2Result <- .filterGarch11Fit(x,parm = c(mu = 0.4, omega = 0.1, alpha = 0.3, beta = 0.5))

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.5381661^2, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = c(0.3), beta = c(0.5), gamma = c(0),delta = 2)

cbind(filter1Result[,1],filter2Result[,1],filter3Result[,1],filter1Result[,2],filter2Result[,2],filter3Result[,2])

# arch(1) == garch(1,0): put p = 1, q = 1 and beta = 0. 
# Notice that the .filterAparchForLoop must be equal 
# to the other filter functions after a couple of steps. 
# This is due to the different starting points they have.
filter1Result <- .filterAparch(data = x, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3), beta = c(0), gamma = c(0),delta = 2)
filter2Result <- .filterGarch11Fit(x,parm = c(mu = 0.4, omega = 0.1, alpha = 0.3, beta = 0))

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.4136527^2, p = 1,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = c(0.3), beta = c(0), gamma = c(0),delta = 2)
cbind(filter1Result[,1],filter2Result[,1],filter3Result[,1],filter1Result[,2],filter2Result[,2],filter3Result[,2])

# garch(2,2)
filter1Result <- .filterAparch(data = x, p = 2,q = 2,mu = 0.4, omega = 0.1, 
                               alpha = c(0.3,0.1), beta = c(0.1,0.2), gamma = c(0,0),delta = 2)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 2,q = 2,mu = 0.4, omega = 0.1, 
                                       alpha = c(0.3,0.1), beta = c(0.1,0.2), gamma = c(0,0),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])


# arch(15) == garch(15,0)
filter1Result <- .filterAparch(data = x, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)
a <- NULL; 
a <- 0.2
is.null(a)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# arch(15) == garch(15,0)
filter1Result <- .filterAparch(data = x, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                               alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 15,q = 1,mu = 0.4, omega = 0.1, 
                                       alpha = rep(0.03,15), beta = c(0), gamma = rep(0,15),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(4,5)
filter1Result <- .filterAparch(data = x, p = 4,q = 5,mu = -5, omega = 5, 
                               alpha = rep(0.05,4), beta = c(0.1,0.03,0.02,0.005,0.02), gamma = rep(0.9,4),delta = 1.1)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 4,q = 5,mu = -5, omega = 5, 
                                       alpha = rep(0.05,4), beta = c(0.1,0.03,0.02,0.005,0.02), gamma = rep(0.9,4),delta = 1.1)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(1,1)
filter1Result <- .filterAparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = rep(0.3), beta = c(0.1), gamma = rep(0),delta = 2)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 1,q =1 ,mu = -5, omega = 5, 
                                       alpha = rep(0.3), beta = c(0.1), gamma = rep(0),delta = 2)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(2,2)
filter1Result <- .filterAparch(data = x, p = 2,q =2 ,mu = -5, omega = 5, 
                               alpha = rep(0.1,2), beta = c(0.1,0.3), gamma = rep(0,2),delta = 1.1)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 2,q =2 ,mu = -5, omega = 5, 
                                       alpha = rep(0.1,2), beta = c(0.1,0.3), gamma = rep(0,2),delta = 1.1)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

# aparch(1,0)
filter1Result <- .filterAparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

filter3Result <- .filterAparchForLoop(data = x, h.init = 0.1, p = 1,q =1 ,mu = -5, omega = 5, 
                                       alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

cbind(filter1Result[,1],filter3Result[,1],filter1Result[,2],filter3Result[,2])
filter1Result[,2]-filter3Result[,2]

###############
# Testing different starting points in our filter function

# arch(1) == garch(1,0)
filter1Result <- .filterAparch(data = x, p = 1,q =1 ,mu = -5, omega = 5, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

filter2Result <- filter.Aparch(data = x, p = 1,q =1 ,mu = -5, omega = 5,init.Value = 100, 
                               alpha = 0.2, beta = 0 , gamma = 0, delta = 20)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]

# aparch(2,3)
filter1Result <- .filterAparch(data = x, p = 2,q =3 ,mu = -2, omega = 0.5, 
                               alpha = rep(0.1,2), beta = rep(0.04,3) , gamma = c(0.3,-0.8), delta = 1.2)

filter2Result <- filter.Aparch(data = x, p = 2,q =3 ,mu = -2, omega = 0.5,init.Value = 100, 
                               alpha = rep(0.1,2), beta = rep(0.04,3) , gamma = c(0.3,-0.8), delta = 1.2)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]

# garch(59,100): Difficult to achieve the same time series after some steps. We need to pay 
# attention that we are not going to start our time series deliberatelly. We may start it 
# at a place that is representative of the conditional variance.
x <- rnorm(100000)
filter1Result <- .filterAparch(data = x, p = 59,q =100 ,mu = -2, omega = 0.5, 
                               alpha = rep(5,59), beta = rep(0.04,100) , gamma = rep(0.5,59), delta = 2)

filter2Result <- .filterAparch(data = x, p = 59,q =100 ,mu = -2, omega = 0.5, init.Value = var((abs(x+2))^2),
                               alpha = rep(5,59), beta = rep(0.04,100) , gamma = rep(0.5,59), delta = 2)

cbind(filter1Result[,2],filter2Result[,2])
filter1Result[,2]-filter2Result[,2]
var((abs(x+2))^2)


################################################################################
# TEST CASES FOF FUNCTIONS:               SPECIFICATION:
#    filter.Ama                           Test how our filter 
#                     
################################################################################

library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1][1:1000]

# Try to enter invalid parameters
result <- .filterArma(data = x,m = 1,n = 1, mu = 0, a = c(3,3), b = 0)
result <- .filterArma(data = x,m = 1,n = 3, mu = 0, a = c(3,3), b = c(3,3))
result <- .filterArma(data = x,m = 1,n = 3, mu = 0, a = c(3,3), b = c(3,3))


# filter with no variables. Must be equal to the original time series
result <- .filterArma(data = x,m = 1,n = 1,mu = 0,a = 0,b = 0)
x-result

# filter with mean = 1. Must be equal to the original time series - 1
x = dem2gbp[, 1][1:1000]
result <- .filterArma(data = x,m = 1,n = 1,mu = 1,a = 0,b = 0)
x - 1 - result

# Create a random arma(1,1) process with random residuals and try 
# to rebuild the residuals vector with function .filterArma 
N = 10000; mu = 100; a = 0.3; b = 0.3; x.init = 0; z.init = 0;
z = rnorm(N)
x = rep(NA,N)
x[1] = mu + a*x.init + b*z.init + z[1]
for (t in 2:N)
  x[t] = mu + a*x[t-1] + b*z[t-1] + z[t]

result <- .filterArma(data = x,m = 1,n = 1,mu = mu,a = a,b = b)
cbind(result,z)
summary((result-z)/result)

# Filter a time series and then try to reconstruct the original
# time series with the result obtained from .filterArma function
data = dem2gbp[, 1]
z = .filterArma(data = data,m = 1,n = 1,mu = mu,a = a,b = b)
N = length(x); mu = 100; a = 0.3; b = 0.3; x.init = 0; z.init = 0;
x = rep(NA,N)
x[1] = mu + a*x.init + b*z.init + z[1]
for (t in 2:N)
  x[t] = mu + a*x[t-1] + b*z[t-1] + z[t]
cbind(data,x)



###############
###############
# Garch Fit function for Garch(1,1) process.
# Code snipet from Wuertz et al. (2006).

garch11Fit = function(x, start.h = 0.1)
{
  flag = TRUE
  # Step 1: Initialize Time Series Globally:
  x <<- x
  # Step 2: Initialize Model Parameters and Bounds:
  Mean = mean(x); Var = var(x); S = 1e-6
  params = c(mu = Mean, omega = 0.1, alpha = 0.1, beta = 0.8)
  lowerBounds = c(mu = -10*abs(Mean), omega = S^2, alpha = S, beta = S)
  upperBounds = c(mu = 10*abs(Mean), omega = 100*Var, alpha = 1-S, beta = 1-S)
  
  # Step 3: Set Conditional Distribution Function:
  garchDist = function(z, hh) { dnorm(x = z/hh)/hh }
  
  # Step 4: Compose log-Likelihood Function:
  garchLLH = function(parm) {
    mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
    z = (x-mu); Mean = mean(z^2)
    
    # Use Filter Representation:
    e = omega + alpha * c(Mean, z[-length(x)]^2)
    #h = filter(e, beta, "r", init = Mean)
    
    
    # Conditional Variance filtering - filter - Wuertz et al. (2006)
    gamma = 0
    delta  = 2
    eps = z
    uv = 1
    h <- rep(start.h, uv)
    edelta = (abs(eps)-gamma*eps)^delta
    edeltat = filter(edelta, filter = c(0, alpha), sides = 1)
    c = omega/(1-sum(beta))
    h = c( h[1:uv], c + filter(edeltat[-(1:uv)], filter = beta,
                               method = "recursive", init = h[uv:1]-c))    
    #     if(alpha == 0.1 && flag == TRUE)
    #     {
    #       print(h)
    #       flag = FALSE
    #     }
    
    
    
    
    hh = sqrt(abs(h))
    llh = -sum(log(garchDist(z, hh)))
    llh 
  }
  print(garchLLH(params))
  # Step 5: Estimate Parameters and Compute Numerically Hessian:
  fit = nlminb(start = params, objective = garchLLH,
               lower = lowerBounds, upper = upperBounds, control = list(trace=3))
  epsilon = 0.0001 * fit$par
  Hessian = matrix(0, ncol = 4, nrow = 4)
  for (i in 1:4) {
    for (j in 1:4) {
      x1 = x2 = x3 = x4 = fit$par
      x1[i] = x1[i] + epsilon[i]; x1[j] = x1[j] + epsilon[j]
      x2[i] = x2[i] + epsilon[i]; x2[j] = x2[j] - epsilon[j]
      x3[i] = x3[i] - epsilon[i]; x3[j] = x3[j] + epsilon[j]
      x4[i] = x4[i] - epsilon[i]; x4[j] = x4[j] - epsilon[j]
      Hessian[i, j] = (garchLLH(x1)-garchLLH(x2)-garchLLH(x3)+garchLLH(x4))/
        (4*epsilon[i]*epsilon[j])
    }
  }
  
  # Step 6: Create and Print Summary Report:
  se.coef = sqrt(diag(solve(Hessian)))
  tval = fit$par/se.coef
  matcoef = cbind(fit$par, se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(matcoef) = list(names(tval), c(" Estimate",
                                          " Std. Error", " t value", "Pr(>|t|)"))
  cat("\nCoefficient(s):\n")
  printCoefmat(matcoef, digits = 6, signif.stars = TRUE)
  fit$par
}





###############
###############
# ARMA(1,1) function.

arma11Fit = function(x)
{
    # Step 1: Initialize Time Series Globally:
    x <<- x
    
    # Step 2: Initialize Model Parameters and Bounds: 
    Mean = mean(x); Var = var(x); S = 1e-6  
    params = c(mu = Mean, a = 0, b = 0, sigma = sqrt(Var))
    lowerBounds = c(mu = -10*abs(Mean), a = -10, b = -10, sigma = S)
    upperBounds = c(mu =  10*abs(Mean), a =  10, b =  10, 20*Var)
    
    # Step 3: Set Conditional Distribution Function:
    armaDist = function(z, hh) { dnorm(x = z/hh)/hh }
    
    # Step 4: Compose log-Likelihood Function:
    armaLLH = function(parm) {
      
      mu = parm[1]; a = parm[2]; b = parm[3]; sigma = parm[4]
      z = (x-mu);
      
      # Use Filter Representation:
      e = z - a*c(0, z[-length(x)])
      h = filter(e, -b, "r", init = 0)
      
      llh = -sum(log(armaDist(h, sigma)))
      
      llh 
    }
  
    #Step 5: Estimate Parameters and Compute Numerically Hessian:
    fit = nlminb(start = params, objective = armaLLH,
                lower = lowerBounds, upper = upperBounds, control = 
                  list(trace=3))
    
    
    fit$par
    -sum(log(armaDist(h, fit$par[4])))
    fit$residuals = .filterArma(x,1,1,mu = fit$par[1], a = fit$par[2],b = fit$par[3])
    fit

}

library(FitARMA)
library(fGarch)
data(dem2gbp)
x = (dem2gbp[, 1]+30)
fitMeu <- arma11Fit(x) 
?arima
fitArima <- arima(x, order = c(1,0,1), optim.control = list(trace = 0))
cbind(fitMeu$residuals,fitArima$residuals)
cbind(fitMeu$par[1:3],fitArima$coef[c(3,1,2)])

.filterArma(x,1,1,mu = 29.9836, a = -0.6229,b = 0.6446)
