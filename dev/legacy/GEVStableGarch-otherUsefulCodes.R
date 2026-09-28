##############
# ARMA ESTIMATION
##############


# ------------------------------------------------------------------------------
# Fitting pure ARMA process 
# ------------------------------------------------------------------------------



# Notes:
# The "sigma" parameter is the square root of the sigma2 parameter estimated
# by function "arima".
# works really well for arma(1,1), arma(m,1), arma(0,n),  arma(1,n), arma(0,n)
library(fGarch)
library(Rsolnp)
library(FitARMA)
data(dem2gbp)
x = dem2gbp[, 1]+10
# arma(1,1)-norm-intercept-nlminb+nm
m <- 2
n <- 2

fit1 <- gsFit(data = x, 
              formula = as.formula(paste ("~ arma(",m,", ",n,")", sep = "", collapse = NULL)),
              cond.dist = "norm", include.mean = TRUE, 
              algorithm = "sqp", DEBUG = FALSE, control = list(trace = 3))
model1 <- arima(x, order = c(m, 0, n))
model2 <- FitARMA(x, order = c(m,0,n))
par.result <- cbind(model1$coef[c(m+n+1,1:(m+n))],model1$coef[c(m+n+1,1:(m+n))],fit1@fit$par[1:(1+m+n)])
colnames(par.result) = c("arima","FitARMA","gsFit")
par.result
absoluteError <- abs((fit1@fit$par[1:(1+m+n)]-model1$coef[c(m+n+1,1:(m+n))])/fit1@fit$par[1:(1+m+n)])
absoluteError
fit1@fit$llh
model1$loglik



# arma(1,1)-std
fit1 <- gsFit(data = x, formula = ~arma(1,1),
              cond.dist = "sstd", include.mean = TRUE, 
              algorithm = "nlminb+nm")

model1 <- gsFit(data = x, formula = ~arma(1,1),
                cond.dist = "std", include.mean = TRUE, 
                algorithm = "sqp")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)

# arma(0,1)-norm
fit1 <- gsFit(data = x, formula = ~arma(0,1),
              cond.dist = "norm", include.mean = TRUE, 
              algorithm = "nlminb+nm")

model1 <- gsFit(data = x, formula = ~arma(0,1),
                cond.dist = "norm", include.mean = TRUE, 
                algorithm = "sqp")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)


# aparch(2,2)-std
fit1 <- gsFit(data = x, formula = ~aparch(2,2),
              cond.dist = "std", include.mean = TRUE, 
              algorithm = "nlminb+nm")

model1 <- gsFit(data = x, formula = ~aparch(2,2),
                cond.dist = "std", include.mean = TRUE, 
                algorithm = "sqp")

model2 <- garchFit(data = x, formula = ~aparch(2,2),
                   cond.dist = "std", include.mean = TRUE, 
                   algorithm = "nlminb+nm")

abs((model1@fit$par-fit1@fit$par)/fit1@fit$par)
abs((model2@fit$par-model1@fit$par)/model1@fit$par)





# ------------------------------------------------------------------------------
# Arma likelihood code inside the llh function of gsFit
# ------------------------------------------------------------------------------
# BEGIN: Log Likelihood for pure ARMA process
#############################################
armaLLH <- function(parm)
{      
  if(DEBUG)
    print(parm)
  # check if some parameters are NAN
  if(sum(is.nan(parm)) != 0) {return(1e99)}
  
  # Getting parameters from parm vector 
  mu <- parm[1];
  a <- parm[(1+1):(2+m-1)]
  b <- parm[(1+m+1):(2+m+n-1)]
  skew <- parm[1+m+n+1]
  shape <- parm[(2+m+n+1):(2+m+n+lengthShape)]
  sigma <- parm[2+lengthShape+m+n+1]
  
  # Setting parameters to accept tapper off MA, AR or GARCH coefficients
  if( AR == TRUE) 
    a <- 0
  if( MA == TRUE) 
    b <- 0
  if (include.mean == FALSE)
    mu <- 0
  
  if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
  {
    if( shape-delta < tolerance$TOLSTABLE || abs(shape) < tolerance$TOLSTABLE ||
          !(abs(skew) < 1) || !((shape - 2) < 0) )
    {
      return(1e99)
    }
  }
  
  # ARMA stationarity condition check
  if(!arCheck(a))
    return(1e99)
  
  # Stop if parametr is not > 0
  if(sigma <= 0)
    return(1e99)
  
  # Filters the Time series to obtain the i.i.d. sequence of 
  # 'innovations' to evaluate the Log-likelihood function
  z <- .filterArma(data = data, size.data = N, m = m, n = n, mu = mu, a = a, b = b)
  
  if(DEBUG)
  {
    print(c("length(z)",length(z)))
    print(c("N",N))
  }
  
  
  # get output Residuals
  if (optim.finished)
  {
    out$residuals <<- as.numeric(z)
    out$sigma.t <<- sigma
    out$h.t <<- sigma
  }
  
  # Return llh function        
  llh.dens <- .armaDist(z = z, sigma = sigma, shape = shape, 
                        skew = skew, cond.dist = cond.dist)
  llh <- llh.dens
  if (is.nan(llh) || is.infinite(llh) || is.na(llh)) 
  {
    llh <- 1e99
  }
  llh
}
# END: Log Likelihood for pure ARMA process
#############################################










# ------------------------------------------------------------------------------
# Old .getStart function
# ------------------------------------------------------------------------------





.getStartOld <- function(data,m,n,p,q, AR = FALSE, MA = FALSE, ARMAonly = FALSE,
                         cond.dist = c("stable", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
                         TOLG = 1e-7, TOLSTABLE = 2e-2)
{    
  # Description:
  #   Get initial values to start the estimation and the bounds for
  #   parameters to be estimated inside GSgarch.Fit function.
  #   Remarks: This function tunes initial parameters to perform optimization
  #   The ARMA coefficients are the ones returned by the "arima" function
  #   adjusted parameters functions from package ("arima" belongs to package "stats" from R)
  #   For GARCH(p,q) Wurtz et al. (2006) suggested
  #   using omega = 0.1, gamma = 0, alpha = 0.1/p and beta = 0.8/q
  #   delta is chosen to be initially equals to 1.5 in almost all the cases 
  #   Keep in mind that delta < alpha for the stable case.
  #   The arma order passed to this function will be 1 even though the 
  #   process does not possess the AR or MA component. Therefore, 
  #   it is not a problem to have AR our MA input parameters equal to FALSE 
  #   even thought the model order says on the contrary.
  
  # Arguments:
  #   data - vector of data
  #   m, n, p, q - model order as in ARMA(m,n)-GARCH/APARCH(p,q)
  #   AR - boolean value that indicates whether we have a model
  #   with th Autoregressive part included
  #   MA - boolean value that indicates whether we have a model
  #   with the Moving Average part included
  #   ARMAonly - Indicates whether we have a pure ARMA model
  #   cond.dist - name of the conditional distribution, one of
  #       gev, stable, norm, std, sstd, ged
  #   TOLSTABLE - boundary tolerance. Should be greater than GSstable.tol
  #   TOLG - pper and lower bounds tolerance. Should be greater than tol
  
  # Return:
  #   Asdf - The asdf 
  
  # FUNCTION:  
  
  # Error control
  if (m < 0 || n < 0 || m %% 1 != 0 || n %% 1 != 0)
    stop("'m' and 'n' need to be integers greater than zero")
  if (m == 0 || n == 0)
    stop("GSgarch.GetStart expects 'm' and 'n' different from zero")
  if ( (m != 1 && AR == TRUE)  || (n != 1 && MA == TRUE) )
    stop("If AR = TRUE, 'm' should be 1 and if MA = TRUE 'n' should be 1")
  if (p < 0 || q < 0 || p %% 1 != 0 || q %% 1 != 0)
    stop("'p' and 'q' need to be integers greater than zero")
  if(p == 0 && q != 0)
    stop("Invalid Garch(p,q) order")
  # Conditional distribution
  cond.dist = match.arg(cond.dist)
  if( !is.numeric(data) || !is.vector(data))
    stop("data set must be a numerical one dimensional vector")
  
  # Initial ARMA parameters
  Mean <- mean(data)
  Var <- var(data)
  Dispersion <- mean(abs(x-Mean))
  arima.fit <- c()
  arima.failed <- FALSE
  arima.fit.try <- "empty"  
  arima.m <- m
  arima.n <- n
  if(AR == TRUE) # we don't have the AR part
    arima.m <- 0
  if(MA == TRUE) # we don't have the MA part 
    arima.n <- 0
  
  if(ARMAonly == FALSE) # try 'arima' function only if we are not in a pure 'arma' model
  {
    # try arima fit with correct order depending on AR and MA
    try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
    if( is.numeric(arima.fit.try) )
      arima.fit <- arima.fit.try
    else
      arima.failed <- TRUE
    
    # Configuring initial ar,ma and mean parameters returned by 'arima' function.
    if (arima.failed == FALSE) 
    {
      if(AR == FALSE && MA == FALSE)
      {
        ar.init <- arima.fit[1:arima.m]
        ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
      }
      if(AR == TRUE && MA == FALSE)
      {
        ar.init <- 0
        ma.init <- arima.fit[(arima.m+1):(arima.m+arima.n)]
      }
      if(AR == FALSE && MA == TRUE)
      {
        ar.init <- arima.fit[1:arima.m]
        ma.init <- 0
      }
      if(AR == TRUE && MA == TRUE)
      {
        ar.init <- 0
        ma.init <- 0
      }
      mean.init <- arima.fit[arima.m+arima.n+1]
    }
    else # arima function failed
    {
      mean.init <- Mean
      ar.init <- rep(0,m)
      ma.init <- rep(0,n)
    }
  } else {
    mean.init <- Mean
    ar.init <- rep(0,m)
    ma.init <- rep(0,n)        
  }
  
  # Initial APARCH and Density parameters
  omega.init <- min(0.1,0.1*Var)
  alpha.init <- rep(0.1/p,p)
  beta.init <- rep(0.8/q,q)
  shape.init <- 1
  skew.init <- 1
  gm.init <- rep(0,p)
  delta.init <- 1.5
  sigma.init <- Var
  
  # Lower Bound
  #mean.lower <- mean.init - 3*abs(mean.init)
  mean.lower <- -10*abs(mean.init)
  arma.lower <- rep(-10+TOLG,m+n)
  omega.lower <- rep(TOLG,1)
  alpha.lower <- rep(TOLG,p)
  beta.lower <- rep(TOLG,q)
  shape.lower <- 0
  skew.lower <- 0
  gm.lower <- rep(-1 + TOLG,p)
  delta.lower <- TOLG
  sigma.lower <- 0
  
  # Upper Bound  
  #mean.upper <- mean.init + 3*abs(mean.init)
  mean.upper <- 10*abs(mean.init)
  arma.upper <- rep(10-TOLG,m+n)
  omega.upper <- rep(1-TOLG,1)
  #omega.upper <- 100*var(data)
  alpha.upper <- rep(1-TOLG,p)
  beta.upper <- rep(1-TOLG,q)
  shape.upper <- 4
  skew.upper <- 4
  gm.upper <- rep(1 - TOLG,p)
  # delta.upper <- 3 + TOLG
  delta.upper <- 30 + TOLG
  sigma.upper <- 10*Var
  
  # Setting skew and shape appropriatelly for other conditional Distributions
  if (cond.dist == "std")
  {   
    shape.init <- 4
    shape.lower <- 2 + TOLG; shape.upper <- 20   
  }
  if (cond.dist == "sstd")
  {   
    shape.init <- 4; skew.init = 1
    shape.lower <- 2 + TOLG; shape.upper <- 20
    skew.lower <- TOLG; skew.upper <- 20 
  }
  if (cond.dist == "skstd")
  {   
    shape.init <- 4; skew.init = 1
    shape.lower <- TOLG; shape.upper <- 100
    skew.lower <- TOLG; skew.upper <- 100   
  }
  if (cond.dist == "ged")
  {   
    shape.init <- 4
    shape.lower <- 0 + TOLG; shape.upper <- 20   
  }
  if (cond.dist == "gat")
  {   
    shape.init <- c(4,1)
    shape.lower <- rep(0 + TOLG,2)
    shape.upper <- rep(100,2)   
  }
  if(cond.dist == "gev")
  {
    mean.init <- 0
    shape.init <- 0.1; 
    alpha.init <- rep(0.05/p,p)
    beta.init <- rep(0.8/q,q)
    # shape.lower <- -0.5 + TOLG; shape.upper <- 0.5 - TOLG
    shape.lower <- -10 + TOLG; shape.upper <- 10 - TOLG
  }
  if(cond.dist == "stable")
  {
    shape.init <- 1.8; delta.init <- 1
    shape.lower <- 1 + TOLSTABLE; shape.upper <- 2 - TOLSTABLE
    skew.init <- 0
    skew.lower <- -1 + TOLSTABLE; skew.upper <- 1 - TOLSTABLE
    delta.upper <- 1.9
  }
  
  # Constructing the initial, lower and upper bound vectors
  arma.init <- c(ar.init,ma.init)
  if(!ARMAonly)
  {
    init <- c(mean.init,arma.init,omega.init,alpha.init,gm.init,
              beta.init,delta.init,skew.init,shape.init)
    lower <- c(mean.lower,arma.lower,omega.lower,alpha.lower,gm.lower,
               beta.lower,delta.lower,skew.lower,shape.lower)
    upper <- c(mean.upper,arma.upper,omega.upper,alpha.upper,gm.upper,
               beta.upper,delta.upper,skew.upper,shape.upper)
    
    namesStart = c("mu", paste("ar", 1:length(ar.init), sep = ""), paste("ma", 1:length(ma.init), sep = ""),
                   "omega", paste("alpha", 1:length(alpha.init), sep = ""),paste("gm", 1:length(gm.init), sep = ""), 
                   paste("beta", 1:length(beta.init), sep = ""), "delta","skew",paste("shape", 1:length(shape.init), sep = ""))        
    
  } else {
    init <- c(mean.init,arma.init,skew.init,shape.init,sigma.init)
    lower <- c(mean.lower,arma.lower,skew.lower,shape.lower,sigma.lower)
    upper <- c(mean.upper,arma.upper,skew.upper,shape.upper,sigma.upper)   
    namesStart = c("mu", paste("ar", 1:length(ar.init), sep = ""), paste("ma", 1:length(ma.init), sep = ""),
                   "skew",paste("shape", 1:length(shape.init), sep = ""),"sigma")
  }
  
  if(arima.failed == TRUE)
    warning("arima function from package stats failed to get initial AR and MA coefficients")
  
  # Create result
  result = rbind(init,lower,upper)
  colnames(result) = namesStart
  
  # Return
  result
}




