#############################################################################
#  gsStableGarch.R
#  Self-contained ARMA-GARCH/APARCH fitting + simulation with stableS0
#
#  Based on GEVStableGarch package v1.1 by Thiago do Rego Sousa.
#  Stable density computed via libstable4u (replaces stabledist::dstable).
#
#  Usage:  source("gsStableGarch.R")
#    then call gsFit(...) or gsSim(...) with cond.dist = "stableS0"
#
#  IMPORTANT: Only cond.dist = "stableS0" is supported.
#  The stableS1 and stableS2 parametrizations appear in the function
#  signatures for backward compatibility, but the density is always
#  computed via libstable4u::stable_pdf, which uses the S0 (pm = 0)
#  parametrization. Using stableS1 or stableS2 will produce incorrect
#  results because the 'param' argument is silently ignored.
#############################################################################

# --- Auto-install required packages if missing ---
.required_pkgs <- c("Rsolnp", "fGarch", "stabledist", "skewt",
                   "timeDate", "timeSeries", "libstable4u", "methods")
.to_install <- .required_pkgs[!sapply(.required_pkgs, requireNamespace, quietly = TRUE)]
if (length(.to_install) > 0) {
  message("Installing packages: ", paste(.to_install, collapse = ", "))
  install.packages(.to_install, dependencies = TRUE)
}
rm(.required_pkgs, .to_install)

# --- Load packages ---
library(methods)
library(Rsolnp)
library(fGarch)
library(stabledist)
library(skewt)
library(timeDate)
library(timeSeries)
library(libstable4u)

#############################################################################
# Source code below was merged from the GEVStableGarch package v1.1 files.
# The .GSgarch.dstable function uses libstable4u::stable_pdf for speed.
#############################################################################




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
# FUNCTION:                           DESCRIPTION:
#  'GEVSTABLEGARCH'                  GEVSTABLEGARCH Class representation
################################################################################


# Class Representation:
setClass("GEVSTABLEGARCH",
         representation(
           call = "call",
           formula = "formula",
           method = "character",
	         convergence = "numeric",
	         messages = "list",
           data = "numeric",
           fit = "list",
           residuals = "numeric",
           h.t = "numeric",
           sigma.t = "numeric",
           title = "character",
           description = "character")
)


################################################################################





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
# FUNCTION:               SPECIFICATION: 
#  GEVSTABLEGARCHSPEC     S4 GEVSTABLEGARCHSPEC Class representation 
################################################################################


# S4 GEVStableGarchSpec Class
setClass("GEVSTABLEGARCHSPEC",
representation(
         call = "call",
        formula = "formula",
         model = "list",
        presample = "matrix",
        distribution = "character",
        rseed = "numeric")
)


################################################################################




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
#  .filterArma             Filter function for Arma, Ar or Ma process
#                          
################################################################################


.filterArma <- function(
  data, size.data,
  m,n, 
  mu, a, b)
{  
    # Arguments
    # data: vector with data
    # m,n: model order arma(m,n)
    # mu,a,b: model parameters. 'a' is the autorregressive coefficient 
    # and 'b' is the moving average coeficient.
    # REMARKS: This function filters, pure 'arma' models, but 
    # it can deals with other models, such as 'ar' and 'ma' models. 
    # All the model parameters must be specified, regardless if they
    # exist. For example:
    # for ar(m) m >= 1 make n = 1 and b = 0.
    # for ma(n) n >= 1 make m = 1 and a = 0.
    
    # Return
    # the residuals 'z'
    
    # FUNCTION:
    
    # error treatment of input parameters
    if( n < 1 || n < 1 || length(a) != m || length(b) != n)
        stop("One or more of these conditions were true:
           n < 1 || n < 1 || length(a) != m || length(b) != n")
    
    # Initial declaration of variables
    x.init <- rep(0,m)
    z.init <- rep(0,n) 
    
    N <- size.data
    x.ZeroMean <- data - mu
    x2 = 0
    V <- c(x.init,x.ZeroMean[1:(N-1)])
    for( i in 1:m)
    {
      x1 <- -a[i]*V
      x2 = x2 +  x1[(m-(i-1)):(m+N-i)]
    }
    x2 <- x.ZeroMean + x2
    
    z <- filter(x2, filter = -b,
                method = "recursive", init = z.init)        
    
    if(length(z) != N)
        stop("Error in filtering function. length(z) != N")    
    
    # return
    return(z)
}


################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparch          Filtering Aparch process and its particular cases
#                          
################################################################################

.filterAparch <- function(
  data,init.Value = NULL,
  p,q, 
  mu, omega, alpha, beta, gamma, delta)
{
  
    # Arguments
    # data: vector with data
    # p,q: model order. aparch(p,q)
    # mu,omega,alpha,beta,gamma,delta: model parameters
    # REMARKS: This function filters, pure 'aparch' models, but 
    # it can deals with other models, such as garch, arch. 
    # All the model parameters must be specified, regardless if they
    # exist. For example:
    # for garch(p,q) with p,q >= 1 make gamma = 0.
    # for arch(p) p >= 1 make q = 1, beta = 0 and gamma = 0.
    
    # Return
    # the series 'z' and 'h'
  
  
    # input 
    # data: vector with data
  
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
        stop("One or more of these conditions were true:
          p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
          || length(beta) != q")
    
    # Initial declaration of variables
    pq = max(p,q)
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    
    # initializing the time series
    if(is.null(init.Value)) {
        edelta.init <- rep(Mean.z,p)
        h.init <- rep(Mean.z,q)
    } else {
        edelta.init <- rep(init.Value,p)
        h.init <- rep(init.Value,q)  
    }  
    
    edeltat = 0
    for( i in 1:p)
    {
      edelta <- alpha[i]*(c(edelta.init,((abs(z)-gamma[i]*z)^delta)[1:(N-1)]))
      edeltat = edeltat +  edelta[(p-(i-1)):(p+N-i)]
    }
    edeltat = omega + edeltat
    
    h <- filter(edeltat, filter = beta,
                method = "recursive", init = h.init)
    
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)
}

################################################################################









################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterGarch11Fit      Filter function for Garch(1,1) process. This code 
#                          snipet was taken from the Wurtz et al. (2006) paper
#                          
################################################################################

.filterGarch11Fit <- function(data, parm)
{
  mu = parm[1]; omega = parm[2]; alpha = parm[3]; beta = parm[4]
  z = (data-mu); Mean = mean(z^2)
  
  # Use Filter Representation:
  e = omega + alpha * c(Mean, z[-length(data)]^2)
  h = filter(e, beta, "r", init = Mean)
  print(h[1])
  hh = sqrt(abs(h))
  cbind(z,hh)  
}

################################################################################













################################################################################
#  FUNCTION:               DESCRIPTION:
#
#  .filterAparchForLoop   Conditional Variance filtering - 
#                          for loop - Wuertz et al. (2006) 
#                          
################################################################################

.filterAparchForLoop <- function(
    data,h.init = 0.1,
    p,q, 
    mu, omega, alpha, beta, gamma, delta)
{
  
    # Return
    # the series 'z' and 'hh'
    
    
    # input 
    # data: vector with data
    
    # error treatment of input parameters
    if( p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
        || length(beta) != q)
      stop("One or more of these conditions were true:
            p < 1 || q < 1 || (length(alpha) != length(gamma)) || length(alpha) != p 
            || length(beta) != q")
  
    # Initial declaration of variables
    pq = max(p,q)
    
    z = (data-mu)
    N = length(data)
    Mean.z = mean(abs(z)^delta)
    h = rep(h.init, pq)
    
    # Calculate h[(pq+1):N] recursively
    for (i in (pq+1):N )
    {
        ed = 0
        for (j in 1:p)
        {        
            ed = ed+alpha[j]*(abs(z[i-j])-gamma[j]*z[i-j])^delta
        }
        h[i] = omega + ed + sum(beta*h[i-(1:q)])
    }
    if(length(z) != length(h))
      stop("Error in filtering function. length(z) != length(h)")  
    
    hh <- (abs(h))^(1/delta)
    
    # return
    cbind(z,hh)    
}

################################################################################







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
#  .getFormula             Gets the formula.mean and formula.variance from
#                          the object formula.
################################################################################


.getFormula <-
  function(
    formula)
  {
    # Description:
    #   This functions reads formula object and converts it into a list 
    #   containing the separeted mean and variance arguments.
    #   Examples from output:
    #     ~arma(1,1)+garch(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(1,1)
    #     ~aparch(1,1): formula.mean = ~arma(0,0); formula.variance = ~aparch(1,1)
    #     ~arch(1): formula.mean = ~arma(0,0); formula.variance = ~arch(1)
    #     ~arma(1,1): formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
    #     ~ar(1): formula.mean = ~ar(1); formula.variance = ~garch(0,0) 
    #     ~ma(1): formula.mean = ~ma(1); formula.variance = ~garch(0,0) 
    
    # Arguments:
    #   formula - ARMA(m,n) + GARCH/APARCH(p,q) mean and variance specification 
    
    # Return:
    #   A list containing two elements, formula.mean and formula.variance, 
    #   formula.order and the boolean formula.isAPARCH
    
    # FUNCTION: 
    
    # Initial variable declaration
    allLabels = attr(terms(formula), "term.labels")
    formulas.mean.allowed = c("arma")
    formulas.variance.allowed = c("garch","aparch")
    formulaOK <- TRUE
    checkFormulaMean <- ""
    checkFormulaVariance <- ""
    isAPARCH = FALSE
    
    # Error treatment of input parameters 
    if( (length(allLabels) != 1 ) && (length(allLabels) != 2 ) )
      formulaOK <- FALSE
    
    # Formula of type: ~ formula1 + formula2
    else if (length(allLabels) == 2)
    {
      formula.mean = as.formula(paste("~", allLabels[1]))
      formula.var = as.formula(paste("~", allLabels[2]))
      checkFormulaMean = rev(all.names(formula.mean))[1]
      checkFormulaVariance = rev(all.names(formula.var))[1]
      if( !any(formulas.mean.allowed == checkFormulaMean) || 
            !any(formulas.variance.allowed == checkFormulaVariance))   
        formulaOK <- FALSE
    }
    # Formula of type: ~formula1
    else if (length(allLabels) == 1) 
    {
      
      # pure 'garch' or 'aparch'
      if(grepl("arch", attr(terms(formula), "term.labels"))) 
      {
        formula.mean = as.formula("~ arma(0, 0)")
        formula.var = as.formula(paste("~", allLabels[1]))
        checkFormulaVariance = rev(all.names(formula.var))[1]
        if(!any(formulas.variance.allowed == checkFormulaVariance))   
          formulaOK <- FALSE
      }
      else # pure 'ar', 'ma' or 'arma' model.
      {
        formula.mean = as.formula(paste("~", allLabels[1]))  
        formula.var = as.formula("~ garch(0, 0)") 
        checkFormulaMean = rev(all.names(formula.mean))[1]
        if(!any(formulas.mean.allowed == checkFormulaMean))   
          formulaOK <- FALSE
      }    
    }
    
    # Check if we are fitting "aparch" model 
    if(checkFormulaVariance == "aparch")
      isAPARCH = TRUE
    
    # Get model order and check if they were specified correctly
    if(formulaOK == TRUE)
    {
      model.order.mean = 
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.mean), 
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      model.order.var = 
        as.numeric(strsplit(strsplit(strsplit(as.character(formula.var), 
                                              "\\(")[[2]][2], "\\)")[[1]], ",")[[1]])
      if( (length(model.order.mean) != 2) || (length(model.order.mean) != 2))
        formulaOK <- FALSE 
    }
    
    # Check if model order was specified correctly.
    if(formulaOK == TRUE)
    {
      m = model.order.mean[1]
      n = model.order.mean[2]
      p = model.order.var[1]
      q = model.order.var[2]        
      if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
           any (c(m,n,p,q) < 0) || (p == 0 && q != 0))
        formulaOK <- FALSE            
    }
    
    
    # Stop if formula was not specified correctly
    if(formulaOK == FALSE)
      stop ("Invalid Formula especification. 
            Formula mean must be 'arma' and 
            Formula Variance must be one of: garch or aparch
            For example:
            ARMA(1,1)-GARCH(1,1):  ~arma(1,1)+garch(1,1),
            AR(1)-GARCH(1,1):      ~arma(1,0)+garch(1,1),
            MA(1)-APARCH(1,0):     ~arma(0,1)+aparch(1,0),
            ARMA(1,1):             ~arma(1,1),
            ARCH(2):               ~garch(1,0),
            For more details just type: ?gsFit")
    
    # Return
    list(formula.mean = formula.mean,formula.var = formula.var, 
         formula.order = c(m,n,p,q), isAPARCH = isAPARCH)
  }

################################################################################











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
#  .getStart               Returns initial and boundary values to 
#                          perform optimization 
################################################################################


.getStart <- function(data,m,n,p,q, AR = FALSE, MA = FALSE,
                      cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
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
    #   Result - A tree columns matrix with columns representing the start, lower and upper
    #     bounds of the parameter set
    
    # FUNCTION:  
    
    # Error control of input parameters
    if (m < 0 || n < 0 || m %% 1 != 0 || n %% 1 != 0)
        stop("'m' and 'n' need to be integers greater than zero")
    
    if (m == 0 || n == 0)
        stop("Expects 'm' and 'n' different from zero")
    
    if ( (m != 1 && AR == TRUE)  || (n != 1 && MA == TRUE) )
        stop("If AR = TRUE, 'm' should be 1 and if MA = TRUE 'n' should be 1")
    
    if (p < 0 || q < 0 || p %% 1 != 0 || q %% 1 != 0)
        stop("'p' and 'q' need to be integers greater than zero")
    
    if(p == 0 && q != 0)
        stop("Invalid Garch(p,q) order")
    
    if( !is.numeric(data) || !is.vector(data))
      stop("data set must be a numerical one dimensional vector")
    
    # Initial variable declaration
    cond.dist = match.arg(cond.dist)
    cond.dist.list = c("stable", "gev", "gat", "norm", "std", "sstd", "skstd", "ged")
    Mean <- mean(data)
    Var <- var(data)
    Dispersion <- mean(abs(data-Mean))
    arima.fit <- c()
    arima.fit.try <- "empty"  
    arima.m <- m
    arima.n <- n
    if(AR == TRUE) # we don't have the AR part
        arima.m <- 0
    if(MA == TRUE) # we don't have the MA part 
        arima.n <- 0
    
    # Try arima fit function to get initial arma parameters
    try(arima.fit.try <- as.vector(arima(data,order = c(arima.m, 0, arima.n))$coef), silent = TRUE)
    if( is.numeric(arima.fit.try) )
    {   
        arima.fit <- arima.fit.try
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
    } else {   
        mean.init <- Mean
        ar.init <- rep(0,m)
        ma.init <- rep(0,n)
    }
 
    # START VALUES
    
    arma.start = c(ar.init, ma.init)
    mu.start = mean.init
    gm.start = rep(0, p)
    
    omega.start = list(
      "stableS0" = 0.1 * Dispersion,
      "stableS1" = 0.1 * Dispersion,
      "stableS2" = 0.1 * Dispersion,
      "gev" = 0.1 * Var,
      "gat" = 0.1 * Var,
      "norm" = 0.1 * Var,
      "std" = 0.1 * Var,
      "sstd" = 0.1 * Var,
      "skstd" = 0.1 * Var,
      "ged" = 0.1 * Var)
    
    alpha.start = list(
      "stableS0" = rep(0.1/p, p),
      "stableS1" = rep(0.1/p, p),
      "stableS2" = rep(0.1/p, p),
      "gev" = rep(0.05/p, p),
      "gat" = rep(0.1/p, p),
      "norm" = rep(0.1/p, p),
      "std" = rep(0.1/p, p),
      "sstd" = rep(0.1/p, p),
      "skstd" = rep(0.1/p, p),
      "ged" = rep(0.1/p, p))
    
    beta.start = list(
      "stableS0" = rep(0.8/q, q),
      "stableS1" = rep(0.8/q, q),
      "stableS2" = rep(0.8/q, q),
      "gev" = rep(0.7/q, q),
      "gat" = rep(0.8/q, q),
      "norm" = rep(0.8/q, q),
      "std" = rep(0.8/q, q),
      "sstd" = rep(0.8/q, q),
      "skstd" = rep(0.8/q, q),
      "ged" = rep(0.8/q, q))
    
    delta.start = list(
      "stableS0" = 1.05,
      "stableS1" = 1.05,
      "stableS2" = 1.05,
      "gev" = 2,
      "gat" = 2,
      "norm" = 2,
      "std" = 2,
      "sstd" = 2,
      "skstd" = 2,
      "ged" = 2)
    
    skew.start = list(
      "stableS0" = 0,
      "stableS1" = 0,
      "stableS2" = 0,
      "gev" = 1,
      "gat" = 1,
      "norm" = 1,
      "std" = 1,
      "sstd" = 1,
      "skstd" = 1,
      "ged" = 1)
    
    shape.start = list(
      "stableS0" = 1.9,
      "stableS1" = 1.9,
      "stableS2" = 1.9,
      "gev" = 0, # numerical tests showed that 0.01 is a good starting parameter.
      "gat" = c(2, 4),
      "norm" = 1,
      "std" = 4,
      "sstd" = 4,
      "skstd" = 4,
      "ged" = 4)
        
    # LOWER BOUNDS
    
    mu.lower = min(- 100 * abs ( mu.start ), -100)
    arma.lower = rep ( - 100, m + n )
    omega.lower = TOLG
    alpha.lower = rep ( TOLG, p )
    beta.lower = rep ( TOLG, q )
    gm.lower = rep( - 1 + TOLG, p)
    
    delta.lower = list(
      "stableS0" = 1,
      "stableS1" = 1,
      "stableS2" = 1,
      "gev" = TOLG,
      "gat" = TOLG,
      "norm" = TOLG,
      "std" = TOLG,
      "sstd" = TOLG,
      "skstd" = TOLG,
      "ged" = TOLG)
    
    skew.lower = list(
      "stableS0" = - 1 + TOLSTABLE,
      "stableS1" = - 1 + TOLSTABLE,
      "stableS2" = - 1 + TOLSTABLE,
      "gev" = 0,
      "gat" = TOLG,
      "norm" = 0,
      "std" = 0,
      "sstd" = TOLG,
      "skstd" = TOLG,
      "ged" = 0)
    
    shape.lower = list(
      "stableS0" = 1 + TOLSTABLE,
      "stableS1" = 1 + TOLSTABLE,
      "stableS2" = 1 + TOLSTABLE,
      "gev" = - 0.5 + TOLG, # to ensure good MLE properties. See Jondeau et al. 
      "gat" = c ( TOLG, TOLG),
      "norm" = 0,
      "std" = 2 + TOLG,
      "sstd" = 2 + TOLG,
      "skstd" = 2 + TOLG, # to ensure finiteness of variance
      "ged" = TOLG)
       
    # UPPER BOUNDS
    
    mu.upper = max( 100 * abs ( mu.start ), 100)
    arma.upper = rep ( 100, m + n )
    omega.upper = 100 * abs ( Dispersion )
    alpha.upper = rep ( 1 - TOLG, p )
    beta.upper = rep ( 1 - TOLG, q )
    gm.upper = rep (1 - TOLG, p)

    delta.upper = list(
      "stableS0" = 2 - TOLSTABLE,
      "stableS1" = 2 - TOLSTABLE,
      "stableS2" = 2 - TOLSTABLE,
      "gev" = 100,
      "gat" = 100,
      "norm" = 100,
      "std" = 100,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 100)
    
    skew.upper = list(
      "stableS0" = 1 - TOLSTABLE,
      "stableS1" = 1 - TOLSTABLE,
      "stableS2" = 1 - TOLSTABLE,
      "gev" = 2,
      "gat" = 100,
      "norm" = 2,
      "std" = 2,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 2)
    
    shape.upper = list(
      "stableS0" = 2 - TOLSTABLE,
      "stableS1" = 2 - TOLSTABLE,
      "stableS2" = 2 - TOLSTABLE,
      "gev" = 0.5 - TOLG, # to ensure finiteness of the variance and mean. 
      "gat" = c ( 100, 100),
      "norm" = 2,
      "std" = 100,
      "sstd" = 100,
      "skstd" = 100,
      "ged" = 100)    

    # CHECK IF THE STARTING MODELS IS STATIONARY ( ONLY FOR THE sqp.restriction ALGORITHM )
    if( ! any ( cond.dist == c("stableS0", "stableS2") ) ){
      
        start.model.persistency = 1
        start.model.persistency = sum( alpha.start[[cond.dist]] * gsMomentAparch( 
              cond.dist = cond.dist, shape = shape.start[[cond.dist]], 
              skew = skew.start[[cond.dist]], delta = delta.start[[cond.dist]], gm = 0) ) + 
              sum ( beta.start[[cond.dist]] )
        # The real condition is 1 but we let the user to set the limit between 0.95 and 0.999...
        if(start.model.persistency >= 0.95)
        {
            print(start.model.persistency)
            print(paste("The starting model with conditional",cond.dist, " is not stationary."))
            stop("Change the starting value of the parameters for the reported model.")
        }
    }
    
    # CONSTRUCT THE RESULT
    
    start = c ( mu.start, arma.start, omega.start[[cond.dist]], alpha.start[[cond.dist]],
               gm.start, beta.start[[cond.dist]], delta.start[[cond.dist]], 
               skew.start[[cond.dist]], shape.start[[cond.dist]])
  
    lower = c ( mu.lower, arma.lower, omega.lower, alpha.lower,
              gm.lower, beta.lower, delta.lower[[cond.dist]], 
              skew.lower[[cond.dist]], shape.lower[[cond.dist]])
  
    upper = c ( mu.upper, arma.upper, omega.upper, alpha.upper,
              gm.upper, beta.upper, delta.upper[[cond.dist]], 
              skew.upper[[cond.dist]], shape.upper[[cond.dist]])
    
    namesStart = c("mu", paste("ar", 1:m, sep = ""), 
                   paste("ma", 1:n, sep = ""),
                   "omega", paste("alpha", 1:p, sep = ""), 
                   paste("gm", 1:p, sep = ""), 
                   paste("beta", 1:q, sep = ""), "delta","skew",
                   paste("shape", 1:length(shape.start[[cond.dist]]), sep = ""))        
    
    # Create result
    result = rbind(start,lower,upper)
    colnames(result) = namesStart
    
    # Return
    result
}



# ------------------------------------------------------------------------------






################################################################################


# Manually load files from GEVStableGarch changing 
#.GSgarch.dstable for the appropriate functions from libstable4u in 
# https://cran.r-project.org/web/packages/libstable4u/index.html

# I found that these are the same
# pars = c(alpha, beta, gamma, delta)
# stable_pdf(x, pars) AND  stabledist::dstable(x, alpha, beta, gamma, delta, pm = 0)



.GSgarch.dstable = function (x, alpha = 1.5, beta = 0, gamma = 1, delta = 0, param = 0) 
{
  # return(stabledist::dstable(x, alpha, beta, gamma, delta, pm = param))
  return(stable_pdf(x, c(alpha, beta, gamma, delta)))
}



.armaGarchDist = function (z, hh, shape = 1.5, skew = 0, cond.dist = c("stableS0", 
                                                      "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                      "skstd", "ged"), TOLG = 1e-08) 
{
  cond.dist = match.arg(cond.dist)
  if (length(z) != length(hh)) 
    stop("Error: Vectors 'z' and 'hh' have different length.")
  if (sum(is.na(hh)) > 0 || min(hh) == 0) {
    return(1e+99)
  }
  if (cond.dist == "norm") 
    return(-sum(log(dnorm(x = z/hh)/hh)))
  if (cond.dist == "std") {
    if (!(shape > 2)) 
      return(Inf)
    return(-sum(log(dstd(x = z/hh, nu = shape)/hh)))
  }
  if (cond.dist == "sstd") {
    if (!(shape > 2) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dsstd(x = z/hh, nu = shape, xi = skew)/hh)))
  }
  if (cond.dist == "skstd") {
    if (!(shape > 2) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dskstd(x = z/hh, nu = shape, xi = skew)/hh)))
  }
  if (cond.dist == "gat") {
    if (!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dgat(x = z/hh, nu = shape[1], d = shape[2], 
                         xi = skew)/hh)))
  }
  if (cond.dist == "ged") {
    if (!(shape > 0)) 
      return(1e+99)
    return(-sum(log(dged(x = z/hh, nu = shape)/hh)))
  }
  if (cond.dist == "gev") {
    if ((shape[1] <= -0.5) || (shape[1] >= 0.5)) 
      return(1e+99)
    if (abs(shape[1]) < 1e-07) {
      llh <- sum(log(hh)) + sum(z/hh) + sum(exp(-z/hh))
      return(llh)
    }
    sig <- hh
    y <- 1 + shape * (z/hh)
    llh <- sum(log(sig)) + sum(y^(-1/shape)) + sum(log(y)) * 
      (1/shape + 1)
    return(llh)
  }
  if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) {
    if (getOption(".stableIsLoaded", default = FALSE) == 
        TRUE) {
      .GSgarch.dstable <- function(x, alpha = 1.5, beta = 0, 
                                   gamma = 1, delta = 0, param = 1) {
        return(stable::dstable.quick(x, alpha, beta, 
                                     gamma, delta, param))
      }
    }
    if (!(shape > 1) || !(shape < 2) || !(abs(skew) < 1)) 
      return(1e+99)
    y <- z/hh
    if (sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y))) 
      return(1e+99)
    if (cond.dist == "stableS0") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 0)/hh))
    if (cond.dist == "stableS1") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 1)/hh))
    if (cond.dist == "stableS2") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 2)/hh))
    return(result)
  }
}



gsFit = function (formula = ~garch(1, 1), data, cond.dist = c("stableS0", 
                                                      "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                      "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
                                                                                                          "sqp.restriction", "nlminb", "nlminb+nm"), control = NULL, 
          tolerance = NULL, title = NULL, description = NULL) 
{
  DEBUG = FALSE
  cond.dist = match.arg(cond.dist)
  algorithm = match.arg(algorithm)
  if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || 
      any(!is.finite(data))) 
    stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")
  CALL = match.call()
  if (is.null(tolerance)) 
    tolerance = list(TOLG = 1e-08, TOLSTABLE = 0.01, TOLSTATIONARITY = 0.001)
  if (tolerance$TOLSTATIONARITY > 0.05) 
    stop("TOLSTATIONARITY can not be > 0.05")
  if (is.null(title)) 
    title = "ARMA-GARCH modelling"
  formula.input <- formula
  formula <- .getFormula(formula)
  m <- formula$formula.order[1]
  n <- formula$formula.order[2]
  p <- formula$formula.order[3]
  q <- formula$formula.order[4]
  APARCH <- formula$isAPARCH
  formula.mean <- ""
  formula.var <- ""
  if (m > 0 || n > 0) 
    formula.mean <- formula$formula.mean
  else formula.mean <- ""
  if (p > 0) 
    formula.var <- formula$formula.var
  else formula.var <- ""
  if (algorithm == "sqp.restriction") {
    if (formula.var == "") 
      stop("sqp.restriction should only be used with GARCH/APARCH models.")
    if (any(cond.dist == c("stableS0", "stableS2"))) 
      stop("sqp.restriction algorithm can only be used with \n stable distribution in S1 parametrization, i.e., cond.dist = 'stableS1'.")
  }
  if (formula.var == "") 
    stop("Pure ARMA model not allowed")
  printRes = TRUE
  ARMAonly <- FALSE
  AR <- FALSE
  MA <- FALSE
  GARCH <- FALSE
  if (m == 0) 
    AR <- TRUE
  if (n == 0) 
    MA <- TRUE
  if (q == 0) 
    GARCH <- TRUE
  optim.finished <- FALSE
  if (AR == TRUE) 
    m <- 1
  if (MA == TRUE) 
    n <- 1
  if ((p == 0) && (q == 0)) 
    ARMAonly = TRUE
  if (GARCH == TRUE && !ARMAonly) 
    q <- 1
  data <- data
  N <- length(data)
  out <- NULL
  messages <- NULL
  out$order <- c(m, n, p, q, include.mean, APARCH)
  TMPvector <- c(if (m != 0 || n != 0) c("arma(", m, ",", n, 
                                         ")-"), if (APARCH == FALSE) c("garch(", p, ",", q, ")"), 
                 if (APARCH == TRUE) c("aparch(", p, ",", q, ")"))
  TMPorder <- paste(TMPvector, collapse = "")
  TMPvectorintercept <- c("include.mean:", if (include.mean == 
                                               TRUE) "TRUE", if (include.mean == FALSE) "FALSE")
  TMPintercept <- paste(TMPvectorintercept, collapse = "")
  out$model <- paste(TMPorder, "##", TMPintercept, collapse = "")
  out$cond.dist <- cond.dist
  out$data <- data
  optim.finished <- FALSE
  if (cond.dist == "gat") 
    lengthShape = 2
  else lengthShape = 1
  if (DEBUG) 
    print(c("lengthShape", lengthShape))
  arCheck <- function(ar) {
    p <- max(which(c(1, -ar) != 0)) - 1
    if (!p) 
      return(TRUE)
    all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
  }
  garchLLH = function(parm) {
    if (sum(is.nan(parm)) != 0) {
      return(1e+99)
    }
    mu <- parm[1]
    a <- parm[(1 + 1):(2 + m - 1)]
    b <- parm[(1 + m + 1):(2 + m + n - 1)]
    omega <- parm[1 + m + n + 1]
    alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
    gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
    beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
                                            p + q - 1)]
    delta <- parm[2 + m + n + 2 * p + q + 1]
    skew <- parm[3 + m + n + 2 * p + q + 1]
    shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
                                                 2 * p + q + lengthShape)]
    if (!APARCH) {
      gm = rep(0, p)
      delta = 2
      if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) 
        delta = 1
    }
    if (AR == TRUE) 
      a <- 0
    if (MA == TRUE) 
      b <- 0
    if (GARCH == TRUE) 
      beta <- 0
    if (include.mean == FALSE) 
      mu <- 0
    cond.general <- FALSE
    cond.normal <- FALSE
    cond.student <- FALSE
    cond.gev <- FALSE
    cond.stable <- FALSE
    parset <- c(omega, alpha, if (!GARCH) beta, delta)
    cond.general <- any(parset < tolerance$TOLG)
    if (cond.general || cond.student || cond.gev || cond.stable) {
      return(1e+99)
    }
    if (!arCheck(a)) 
      return(1e+99)
    if (AR == TRUE && MA == TRUE) {
      filteredSeries <- .filterAparch(data = data, p = p, 
                                      q = q, mu = mu, omega = omega, alpha = alpha, 
                                      beta = beta, gamma = gm, delta = delta)
      z <- filteredSeries[, 1]
      hh <- filteredSeries[, 2]
    }
    if (AR == FALSE || MA == FALSE) {
      filteredArma <- .filterArma(data = data, size.data = N, 
                                  m = m, n = n, mu = mu, a = a, b = b)
      filteredSeries <- .filterAparch(data = filteredArma, 
                                      p = p, q = q, mu = 0, omega = omega, alpha = alpha, 
                                      beta = beta, gamma = gm, delta = delta)
      z <- filteredSeries[, 1]
      hh <- filteredSeries[, 2]
    }
    if (optim.finished) {
      out$residuals <<- as.numeric(z)
      out$sigma.t <<- as.numeric(hh)
      out$h.t <<- as.numeric(hh^delta)
    }
    llh.dens <- .armaGarchDist(z = z, hh = hh, shape = shape, 
                               skew = skew, cond.dist = cond.dist)
    llh <- llh.dens
    if (is.nan(llh) || is.infinite(llh) || is.na(llh)) {
      llh <- 1e+99
    }
    llh
  }
  start <- .getStart(data = data, m = m, n = n, p = p, q = q, 
                     AR = AR, MA = MA, cond.dist = cond.dist, TOLG = tolerance$TOLG, 
                     TOLSTABLE = tolerance$TOLSTABLE)
  if (DEBUG) {
    print("start")
    print(start)
  }
  garch.stationarity <- function(parm) {
    omega <- parm[1 + m + n + 1]
    alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
    gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
    beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
                                            p + q - 1)]
    delta <- parm[2 + m + n + 2 * p + q + 1]
    skew <- parm[3 + m + n + 2 * p + q + 1]
    shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
                                                 2 * p + q + lengthShape)]
    model = list(omega = omega, alpha = alpha, gm = gm, beta = beta, 
                 delta = delta, skew = skew, shape = shape)
    .stationarityAparch(model = model, formula = formula, 
                        cond.dist = cond.dist)
  }
  modelLLH <- garchLLH
  if (algorithm == "nlminb") {
    fit1 <- nlminb(start[1, ], objective = modelLLH, lower = start[2, 
    ], upper = start[3, ], control = control)
    out$llh <- fit1$objective
    out$par <- fit1$par
    out$hessian <- fit1$hessian
    out$convergence <- fit1$convergence
  }
  if (algorithm == "nlminb+nm") {
    fit1.partial <- nlminb(start[1, ], objective = modelLLH, 
                           lower = start[2, ], upper = start[3, ], control = control)
    fit1 <- optim(par = fit1.partial$par, fn = modelLLH, 
                  method = "Nelder-Mead", hessian = TRUE)
    out$llh <- fit1$value
    out$par <- fit1$par
    out$hessian <- fit1$hessian
    out$convergence <- fit1$convergence
  }
  if (algorithm == "sqp") {
    fit1 <- solnp(pars = start[1, ], fun = modelLLH, LB = start[2, 
    ], UB = start[3, ], control = control)
    out$llh <- fit1$values[length(fit1$values)]
    out$par <- fit1$pars
    out$hessian <- fit1$hessian
    out$convergence <- fit1$convergence
  }
  if (algorithm == "sqp.restriction") {
    fit1 <- solnp(pars = start[1, ], fun = modelLLH, ineqfun = garch.stationarity, 
                  ineqLB = 0, ineqUB = 1 - tolerance$TOLSTATIONARITY, 
                  LB = start[2, ], UB = start[3, ], control = control)
    out$llh <- fit1$values[length(fit1$values)]
    out$par <- fit1$pars
    out$hessian <- fit1$hessian
    sizeHessian = length(out$hessian[1, ]) - 1
    out$hessian = out$hessian[1:sizeHessian, 1:sizeHessian]
    out$convergence <- fit1$convergence
  }
  if ((any(c("sqp", "nlminb") == algorithm)) && !is.numeric(try(sqrt(diag(solve(out$hessian))), 
                                                                silent = TRUE))) {
    fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                          method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
    if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
      out$hessian = fit1.partial$hessian
      if (DEBUG) 
        print(abs(fit1.partial$par - fit1$par))
    }
  }
  if ((algorithm == "sqp.restriction")) {
    fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                          method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
    if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
      out$hessian = fit1.partial$hessian
      if (DEBUG) 
        print(abs(fit1.partial$par - fit1$par))
    }
  }
  if (DEBUG) 
    print(fit1)
  if (out$llh == 1e+99) 
    out$convergence = 1
  optim.finished = TRUE
  modelLLH(out$par)
  if (!ARMAonly) {
    outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
                                                           1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
                                                                                                           m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
                                                                                                                                                                         m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                        m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
                                                                                                                                                                                                                                                                               m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                                                                                                                                              m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
                                                                                                                                                                                                                                                                                                                                                                                "stableS0", "stableS1", "stableS2", "gat") == cond.dist)) (3 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                             m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "stableS1", "stableS2", "sstd", "skstd", "ged", "gat") == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             cond.dist)) (4 + m + n + 2 * p + q + 1):(4 + m + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n + 2 * p + q + lengthShape))
  }
  else {
    outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
                                                           1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
                                                                                                           m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
                                                                                                                                                                         m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                        m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
                                                                                                                                                                                                                                                                               m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                                                                                                                                              m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
                                                                                                                                                                                                                                                                                                                                                                                "stableS0", "stableS1", "stableS2", "gat") == cond.dist)) (1 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                             m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "stableS1", "stableS2", "sstd", "skstd", "ged", "gat") == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             cond.dist)) (2 + m + n + 2 * p + q + 1):(2 + m + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n + 2 * p + q + lengthShape), length(out$par))
  }
  outnames <- c(if (include.mean) "mu", if (AR == FALSE) paste("ar", 
                                                               1:m, sep = ""), if (MA == FALSE) paste("ma", 1:n, sep = ""), 
                if (!ARMAonly) "omega", if (!ARMAonly) paste("alpha", 
                                                             1:p, sep = ""), if (APARCH) paste("gamma", 1:p, sep = ""), 
                if (!GARCH) paste("beta", 1:q, sep = ""), if (APARCH) "delta", 
                if (any(c("sstd", "stableS0", "stableS1", "stableS2", 
                          "gat", "skstd") == cond.dist)) "skew", if (any(c("std", 
                                                                           "gev", "stableS0", "stableS1", "stableS2", "sstd", 
                                                                           "ged", "gat", "skstd") == cond.dist)) paste("shape", 
                                                                                                                       1:lengthShape, sep = ""), if (ARMAonly) "sigma")
  if (DEBUG) {
    print(c("out", out))
    print(c("outindex", outindex))
  }
  out$par <- out$par[outindex]
  names(out$par) <- outnames
  out$hessian <- out$hessian[outindex, outindex]
  nParam <- length(out$par)
  out$aic = 2 * out$llh + 2 * nParam
  out$aicc = 2 * out$llh + 2 * (nParam + 1) * N/(N - nParam - 
                                                   2)
  out$bic = 2 * out$llh + nParam * log(N)
  out$ics = c(out$aic, out$bic, out$aicc)
  names(out$ics) <- c("AIC", "BIC", "AICc")
  if (printRes) {
    solveHessianFailed = FALSE
    out$se.coef <- 0
    out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
    if (!is.numeric(out$se.coef)) {
      solveHessianFailed = TRUE
      messages$solving.hessian.matrix = "Error inverting Hessian Matrix"
      out$matcoef <- cbind(out$par, rep(NA, length(out$par)), 
                           rep(NA, length(out$par)), rep(NA, length(out$par)))
      dimnames(out$matcoef) = dimnames(out$matcoef) = list(names(out$par), 
                                                           c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
    }
    else {
      out$tval <- try(out$par/out$se.coef, silent = TRUE)
      out$matcoef = cbind(out$par, if (is.numeric(out$se.coef)) 
        out$se.coef, if (is.numeric(out$tval)) 
          out$tval, if (is.numeric(out$tval)) 
            2 * (1 - pnorm(abs(out$tval))))
      dimnames(out$matcoef) = list(names(out$tval), c(" Estimate", 
                                                      " Std. Error", " t value", "Pr(>|t|)"))
    }
    cat("\nFinal Estimate of the Negative LLH:\n")
    cat("-LLH:", out$llh)
    if (out$convergence == 0) 
      messages$optimization.algorithm = "Algorithm achieved convergence"
    else messages$optimization.algorithm = "Algorithm did not achieved convergence"
    cat("\nCoefficient(s):\n")
    printCoefmat(round(out$matcoef, digits = 6), digits = 6, 
                 signif.stars = TRUE)
  }
  out$order <- c(formula$formula.order[1], formula$formula.order[2], 
                 formula$formula.order[3], formula$formula.order[4])
  names(out$order) <- c("m", "n", "p", "q")
  fit <- list(par = out$par, llh = out$llh, hessian = out$hessian, 
              ics = out$ics, order = out$order, cond.dist = cond.dist, 
              se.coef = out$se.coef, tval = out$tval, matcoef = out$matcoef)
  new("GEVSTABLEGARCH", call = as.call(match.call()), formula = formula.input, 
      method = "Max Log-Likelihood Estimation", convergence = out$convergence, 
      messages = messages, data = data, fit = fit, residuals = out$residuals, 
      h.t = out$h.t, sigma.t = as.vector(out$sigma.t), title = as.character(title), 
      description = as.character(description))
}


gsSelect = function (data, order.max = c(1, 1, 1, 1), selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, cond.dist = c("stableS0", 
                                                                                                                   "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                                                                                   "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
                                                                                                                                                                       "sqp.restriction", "nlminb", "nlminb+nm"), ...) 
{
  cond.dist = match.arg(cond.dist)
  algorithm = match.arg(algorithm)
  selection.criteria = match.arg(selection.criteria)
  if (!is.numeric(data) || !is.vector(data)) 
    stop("data set must be a numerical one dimensional vector")
  goodness.of.fit.min <- 1e+99
  goodness.of.fit.current <- 1e+99
  fit.min <- list()
  fit <- list()
  order.list <- .getOrder(order.max = order.max)
  order.list.size <- length(order.list[, 1])
  for (i in 1:order.list.size) {
    m = order.list[i, 1]
    n = order.list[i, 2]
    p = order.list[i, 3]
    q = order.list[i, 4]
    if (is.aparch == TRUE) {
      formula = as.formula(paste("~ arma(", m, ", ", n, 
                                 ") + aparch(", p, ", ", q, ")", sep = "", collapse = NULL))
    }
    else {
      formula = as.formula(paste("~ arma(", m, ", ", n, 
                                 ") + garch(", p, ", ", q, ")", sep = "", collapse = NULL))
    }
    cat("\n------------------------------------------------------------------------------------------\n")
    cat(paste("Model: ", paste(formula)[2], " with '", cond.dist, 
              "' conditional distribution", sep = ""))
    cat("\n------------------------------------------------------------------------------------------\n")
    fit = gsFit(data = data, formula = formula, cond.dist = cond.dist, 
                algorithm = algorithm, include.mean = include.mean, 
                ...)
    goodness.of.fit.current = fit@fit$ics[selection.criteria]
    if (goodness.of.fit.current < goodness.of.fit.min) {
      fit.min <- fit
      goodness.of.fit.min <- fit@fit$ics[selection.criteria]
    }
  }
  cat("\n------------------------------------------------------------------------------------------\n")
  cat(paste("Best Model: ", paste(fit.min@formula)[2]))
  cat("\n------------------------------------------------------------------------------------------\n")
  return(fit.min)
}



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
#  .getOrder        Return a matrix with parameter order to be used 
#                          inside function GSgarch.FitAIC 
################################################################################


.getOrder <- 
    function(order.max = c(1,1,1,1))
{
      
    # Description:
    #   Iterates over the parameters to create a vector with parameter
    #   orders (like (1,1,0,1)) to use inside function GSgarch.FitAIC
    
    # Arguments:
    #   nMAX, mMAX, pMAX, qMAX - maximum order to be estimated
    
    # Return:
    #   arma.garch.order - A matrix with several model orders in the 
    #   format [m,n,p,q] 
    
    # FUNCTION:            
      
    # error treatment on input parameters
    m = order.max[1]; n = order.max[2]; p = order.max[3]; q = order.max[4] 
    if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
        any (c(m,n,p,q) < 0) || (p == 0 && q != 0) ||
        any (c(m,n,p,q) > 10) ) 
        stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
              The order of the parameters could be set up to 10.")
    arma.garch.order <- c()
    for(i1 in 0:m)
    {
        for(i2 in 0:n)
        {
            for(i3 in 1:p)
            {
                for(i4 in 0:q)
                { 
                    ord <- c(i1,i2,i3,i4)
                    arma.garch.order <- rbind(arma.garch.order,c(i1,i2,i3,i4))
                }
            }
        }
    }
    return(arma.garch.order)
}


################################################################################




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


# Copyrights (C) 2015, Thiago do Rego Sousa <thiagoestatistico@gmail.com>
# This is a modified version of the code contained inside file 
# garch-Spec.R from package fGarch, version 3010.82.

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               SPECIFICATION:
#  gsSpec               Creates a 'GEVSTABLEGARCHSPEC' object from scratch
################################################################################


gsSpec <-
    function (model = list(), presample = NULL,
              cond.dist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
              rseed = NULL)
{
      
    # A function originally implemented by Diethelm Wuertz and modified
    # to be used inside package GEVStableGarch. See the latest copyright notice.       
      
    # Description:
    #   Creates a "gsSpec" object.
    
    # Arguments:
    #   model - a list with the model parameters as entries
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the mean value for ARMA specification,
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     delta - the exponent value used in the variance equation.
    #     skew - a numeric value listing the distributional
    #        skewness parameter.
    #     shape - a numeric value listing the distributional
    #        shape parameter.
    #   presample - a numeric "matrix" with 3 columns and
    #       at least max(m,n,p,q) rows. 
    #       The first culumn are the innovations, the second
    #       the conditional variances, and the last the time series.
    #       When presample is missing, we construct our presample matrix
    #       as [z,h,y] where z = rnorm(0,1), h = "uev" recursion 
    #       initialization described in Wuertz et al. (2006) and
    #       y = mu. Note that the conditional variance column
    #       can contain only strictly positive numbers.
    #       If the model is a pure AR or MA the presample 
    #       matrix will become a simple array.
    #   cond.dist - a character string naming the distribution
    #       function.
    #   rseed - optional random seed. The default seed is '0'.
     
    # Return: An object of class GEVSTABLEGARCHSPEC   
        # Slots:
        #   call - the function call.
        #   formula - a formula object describing the model, e.g.
        #       ARMA(m,n) + GARCH(p,q). ARMA can be missing or
        #       specified as AR(m) or MA(n) in the case of pure
        #       autoregressive or moving average models. GARCH may
        #       alternatively be specified as ARCH(p) or APARCH(p,q).
        #   model - as declared in the input.

    # FUNCTION
      
    # Skewness Parameter Settings:
    skew = list(
        "stableS0" = 0,
        "stableS1" = 0,
        "stableS2" = 0,
        "gev" = NULL,
        "gat" = 1,
        "norm" = NULL,
        "std" = NULL,
        "sstd" = 0.9,
        "skstd" = 1,
        "ged" = NULL)

    # Shape Parameter Settings:
    shape = list(
        "stableS0" = 1.7,
        "stableS1" = 1.7,
        "stableS2" = 1.7,
        "gev" = 0.3,
        "gat" = c(3,1),
        "norm" = NULL,
        "std" = 4,
        "sstd" = 4,
        "skstd" = 4,
        "ged" = 2)

    # Conditional distribution
    cond.dist = match.arg(cond.dist)

    # Default Model (AR(1)):
    initialDelta = NULL
    if(!is.null(model$alpha) && is.null(model$delta)) # Garch model
    {
        if( any ( cond.dist == c("stableS0", "stableS1", "stableS2") ) )
            initialDelta = 1
        else 
            initialDelta = 2
    }
    	
    control = list(
        omega = 1,
        alpha = NULL,
        gamma = NULL,
        beta = NULL,
        mu = NULL,
        ar = NULL,
        ma = NULL,
        delta = initialDelta,
        skew = skew[[cond.dist]],
        shape = shape[[cond.dist]]
        )

    # Update Control:
    control[names(model)] <- model
    model <- control

    # Model Orders:
    order.ar = length(model$ar)
    order.ma = length(model$ma)
    order.alpha = length(model$alpha)
    order.beta = length(model$beta)
    order.max = max(order.ar, order.ma, order.alpha, order.beta)
    
    
    # Error treatment of input parameters:   
    if((order.alpha == 0 && order.beta != 0))
    	stop("In a Garch(p,q)/Aparch(p,q) model we must have p > 0")  
       
    if(!is.null(model$delta)){
    	if(length(model$delta) > 1)
       stop("The parameter 'delta' must be a single number")
    	 if(!(model$delta > 0))
      	stop("The parameter 'delta' must be > 0.")
    }
    
    if( (length(model$gamma) != 0) && (length(model$alpha) != 0)) # means aparch model
    {
      if(length(model$alpha) != length(model$gamma))
        stop("'alpha' and 'gamma' must have the same size for APARCH models")    
    }
   
    if(!is.null(model$gamma)){
    	if(sum(!(abs(model$gamma)<1)) > 0) # all gamma in (-1,1)
        stop("The parameter 'gamma' must be in the range -1 < gamma < 1")     	
    }
             
    if(sum(model$alpha < 0) > 0) # all alpha in [0,+infty)
        stop("The parameter 'alpha' must be greater than or equal 0")          
        
        
    if(sum(model$beta < 0) > 0) # all beta in [0,+infty)
        stop("The parameter 'beta' must be greater than or equal 0")
        
    if(!(model$omega > 0)) # omega > 0
        stop("The parameter 'omega' must be > 0")        
            
    
    if(is.null(model$ar) && is.null(model$ma) && is.null(model$alpha))
    	stop("The model parameters were not specified correctly.")
    	
    if(!is.null(model$delta) && is.null(model$alpha)) 	
       stop("The parameter delta should be only specified for a GARCH or APARCH model.")
       
        # if we have a presample check if it has the correct range.       
    if(is.matrix(presample) && !is.null(presample)){
    	if( dim(presample)[2] != 3 || dim(presample)[1] < order.max)
    	   stop(cat("The presample object should be a matrix with three columns formated as: \n [Innovations, Conditional Variance, Time Series] with dimensions \n l x 3, where l = max(m,n,p,q). "))
    	if( dim(presample)[1] != order.max )
    	{
          warning(cat("The number of columns of the Presample matrix is \n bigger than l = max(m,n,p,q). The simulated series \n will use only the first 'l' columns"))
          presample = as.matrix ( presample[1:order.max,] )
          if(order.max == 1)
              presample = t(presample)
    	}
    }
      	
    # Compose Mean Formula Object:
    formula.mean = ""
    if (order.ar == 0 && order.ma == 0) {
        formula.mean = ""
    }
    else {
        formula.mean = paste ("arma(", as.character(order.ar), ", ",
            as.character(order.ma), ")", sep = "")
    }
    
    # Compose Variance Formula Object:
    formula.var = ""
    if (order.alpha > 0) formula.var = "garch" 

    if(!is.null(model$alpha)){ # decide if we have aparch model
    	if (!is.null(model$gamma)){
    		if(sum(model$gamma == 0) != length(model$gamma))
    		   formula.var = "aparch"
    	} else {
    	   if (model$delta != 2 && ! any ( cond.dist == c("stableS0", "stableS1", "stableS2") )) 
            formula.var = "aparch" # gamma = 0 and delta != 0 we get powergarch model   	
    	   if (model$delta != 1 && any ( cond.dist == c("stableS0", "stableS1", "stableS2") )) 
    	      formula.var = "aparch" 
      }
    }
   
    if (order.alpha == 0 && order.beta == 0) {
        formula.var = formula.var
    }
    if (order.alpha > 0) {
        formula.var = paste(formula.var, "(", as.character(order.alpha),
                            ", ", as.character(order.beta), ")", sep = "")
    }

    # Compose Mean-Variance Formula Object:
    if (formula.mean == "") {
        formula = as.formula(paste("~", formula.var))
    } 
    if (formula.var == "") {
        formula = as.formula(paste("~", formula.mean))
    }
    if ((formula.mean != "") && (formula.var != "")){
        formula = as.formula(paste("~", formula.mean, "+", formula.var))
    }


    # Stop if the user specified a pure arma model
    if(formula.var == "")
        stop("Pure ARMA model not allowed")


    # Add NULL default entries:
    if (is.null(model$mu)) model$mu = 0
    if (is.null(model$ar)) model$ar = 0
    if (is.null(model$ma)) model$ma = 0
    if (is.null(model$gamma)) model$gamma = rep(0, times = order.alpha)

    # Seed:
    if (is.null(rseed)) {
      rseed = 0
    }
    else {
      set.seed(rseed)
    }

   # Define Missing Presample:
	 persistency = 1-sum(model$alpha)-sum(model$beta)
	 if(persistency*(1-persistency) < 0) # avoid to construct a presample with negative conditional variance.
	    persistency = 0.1 
	 if(!is.null(presample) && is.matrix(presample)){
	     z = presample[, 1]
	     h = presample[, 2]
	     y = presample[, 3]
	     if(sum(!(h > 0)))
	        stop("Conditional Variance column can have only strictly positive numbers")
	 }else{
	     z = rnorm(n = order.max)
	     h = rep(model$omega*(1 + persistency*(1-persistency)), times = order.max)
	     y = rep(model$mu, times = order.max)
	 }
	 presample = cbind(z, h, y)

    
    # Result: 
    new("GEVSTABLEGARCHSPEC",
        call = match.call(),
        formula = formula,
        model = list(omega = model$omega, alpha = model$alpha,
            gamma = model$gamma, beta = model$beta, mu = model$mu,
            ar = model$ar, ma = model$ma, delta = model$delta,
            skew = model$skew, shape = model$shape),
        presample = as.matrix(presample),
        distribution = as.character(cond.dist),
        rseed = as.numeric(rseed)
    )
}


################################################################################




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


# Copyrights (C) 2015, Thiago do Rego Sousa <thiagoestatistico@gmail.com>
# This is a modified version of the code contained inside file 
# garch-Spec.R from package fGarch, version 3010.82.

# Copyrights (C)
# for this R-port:
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file



################################################################################
# FUNCTION:               SIMULATION:
#  gsSim            Simulates a GARCH/APARCH process with GEV or stable
#						        conditional distribution
################################################################################


gsSim <-
    function(spec = gsSpec(), n = 100, n.start = 100)
{
    # A function originally implemented by Diethelm Wuertz and modified
    # to be used inside package GEVStableGarch. See the latest copyright notice. 
      
    # Description:
    #   Simulates a time series process from the GARCH family

    # Arguments:
    #   model - a specification object of class 'GEVSTABLEGARCHSPEC' as
    #     returned by the function \code{gsSpec}:
    #     ar - a vector of autoregressive coefficients of
    #       length m for the ARMA specification,
    #     ma - a vector of moving average coefficients of
    #       length n for the ARMA specification,
    #     omega - the variance value for GARCH/APARCH
    #       specification,
    #     alpha - a vector of autoregressive coefficients
    #       of length p for the GARCH/APARCH specification,
    #     gamma - a vector of leverage coefficients of
    #       length p for the APARCH specification,
    #     beta - a vector of moving average coefficients of
    #       length q for the GARCH/APARCH specification,
    #     mu - the intercept for ARMA specification (mean=mu/(1-sum(ar))),
    #     delta - the exponent value used in the variance
    #       equation.
    #     skew - a numeric value for the skew parameter.
    #     shape - a numeric value for the shape parameter.
    #   n - an integer, the length of the series
    #   n.start - the length of the warm-up sequence to reduce the
    #     effect of initial conditions.
     
    # Return:
    #   ans - An object returned by function timeDate with the simulated sample
    #     path of the specified ARMA-GARCH/APARCH model.
    
    # FUNCTION: 

	# Error treatment of input parameters
	if(n < 2)
	   stop("The parameter 'n' must be > 2")

    # Specification:
    stopifnot(class(spec) == "GEVSTABLEGARCHSPEC")
    model = spec@model

    # Random Seed:
    if (spec@rseed != 0) set.seed(spec@rseed)

    # Enlarge Series:
    n = n + n.start

    # Create Innovations:
	  if (spec@distribution == "stableS0")
	      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 0)
  
	  if (spec@distribution == "stableS1")
	      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 1)
  
	  if (spec@distribution == "stableS2")
	      z = stabledist::rstable(n = n, alpha = model$shape, beta = model$skew, pm = 2)
  
    if (spec@distribution == "gev")
        z = rgev(n, xi = model$shape)
  
	  if (spec@distribution == "gat") 
	      z = rgat(n, nu = model$shape[1], d = model$shape[2], xi = model$skew)
  
    if (spec@distribution == "norm")
        z = rnorm(n)
  
    if (spec@distribution == "std")
        z = rstd(n, nu = model$shape)
  
	  if (spec@distribution == "sstd") 
	      z = rsstd(n, nu = model$shape, xi = model$skew)
  
	  if (spec@distribution == "skstd") 
	      z = rskstd(n, nu = model$shape, xi = model$skew)
  
	  if (spec@distribution == "ged") 
	      z = rged(n, nu = model$shape)


    # Expand to whole Sample:    NAO ENTENDI PORQUE USAR A FUNCAO rev()???
    delta = model$delta
  
	  z = c(rev(spec@presample[, 1]), z)
	  h = c(rev(spec@presample[, 2]), rep(NA, times = n))
	  y = c(rev(spec@presample[, 3]), rep(NA, times = n))
	  m = length(spec@presample[, 1])
	  names(z) = names(h) = names(y) = NULL
    
    # Determine Coefficients:
    mu = model$mu
    ar = model$ar
    ma = model$ma
    omega = model$omega
    alpha = model$alpha
    gamma = model$gamma
    beta = model$beta
    deltainv = 1/delta

    # Determine Orders:
    order.ar = length(ar)
    order.ma = length(ma)
    order.alpha = length(alpha)
    order.beta = length(beta)

    # Iterate GARCH / APARCH Model and create Sample:
	  # print(c(omega,alpha,gamma,beta,delta))
  	eps = h^deltainv*z   # here the variable 'h' represents the process '(sigma_t)^delta'
  	for (i in (m+1):(n+m)) {
     	 	h[i] =  omega +
          	sum(alpha*(abs(eps[i-(1:order.alpha)]) -
              gamma*(eps[i-(1:order.alpha)]))^delta) +
          	sum(beta*h[i-(1:order.beta)])
        
          
          
      	eps[i] = h[i]^deltainv * z[i]
      	y[i] =
          	sum(ar*y[i-(1:order.ar)]) +
          	sum(ma*eps[i-(1:order.ma)]) + eps[i]
  	}
	y = y +  mu
  	# Sample:
  	data = cbind(
     	 	z = z[(m+1):(n+m)],
      	sigma = h[(m+1):(n+m)]^deltainv,
      	y = y[(m+1):(n+m)])    	
    
    rownames(data) = as.character(1:n)
    if(n.start > 0)
    	  data = data[-(1:n.start),]


    # Return Values:
    from <-
        timeDate(format(Sys.time(), format = "%Y-%m-%d")) - NROW(data)*24*3600
    charvec  <- timeSequence(from = from, length.out = NROW(data))
    ans <- timeSeries(data = data[, c(3,2,1)], charvec = charvec)
    colnames(ans) <- c("Series", "Volatility", "Innovations")    
    attr(ans, "control") <- list(gsSpec = spec)

    # Return Value:
    ans
}



################################################################################
