
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
#  .armaDist               Computes the ARMA log likelihood for
#                          several conditional distributions
################################################################################


.armaDist <- 
  function(z, sigma = 1, shape = 1.5, skew = 0, 
           cond.dist = c("stable", "gev", "gat", "norm"), 
           TOLG = 1e-8) 
  {
    # Description:
    #   Calculates the likelihood function for a vector of points (z)
    #   according to the specified distribution.
    #   REMARKS: We choose not to avoid calling this function with out of bound parameters
    #   because in these cases we simply return a huge likelihood to guide 
    #   the optimization algorithm.
    
    # Arguments:
    #   z - vector of points to calculate the llh.
    #   h - the sclale parameter
    #   cond.dist - name of the conditional distribution, one of
    #       "stable", "gev", "gat", "norm" (for testing only)
    #   TOLG - general tolerance for arma-garch parameters. 
    #   In the beggining it was set to 1e-5
    
    # FUNCTION:  
    
    # Error treatment of input parameters
    cond.dist = match.arg(cond.dist)
    if(sum(is.na(sigma)) > 0 || min(sigma) == 0)
    {
      warning("NA or zero element found in vector sigma")
      return(1e99)
    }
    
    # normal conditional distribution
    if(cond.dist == "norm")
      return(-sum(log(dnorm(x = z, sd = sigma))))
    
    # GAt distribution
    if(cond.dist == "gat")
    {
      if(!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0))
      {
        return(1e99)
      }
      
      nu = shape[1]
      d = shape[2]
      xi = skew
      
      return(-sum(log(dgat(x = z/sigma, nu = nu, d = d, xi = xi)/sigma)))        
    }
    
    # GEV conditional distribution
    if(cond.dist == "gev")
    {
      if(abs(shape) < 1e-6)
        stop("shape parameter from GEV is to small. abs(shape) < 1e-6")
      sig <- sigma
      xi <- shape
      arg <- z/sig
      y <- 1 + xi * arg
      if(sum(is.na(y)) || sum(is.nan(y)) ||
           sum(is.infinite(y))) 
        return(1e99)
      gev.cond.gD <- any(y < TOLG ) || abs(xi) < 1e-6
      if(gev.cond.gD)
      {
        return(1e99)
      }
      llh <- sum(log(sig)) + sum(y^(-1/xi)) + sum(log(y))*(1/xi + 1)
      return(llh)
    }
    
    # stable conditional distribution
    if(cond.dist == "stable")
    {
      if( !(shape > 1) || !(shape < 2) || !(abs(skew) < 1))
        stop("Invalid shape or skew in stable parameters. 1 < shape < 2 and |skew| < 1")
      sig <- sigma
      alpha.stable <- shape
      beta.stable <- skew
      arg <- z/sig
      y <- arg
      if(sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y)))
        return(1e99)
      dens.stable <- stable::dstable.quick(y,alpha = alpha.stable, 
                                           beta = beta.stable, gamma = 1, delta = 0, param = 1)
      #dens.stable <- GSgarch.dstable(y,alpha = alpha.stable,
      #               beta = beta.stable, gamma = 1,delta = 0, param = 2)
      llh <- sum(log(sig[sig>0])) - 
        sum(log(dens.stable[dens.stable>0]))
      return(llh)
    }
  }



################################################################################

