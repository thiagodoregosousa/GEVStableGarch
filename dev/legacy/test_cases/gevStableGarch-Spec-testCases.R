
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
# TEST CASES FOF FUNCTION:               SPECIFICATION:
#  gsSpec               Useful for finding erros. Please add the expected 
#							               output of the test case whenever possible. 
#							               This will help to trace erros during debugging.
################################################################################


# General Testing
spec <- gsSpec(model = list(), presample = a,cond.dist = ,rseed = NULL)
spec <- gsSpec(model = list(delta = 2), presample = a,cond.dist = ,rseed = NULL)
spec <- gsSpec(model = list(delta = 2), presample = NULL,cond.dist = ,rseed = NULL) 
spec <- gsSpec(model = list(ar = 1, delta = 2),cond.dist = ,rseed = NULL)  

# Testing the FORMULA (garch, arma, arch, etc)
spec <- gsSpec(model = list(ar = 1, alpha = 4, beta = 3, delta = 2), presample = NULL,cond.dist = ,rseed = NULL)   
spec <- gsSpec(model = list(ar = 1, alpha = -4, beta = 3, delta = 2), 
presample = NULL,cond.dist = ,rseed = NULL)    
spec <- gsSpec(model = list(ar = 1, alpha = 4, beta = 3, delta = 2), 
presample = NULL,cond.dist = NULL,rseed = NULL)      
spec <- gsSpec(model = list(mu = 3,ar = 1, alpha = 4, beta = 3, delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)         
# ARMA(1,1)-ARCH(1)
spec <- gsSpec(model = list(mu = 3,ar = 1, ma = 4,alpha = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)   
# AR(2)-ARCH(1)
spec <- gsSpec(model = list(mu = 3,ar = c(1,2),alpha = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)       
# MA(2)-ARCH(1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)     
# MA(2)-APARCH(1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = 3,delta = 1), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# MA(2)-APARCH(0,1): expect error because alpha was not specified.
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),beta = 3,delta = 1), 
presample = NULL,cond.dist = c("norm"),rseed = 3)    
# MA(2)-APARCH(1,1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = 3,beta = 3,delta = 1), 
presample = NULL,cond.dist = c("norm"),rseed = 3)  
# MA(2)-APARCH(2,1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = c(3,3),
gamma = c(-0.5,0),beta = 3,delta = 1), 
presample = NULL,cond.dist = c("norm"),rseed = 3)  
# MA(2)-APARCH(2,1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = c(3,3),
gamma = c(-0.5,0),beta = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# MA(2)-GARCH(2,1)
spec <- gsSpec(model = list(mu = 3,ma = c(1,2),alpha = c(3,3),
gamma = c(0,0),beta = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# AR(2)-GARCH(2,1)
spec <- gsSpec(model = list(mu = 3,ar = c(1,2),alpha = c(3,3),
gamma = c(0,0),beta = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# AR(2)-APARCH(2,1)
spec <- gsSpec(model = list(mu = 3,ar = c(1,2),alpha = c(3,3),
gamma = c(0,0.4),beta = 3,delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARMA(2,3)-APARCH(2,2)
spec <- gsSpec(model = list(ar = c(1,2),ma = c(3,3,3), alpha = c(3,3),
gamma = c(0,0.4),beta = c(3,3),delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# ARCH(1): expects to have an arch(1) model with delta = 2.
spec <- gsSpec(model = list(alpha = c(3), delta = 2), 
presample = NULL,cond.dist = c("norm"),rseed = 3)
# APARCH(1): expects an aparch(1) model with delta = 0.4.
spec <- gsSpec(model = list(alpha = c(3), delta = 0.4), 
presample = NULL,cond.dist = c("gev"),rseed = 3)    

# Testing the Presample matrix
spec <- gsSpec(model = list(alpha = c(3,3,3), delta = 0.4), 
presample = NULL,cond.dist = c("gev"),rseed = 3)   
spec <- gsSpec(model = list(alpha = c(0.99), delta = 0.4), 
presample = NULL,cond.dist = c("gev"),rseed = 4)   
# size of presample greater than needed. Expects error about dimension
pre <- matrix(4,nrow = 3,ncol = 3)
spec <- gsSpec(model = list(alpha = c(0.99), delta = 0.4), 
presample = pre, cond.dist = c("gev"),rseed = 4)
# size of smaller than needed. Expects error about dimension
pre <- matrix(4,nrow = 1,ncol = 3)
spec <- gsSpec(model = list(alpha = c(0.4,4), delta = 0.4), 
presample = pre, cond.dist = c("gev"),rseed = 4)
# exact presample size. OK
pre <- matrix(4,nrow = 1,ncol = 3)
spec <- gsSpec(model = list(alpha = c(0.4), delta = 0.4), 
presample = pre, cond.dist = c("gev"),rseed = 4)
# negattive conditional variance
pre <- matrix(4,nrow = 3,ncol = 3)
pre[1,2] = 0
spec <- gsSpec(model = list(alpha = c(0.4,3,4), delta = 0.4), 
presample = pre, cond.dist = c("gev"),rseed = 4)


################################################################################

