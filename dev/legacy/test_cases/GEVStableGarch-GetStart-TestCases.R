
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
#  .getStart           Test on input parameters for different models 
#    					                including pure ARMA models							               
################################################################################


#######
# Test the output for arma(1,1)-aparch(1,1) model
library(fExtremes); library(skewt)
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
mylist = c("stableS0", "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
    start = .getStart(data = x,m = 1,n = 1,p = 1,q = 1, 
                 AR = FALSE, MA = FALSE,
                 cond.dist = mylist[i]) 
    print(mylist[i])
    print(start)
    cat("\n\n")
} 

# Test the output for the arma(3,3)-aparch(3,3) model
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
mylist = c("stable", "gev", "gat", "norm", "std", "sstd", "skstd", "ged")
for( i in 1:length(mylist) )
{
  start = .getStart(data = x,m = 3,n = 3,p = 3,q = 2, 
            AR = FALSE, MA = FALSE,
            cond.dist = mylist[i]) 
  print(mylist[i])
  print(start)
} 
 

# Test start parameters for gev distribution
library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
.getStart(data = x,m = 1,n = 1,p = 1,q = 1, 
          AR = FALSE, MA = FALSE,
          cond.dist = "gev")  

