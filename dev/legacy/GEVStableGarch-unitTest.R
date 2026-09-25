
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
# TEST THE MAIN FUNCTIONS OF PACKAGE GEVStableGarch: 
# 
# Note: There is also the possibility to enforce the tests when 
# checking the package with R CMD using package testthat.
################################################################################

library(Rsolnp)
library(stable)
library(stabledist)
library(skewt)


# Required packages
library(testthat)



# Create file to save the result of the test
log.file.directory = '/Users/thiago/downloads'
log.file.name = format(Sys.time(), "DATE-%m-%d-%y---TIME-%Hh%Mm%Ssec.md")
log.file.adress = paste(log.file.directory, "/", log.file.name, sep = '')

write.log.function = function(x) {
  write(x, file = log.file.adress,
        ncolumns = if(is.character(x)) 1 else 5,
        append = TRUE, sep = "")
}
write.log.function("===================================================")
write.log.function("Testing functions from package GEVStableGarch.")
write.log.function(paste("\nDate and time: ",log.file.name))
write.log.function("===================================================")



# ------------------------------------------------------------------------------



# GARCH(1,1)-NORMAL
# Comparing the estimated parameters of a 
# Garch(1,1) with normal distribution and the fGarch package
# Compare the results obtained from 3 different algorithms 
# againts the reported results from garchFit, G@RCH and Splus FinMetrics softwares.

data(dem2gbp)
x = dem2gbp[,1]

garch11.garchFit = c(-0.006190314, 0.010761384, 0.153134059, 0.805973748)
garch11.garch.software = c(-0.006184, 0.010760, 0.153407, 0.805874)
garch11.finmetrics = c(-0.006053, 0.010896, 0.15421, 0.80445)

garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~garch(1,1),
                          cond.dist = "norm", algorithm = "sqp")@fit$par

garch11.gevstablegarch.nlminb = gsFit(data = x, formula = ~garch(1,1),
                                         cond.dist = "norm", algorithm = "nlminb")@fit$par

garch11.gevstablegarch.nlminb.nm = gsFit(data = x, formula = ~garch(1,1),
                                   cond.dist = "norm", algorithm = "nlminb+nm")@fit$par

garch11.gevstablegarch.sqp.restriction = gsFit(data = x, formula = ~garch(1,1),
                                            cond.dist = "norm", algorithm = "sqp.restriction")@fit$par

write.log.function("\n\n-----\nDescription: GARCH(1,1)-NORMAL") 
write.log.function("Details: Testing the estimated parameters from the 
garch(1,1)-norm against garchFit, G@RCH and SPLUS FinMetrics...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garchFit)/garch11.gevstablegarch.sqp), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garch.software)/garch11.gevstablegarch.sqp), 1e-2)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.finmetrics)/garch11.gevstablegarch.sqp), 0.03)
write.log.function("Result: OK")

write.log.function("\n\n-----\nDescription: GARCH(1,1)-NORMAL-ALGORITHMS") 
write.log.function("Details: Testing the estimated parameters from the 
garch(1,1)-norm with several algorithms algorithms: sqp, nlminb, nlminb+nm and sqp.restriction...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb.nm), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.nlminb - garch11.gevstablegarch.nlminb)/garch11.gevstablegarch.nlminb.nm), 1e-4)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.sqp.restriction)/garch11.gevstablegarch.sqp), 1e-4)
write.log.function("Result: OK")



# ------------------------------------------------------------------------------



# GARCH(1,1)-STUDENT'S T
# Comparing the estimated parameters of a 
# Garch(1,1) with student's t conditional distribution
# againts the reported results from garchFit, G@RCH and Splus FinMetrics softwares.


data(dem2gbp)
x = dem2gbp[,1]

garch11.garchFit = c(0.002248696, 0.002319032, 0.124438126, 0.884653077, 4.118425040)
garch11.garch.software = c(0.002251,0.002322,0.124874,0.884481,4.112100)
garch11.finmetrics = c(0.001175, 0.001981, 0.120837, 0.887464, 4.163737)

garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~garch(1,1),
                                   cond.dist = "std", algorithm = "sqp")@fit$par

write.log.function("\n\n-----\nDescription: GARCH(1,1)-STUDENTS'T") 
write.log.function("Details: Testing the estimated parameters from the 
garch(1,1)-std against garchFit, G@RCH and SPLUS FinMetrics...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garchFit)/garch11.gevstablegarch.sqp), 1e-3)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garch.software)/garch11.gevstablegarch.sqp), 1e-2)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.finmetrics)/garch11.gevstablegarch.sqp), 0.5)
write.log.function("Result: OK")



# ------------------------------------------------------------------------------



# GARCH(1,1)-GED
# Comparing the estimated parameters of a 
# Garch(1,1) with GED conditional distribution
# againts the reported results from garchFit, G@RCH and Splus FinMetrics softwares.


data(dem2gbp)
x = dem2gbp[,1]

garch11.garchFit = c(0.001692458, 0.004478992, 0.130834320, 0.859286305, 1.149397930)
garch11.garch.software = c(0.001699, 0.004479, 0.131135, 0.859152, 1.149174)
garch11.finmetrics = c(0.000411, 0.004144, 0.131326, 0.860243, 1.147371)

garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~garch(1,1),
                                   cond.dist = "ged", algorithm = "sqp")@fit$par

write.log.function("\n\n-----\nDescription: GARCH(1,1)-GED") 
write.log.function("Details: Testing the estimated parameters from the 
garch(1,1)-ged against garchFit, G@RCH and SPLUS FinMetrics...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garchFit)/garch11.gevstablegarch.sqp), 1e-3)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.garch.software)/garch11.gevstablegarch.sqp), 1e-2)
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.finmetrics)/garch11.gevstablegarch.sqp), 0.8)
write.log.function("Result: OK")



# ------------------------------------------------------------------------------



# AR(1)-GARCH(1,1)-NORMAL-DING-GRANGER-ENGLE-DATASET
# Comparing the estimated parameters of a 
# Ar(1)-Garch(1,1) with normal conditional distribution for 
# the sp500 dataset from Ding, Granger and Engle (1983)
# againts the reported results from garchFit, G@RCH and Splus FinMetrics softwares.

x = 100*sp500dge[, 1]
ar1garch11.dge = c(0.00021, 0.145, 0.000014, 0.083, 0.373, 0.920, 1.43)

ar1garch11.garchFit = c(2.045478e-04, 1.446432e-01, 1.512069e-05, 8.362974e-02, 3.767485e-01, 9.199221e-01, 1.418821e+00  ) 
ar1garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~arma(0,1)+aparch(1,1),
                       cond.dist = "norm", algorithm = "sqp")@fit$par
ar1garch11.gevstablegarch.sqp[1] = 0.01*ar1garch11.gevstablegarch.sqp[1]
ar1garch11.gevstablegarch.sqp[3] = 0.01^(2/ar1garch11.gevstablegarch.sqp[7])*ar1garch11.gevstablegarch.sqp[3]
ar1garch11.finmetrics = c(0.000208, 0.1447, 0.0000159, 0.08375, 0.3710, 0.9195, 1.429)
ar1garch11.garch.software = c(0.000204, 0.1446, 0.0000150, 0.08377, 0.3765, 0.9199, 1.416)

write.log.function("\n\n-----\nDescription: AR(1)-GARCH(1,1)-NORMAL-DING-GRANGER-ENGLE-DATASET") 
write.log.function("Details: Testing the estimated parameters from the 
Ar(1)-garch(1,1)-normal against garchFit, G@RCH and SPLUS FinMetrics and the Ding et al. paper...")
expect_less_than( max(abs(ar1garch11.gevstablegarch.sqp - ar1garch11.dge)/ar1garch11.gevstablegarch.sqp), 0.1)
expect_less_than( max(abs(ar1garch11.gevstablegarch.sqp - ar1garch11.garchFit)/ar1garch11.gevstablegarch.sqp), 1e-5)
expect_less_than( max(abs(ar1garch11.gevstablegarch.sqp - ar1garch11.finmetrics)/ar1garch11.gevstablegarch.sqp), 0.5)
expect_less_than( max(abs(ar1garch11.gevstablegarch.sqp - ar1garch11.garch.software)/ar1garch11.gevstablegarch.sqp), 1e-2)
write.log.function("Result: OK")



# ------------------------------------------------------------------------------



# GARCH(1,1)-GEV
# Garch(1,1) compare the estimation of the GEV model for the dem2gbp
# with the estimated values from the previous version of package GEVStableGarch

data(dem2gbp)
x = dem2gbp[,1]

garch11.gevstablegarch.sqp.previous = c(-0.14581498, 0.05386922, 0.36623600, 0.48890513, -0.22038572)
garch11.gevstablegarch.sqp = gsFit(data = x, formula = ~garch(1,1),
                                   cond.dist = "gev", algorithm = "sqp")@fit$par

write.log.function("\n\n-----\nDescription: GARCH(1,1)-GEV") 
write.log.function("Details: Compare with the estimated parameters from the previous version of
                   package GEVStableGarch...")
expect_less_than( max(abs(garch11.gevstablegarch.sqp - garch11.gevstablegarch.sqp.previous)/garch11.gevstablegarch.sqp), 1e-5)
write.log.function("Result: OK")




