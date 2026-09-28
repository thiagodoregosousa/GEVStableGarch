
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
# TEST CASES FOF FUNCTION: 
# 
#  gsFit                
################################################################################
library(fGarch)
library(GEVStableGarch)
library(Rsolnp)
library(stable)
library(stabledist)
library(skewt)
data(dem2gbp)
x = dem2gbp[,1]


# ------------------------------------------------------------------------------
# Testing different type of datasets as input
# ------------------------------------------------------------------------------



x <- c("asdf",rnorm(1000))
x <- c(rnorm(100),NA,rnorm(200))
x <- c(rnorm(100),-Inf)
x <- c(NULL)
x <- rnorm(100)
gsFit(data = x , formula = ~arma(1,1)+garch(1,1),
            cond.dist = "norm", include.mean = TRUE, 
            algorithm = "sqp")



# ------------------------------------------------------------------------------
# Comparison with package GEVStableGarch from CRAN
# ------------------------------------------------------------------------------



# with dem2gbp[, 1] dataset
# Instructions: Clear the Workspace, 
# Load the packages and run 'model1', 2 and so on.
# Then load functions from the new version of the package 
# and run 'fit1', 2 and so on. Finally, compare the estimated
# parameters. 
library(fGarch)
library(GEVStableGarch)
data(dem2gbp)
x = dem2gbp[,1]

# garch(1,1)-norm-intercept
model1 = gsFit(data = x, 0,0,1,1, cond.dist = "norm", 
                     intercept = TRUE, algorithm = "nlminb+nm")

fit1 <- gsFit(data = x , formula = ~garch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

# garch(1,1)-std-intercept
model2 = gsFit(data = x, 0,0,1,1, cond.dist = "std", 
                     intercept = TRUE, algorithm = "nlminb+nm")

fit2 <- gsFit(data = x , formula = ~garch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

# garch(1,1)-gev-intercept
model3 = gsFit(data = x, 0,0,1,1, cond.dist = "gev", 
                     intercept = TRUE, algorithm = "nlminb+nm")

fit3 <- gsFit(data = x , formula = ~garch(1,1),
                    cond.dist = "gev", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

# arma(1,1)-aparch(1,1)-std-intercept
model4 = gsFit(data = x, 1,1,1,1, cond.dist = "std", 
                     intercept = TRUE, algorithm = "nlminb+nm",APARCH = TRUE)

fit4 <- gsFit(data = x , formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1@fit$par
fit1@fit$par
model2@fit$par
fit2@fit$par
model3@fit$par
fit3@fit$par
model4@fit$par
fit4@fit$par



# ------------------------------------------------------------------------------
# Comparison between GEVStableGarch and fGarch fit functions
# ------------------------------------------------------------------------------



library(fGarch)
data(dem2gbp)
x = dem2gbp[, 1]
library(Rsolnp)
# garch(1,1)-norm-intercept
fgarch_garch11_norm <- garchFit(data = x, formula = ~garch(1,1),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb+nm")@fit$params$params
gevstablegarch_garch11_norm <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "sqp.restriction")@fit$par
fgarch_garch11_norm["mu"]
gevstablegarch_garch11_norm["mu"]





# garch(1,1)-std-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                      cond.dist = "std", include.mean = TRUE,
                      algorithm = "nlminb+nm")
model1 <- gsFit(data = x , formula = ~garch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par

# garch(2,2)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(2,2),
                 cond.dist = "std", include.mean = TRUE,
                 algorithm = "nlminb+nm")


model1 <- gsFit(data = x , formula = ~garch(2,2),
                      cond.dist = "std", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par

# garch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb+nm")


model1 <- gsFit(data = x , formula = ~garch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par


# garch(1,1)-sstd-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                 cond.dist = "sstd", include.mean = TRUE,
                 algorithm = "nlminb+nm")


model1 <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "sstd", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par

# garch(1,1)-skstd-intercept

model1 <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "skstd", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
model1@fit$par


# aparch(1,1)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,1),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb+nm")


model1 <- gsFit(data = x , formula = ~aparch(1,1),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par

# aparch(1,0)-norm-intercept
fit1 <- garchFit(data = x, formula = ~aparch(1,0),
                 cond.dist = "norm", include.mean = TRUE,
                 algorithm = "nlminb+nm")


model1 <- gsFit(data = x , formula = ~aparch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "nlminb+nm", DEBUG = FALSE)
fit1@fit$par-model1@fit$par


# garch(1,1)-ged-intercept
fit1 <- garchFit(data = x, formula = ~garch(1,1),
                 cond.dist = "ged", include.mean = TRUE,
                 algorithm = "nlminb+nm")
model1 <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "ged", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par
fit1@fit$llh
model1@fit$llh

# aparch(3,2)-ged-intercept
fit1 <- garchFit(data = x, formula = ~aparch(3,2),
                 cond.dist = "ged", include.mean = TRUE,
                 algorithm = "nlminb+nm")
model1 <- gsFit(data = x , formula = ~aparch(3,2),
                      cond.dist = "ged", include.mean = TRUE, 
                      algorithm = "nlminb+nm")
fit1@fit$par-model1@fit$par
fit1@fit$llh
model1@fit$llh


# garch(1,1)-GAt-intercept
model1 <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "gat", include.mean = TRUE,
                      algorithm = "nlminb+nm")
model1@fit$par
model1@fit$llh


# aparch(1,1)-GAt-intercept-nlminb+nm
model1 <- gsFit(data = x , formula = ~aparch(1,1),
                      cond.dist = "gat", include.mean = TRUE, DEBUG = TRUE,
                      algorithm = "nlminb+nm")

# garch(1,1)-GAt-intercept-sqp
model1 <- gsFit(data = x , formula = ~garch(1,1),
                      cond.dist = "gat", include.mean = TRUE, DEBUG = TRUE,
                      algorithm = "sqp")

model1@fit$par
model1@fit$llh


# arma(1,1)-garch(1,1)-std-intercept
fit1 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par

# arma(5,0)-garch(1,0)-std-intercept
fit1 <- gsFit(data = x, formula = ~arma(5,0)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(5,0)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par

# arma(0,5)-garch(1,0)-std-intercept
fit1 <- gsFit(data = x, formula = ~arma(0,5)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(0,5)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par
fit1@fit$par
model1@fit$par-fit1@fit$par


# arma(1,1)-garch(1,0)-norm-intercept
fit1 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(1,1)+garch(1,0),
                   cond.dist = "norm", include.mean = TRUE)
model1@fit$par-fit1@fit$par


# arma(1,1)-aparch(1,1)-std-intercept
fit1 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "std", include.mean = TRUE, 
                    algorithm = "sqp")

model1 <- garchFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                   cond.dist = "std", include.mean = TRUE)
(model1@fit$par-fit1@fit$par)/fit1@fit$par


# arma(1,1)-aparch(1,1)-sstd-intercept
fit1 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                    cond.dist = "sstd", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                   cond.dist = "sstd", include.mean = TRUE)
(model1@fit$par-fit1@fit$par)/fit1@fit$par
cbind(model1@fit$par,fit1@fit$par)


# arma(2,2)-aparch(2,2)-sstd-intercept
fit1 <- gsFit(data = x, formula = ~arma(2,2)+aparch(2,2),
              cond.dist = "sstd", include.mean = TRUE, 
              algorithm = "nlminb+nm")

model1 <- garchFit(data = x, formula = ~arma(2,2)+aparch(2,2),
                   cond.dist = "sstd", include.mean = TRUE)
(model1@fit$par-fit1@fit$par)/fit1@fit$par
cbind(model1@fit$par,fit1@fit$par)

# arma(2,2)-aparch(2,2)-skstd-intercept
fit1 <- gsFit(data = x, formula = ~arma(2,2)+aparch(2,2),
              cond.dist = "skstd", include.mean = TRUE, 
              algorithm = "nlminb+nm")
model1 <- garchFit(data = x, formula = ~arma(2,2)+aparch(2,2),
                   cond.dist = "sstd", include.mean = TRUE)

comparison = cbind(fit1@fit$par,model1@fit$par); colnames(comparison) = c("gsFit","garchFit")
comparison
fit1@fit$llh
model1@fit$llh


# arma(1,1)-garch(1,1)-GAt-intercept-nlminb+nm
model1 <- gsFit(data = x , formula = ~arma(1,1)+garch(1,1),
                      cond.dist = "gat", include.mean = TRUE, DEBUG = TRUE,
                      algorithm = "nlminb+nm")


# arma(0,1)-garch(1,1)-norm-intercept
data(sp500dge)
x = 100*sp500dge[, 1]
fit1 <- gsFit(data = x, formula = ~arma(0,1)+aparch(1,1),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "sqp")
model1 <- garchFit(data = x, formula = ~arma(0,1)+aparch(1,1),
                   cond.dist = "norm", include.mean = TRUE, 
                   algorithm = "nlminb+nm")
model1@fit$par-fit1@fit$par



# ------------------------------------------------------------------------------
# Fitting models using different Algorithms ('nlminb+nm' and 'sqp')
# ------------------------------------------------------------------------------



# arma(1,1)-aparch(1,1)-std
model1 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "std", include.mean = TRUE, 
                algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "std", include.mean = TRUE, 
                algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "std", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)


# arma(1,0)-aparch(1,0)-norm
model1 <- gsFit(data = x, formula = ~arma(0,1)+aparch(1,0),
                    cond.dist = "norm", include.mean = TRUE, 
                    algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(0,1)+aparch(1,0),
                      cond.dist = "norm", include.mean = TRUE, 
                      algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(0,1)+aparch(1,0),
                cond.dist = "norm", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)


# arma(1,1)-aparch(1,1)-sstd
model1 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(1,1)+aparch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)


# arma(1,1)-garch(1,1)-GAt
model1 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "gat", include.mean = TRUE, 
                algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "gat", include.mean = TRUE, 
                algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "gat", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)


# arma(1,1)-garch(1,1)-ged
model1 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "ged", include.mean = TRUE, 
                algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "ged", include.mean = TRUE, 
                algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "ged", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)


# arma(1,1)-garch(1,1)-sstd
model1 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "nlminb+nm")

model2 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "sqp")

model3 <- gsFit(data = x, formula = ~arma(1,1)+garch(1,1),
                cond.dist = "sstd", include.mean = TRUE, 
                algorithm = "sqp.restriction")
cbind(model1@fit$par,model2@fit$par,model3@fit$par)



# ------------------------------------------------------------------------------
# Testing the sqp.restriction algorithm
# ------------------------------------------------------------------------------


library(stable)
library(fGarch)
data(dem2gbp)
library(skewt)
x = dem2gbp[, 1]
# c("stableS1", "gev", "gat", "norm", "std", "sstd", "skstd", "ged")


# garch(1,1)-stableS1-intercept-NAO FIZ AINDA
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "stableS2", include.mean = TRUE, 
              algorithm = "sqp", control = list( trace = 3, tol = 1e-5))


# garch(1,1)-gev-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "gev", include.mean = TRUE, 
              algorithm = "sqp.restriction",
              tolerance = list( TOLG = 1e-7, TOLSTABLE = 1e-2, TOLSTATIONARITY = 1e-3))

fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "gev", include.mean = TRUE, 
              algorithm = "sqp")


# garch(1,1)-GAt-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "gat", include.mean = TRUE, 
              algorithm = "sqp.restriction")

fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "gat", include.mean = TRUE, 
              algorithm = "sqp", control = list( trace = 3, tol = 1e-5))
gsMomentAparch(cond.dist = "gat", shape = fit1@fit$par[6:7], 
               skew = fit1@fit$par[5], gm = 0, delta = 2)*fit1@fit$par[3] + fit1@fit$par[4]


# garch(1,1)-norm-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "norm", include.mean = TRUE, 
              algorithm = "sqp.restriction", control = list( trace = 3, tol = 1e-5))
fit1@fit$par[3] + fit1@fit$par[4]


# garch(1,1)-skstd-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "skstd", include.mean = TRUE, 
              algorithm = "sqp.restriction", control = list( trace = 3, tol = 1e-5))

fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "skstd", include.mean = TRUE, 
              algorithm = "sqp", control = list( trace = 3, tol = 1e-5))

.stationarityAparch(model = list(alpha = fit1@fit$par[3], beta = fit1@fit$par[4],
                                 skew = fit1@fit$par[5], shape = fit1@fit$par[6], gm = 0,
                                 delta = 2), formula = .getFormula(~garch(1,1)),
                    cond.dist = "skstd")

# garch(1,1)-aparch-intercept
fit1 <- gsFit(data = x, formula = ~aparch(1,1),
              cond.dist = "skstd", include.mean = TRUE, 
              algorithm = "sqp.restriction")

fit1 <- gsFit(data = x, formula = ~aparch(1,1),
              cond.dist = "skstd", include.mean = TRUE, 
              algorithm = "sqp")

.stationarityAparch(model = list(alpha = fit1@fit$par[3], beta = fit1@fit$par[5],
                                 skew = fit1@fit$par[7], shape = fit1@fit$par[8], 
                                 gm = fit1@fit$par[4],
                                 delta = fit1@fit$par[6]), formula = .getFormula(~aparch(1,1)),
                    cond.dist = "skstd")


# garch(1,1)-ged-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "ged", include.mean = TRUE, 
              algorithm = "sqp.restriction")

fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "ged", include.mean = TRUE, 
              algorithm = "sqp")

fit1@fit$par[3] + fit1@fit$par[4]

# garch(1,1)-sstd-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "sstd", include.mean = TRUE, 
              algorithm = "sqp.restriction")

model1 <- garchFit(data = x, formula = ~garch(1,1),
         cond.dist = "sstd", include.mean = TRUE)

fit1@fit$par[3]+fit1@fit$par[4]
model1@fit$par[3]+model1@fit$par[4]



# ------------------------------------------------------------------------------
# Testing the tolerance parameter
# ------------------------------------------------------------------------------



library(fGarch)
data(dem2gbp)
x = dem2gbp[,1]
x = rnorm(1000)
# garch(1,1)-norm-intercept
fit1 <- gsFit(data = x, formula = ~garch(1,1),
              cond.dist = "sstd", include.mean = TRUE, 
              algorithm = "sqp.restriction", 
              tolerance = list (TOLG = 1e-8, TOLSTABLE = 1e-2, TOLSTATIONARITY = 0.035))
0.124833 + 0.883072
0.117935 + 0.881065
0.143533 + 0.816467
0.138609 + 0.826391






# ------------------------------------------------------------------------------
# Fitting models with GEV conditional distribution and several datasets
# ------------------------------------------------------------------------------
# SET CURRENT DATASET FOLDER
mac.folder = "/Users/thiago/Desktop/dataset/"
windows.folder = "c:/dataset/"
data.set.folder = mac.folder

# JOSE DATASETS
returnCACjose = read.csv(paste(data.set.folder,"jose/DAX.csv", sep = ""), sep = ";")[,2]
returnDJIAjose = read.csv(paste(data.set.folder,"jose/DJIA.csv", sep = ""), sep = ";")[,2]
returnPSI20jose = read.csv(paste(data.set.folder,"jose/PSI20.csv", sep = ""), sep = ";")[,2]

# CFRAIN DATASETS
returnDAXcfrain = read.table(file=paste(data.set.folder,"cfrain/DAX.dat", sep = ""),skip = 1)[1:11502,4]
returnFTSEcfrain = read.table(file=paste(data.set.folder,"cfrain/FTSE.dat", sep = ""),skip = 1)[1:4982,4]



data(dem2gbp)
x = dem2gbp[, 1]
x = returnCACjose
x = returnDJIAjose
x = returnPSI20jose
x = returnDAXcfrain
x = returnFTSEcfrain
x = (x^2*sign(x))
x = atan(x)
x = cos(x)
x = exp(x)
# garch(1,1)-gev-intercept

fit1 = gsFit(data = sp500[,1]+0.2, formula = ~garch(1,1),
              cond.dist = "gev", include.mean = TRUE, 
              algorithm = "sqp", DEBUG = TRUE)

fit2 = garchFit(data = x, formula = ~garch(1,1),
         cond.dist = "ged", include.mean = TRUE, algorithm = "nlminb")








################################################################################












