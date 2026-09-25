library(GEVStableGarch) # install manually skewt and Nlpsolv then download version 1.1. from here 
#  https://cran.r-project.org/src/contrib/Archive/GEVStableGarch/ and install manually
library(readxl)

# This examples uses the dem2gbp dataset to estimate
# an ARMA(1,1)-GARCH(1,1) with GEV conditional distribution.
data(dem2gbp)
x = dem2gbp[, 1]
gev.model = gsFit(data = x , formula = ~garch(1,1), cond.dist = "gev")

retornos_10anos <- read_excel("Documents/projetos-divcri/bgev paper/data_analysis/retornos_10anos.xlsx", 
                   col_types = c("date", "numeric", "numeric"))

petro = retornos_10anos$PETR4.SA.Adjusted
vale = retornos_10anos$VALE3.SA.Adjusted

gsFit(data = vale , formula = ~arma(1,1)+garch(1,1), cond.dist = "stableS0")




petro_best = gsSelect(data = petro, order.max = c(2,2,2,2), selection.criteria = "AIC", cond.dist = "stableS0", include.mean = TRUE)
                                                                                                                                                                                         
vale_best = gsSelect(data = vale, order.max = c(2,2,2,2), selection.criteria = "AIC", cond.dist = "stableS0", include.mean = TRUE)

