#install.packages("pacman")
library(pacman)
pacman::p_load(tidyverse, hms, scales, ggplot2, lubridate, dplyr, readr, gridExtra,
               copula, fBasics, StableEstim, stabledist, PerformanceAnalytics,
               extRemes, ismev, evmix, extremis, DT, gdtools, kableExtra, VGAM, evd, 
               fExtremes, graphics, ggExtra, patchwork, quantmod, xtable, VineCopula, skewt,
               copula, rugarch, fGarch, GEVStableGarch, moments, writexl, rmgarch, tictoc, parallel)

source("gevstable_garch_fast_stable.R")

# Importa série
series <- new.env()
getSymbols(Symbols = c("PETR4.SA", "VALE3.SA", "ABEV3.SA"), 
           env = series, src = "yahoo", 
           from =as.Date("2015-09-01"), to = as.Date("2025-08-29"))



# Ordem do modelo
arma_p = 1; arma_q = 0; 
aparch_p = 1; aparch_q = 1

formula_modelo <- as.formula(paste0("~ arma(", arma_p, ",", arma_q, ") + aparch(", aparch_p, ",", aparch_q, ")"))

#### Exemplo com pvalores consistentes ####
# Série completa de retornos
abev_ret_hi <- diff(log(Hi(series$ABEV3.SA))) %>% na.omit()
abev_ret_hi_num <- as.numeric(abev_ret_hi)

plot(abev_ret_hi, type = "l")

fit_stableS0_abev <- gsFit(data = abev_ret_hi_num, 
                      formula = formula_modelo, 
                      cond.dist = "stableS0")

fit_stableS0_abev@fit[["matcoef"]] #pvalores consistentes


#### Exemplo com pvalores inconsistentes ####
preco_window = Hi(series[["ABEV3.SA"]])[177:2166]
retorno_window <- diff(log(preco_window))
retorno_window_num <- as.numeric(na.omit(retorno_window))

plot(preco_window, type = "l")
plot(retorno_window, type = "l")


fit_stableS0_abev2 <- gsFit(data = retorno_window_num, 
                           formula = formula_modelo, 
                           cond.dist = "stableS0")

fit_stableS0_abev2@fit[["matcoef"]] #pvalores inconsistentes


# vamos estimar um modelo simplificado ar(1)-garch(1,1). Aqui também customizei 
# o solver utilizado com outros parâmetros de tolerância. Muuito importante aqui notar 
# estamos diante de problemas com existências de vários mínimos locais e ao estimar uma vez
# garantimos uma 'boa solução', mas não uma global. A qualidade dela avaliamos nos diagnósticos. 

fit_stableS0_abev3 <- gsFit(data = retorno_window_num, 
                            formula = ar(1)~garch(1,1), 
                            cond.dist = "stableS0", include.mean = FALSE, 
                            control = list(delta = 1e-10, 
                                          tol = 1e-10))

fit_stableS0_abev3@fit[["matcoef"]] #pvalores inconsistentes








# comparando com a normal
garchspec.norm <- ugarchspec(mean.model = list(armaOrder = c(arma_p,arma_q), 
                                               include.mean = T),
                             variance.model = list(model = "apARCH", 
                                                   garchOrder = c(aparch_p, aparch_q)),
                             distribution.model = "norm",
                             fixed.pars = list(delta = 2))

fit_norm_abev2 <- ugarchfit(data = retorno_window_num, 
                            spec = garchspec.norm, 
                            solver = "hybrid")

fit_norm_abev2@fit$matcoef

