
###############

p_default_saltos <- read_excel("Documents/projetos-divcri/bgev paper/data_analysis/p_default_saltos.xlsx")
ret = p_default_saltos$p_naive_saltos_diff

model_fit = gsFit(data = ret , formula = ~arma(1,1)+garch(1,1), cond.dist = "stableS0")
residuals = model_fit@residuals/model_fit@sigma.t #model_fit@residuals

x = seq(min(residuals),max(residuals), 0.01)
y = stable_pdf(x, pars = c(model_fit@fit$par['shape1'], model_fit@fit$par['skew']))

hist(residuals, freq = FALSE, ylim = c(0,2*max(y)), breaks = 20, main = c(round(model_fit@fit$par['shape1'],3), round(model_fit@fit$par['skew'],3)))
lines(x,y)






###############
data(dem2gbp)
x = dem2gbp[, 1]
model_fit = gsFit(data = x , formula = ~garch(1,1), cond.dist = "stableS0")
residuals = model_fit@residuals/model_fit@sigma.t #model_fit@residuals

x = seq(-3,3,0.01)
y = stable_pdf(x, pars = c(model_fit@fit$par['shape1'], model_fit@fit$par['skew']))

hist(residuals, freq = FALSE, ylim = c(0,2*max(y)), breaks = 20, main = c(round(model_fit@fit$par['shape1'],3), round(model_fit@fit$par['skew'],3)))
lines(x,y)





