library(stabledist)
library(libstable4u)
alpha = runif(1,0.1,1.99)
beta = runif(1,-0.99,0.99)
gamma = runif(1,0.01,10)
delta =runif(1,-10,10)
x = runif(1,-10,10)



pars = c(alpha, beta, gamma, delta)
(stable_pdf(x, pars) - stabledist::dstable(x, alpha, beta, gamma, delta, pm = 0)) < 0.00001
