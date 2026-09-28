#' gsSim
#' 
#' Simulate ARMA-GARCH/APARCH model
#' 
#' details
#' 
#' @param model a list, specification of ARMA-GARCH/APARCH model
#' @param n integer, length of time series
#' @param burnin integer, burn-in period used to reduce dependence on initial values. Depends on how dependent the process is.  
#' @return a list, including tree time series: "Series", "Volatility", "Innovations"
#' @references 
#' 
#' Definition 3.1.2(The ARMA (p,q) Process). TSTM Brockwell pg 78
#' 
#' Eq. 7.2.1 Francq Zarqoian, Def. of ARMA-GARCH models. Simulate independently the process (ht,eta,zt) and Xt. 
#' 
#' Same Francq book, Theorem 7.4 (Consistency of the QMLE). 
#' 
#' @export 
gsSim <- function(model, n = 100, burnin = 1000, seed)
{

	  # error treatment of input parameters

    # set seed
    if(!missing(seed))
      set.seed(rseed)
  
    # alongate series

    # innovations
	  if (model$spec$cond_dist == "stable")
	      z = stabledist::rstable(n = n, alpha = model$params$dist_params$alpha, beta = model$params$dist_params$beta, pm = 0)
  
    if (model$spec$cond_dist == "gev")
        z = rgev(n, xi = model$params$dist_params$xi)
  
	  if (spec@distribution == "gat") 
	      z = rgat(n, nu = model$params$dist_params$nu, d = model$params$dist_params$d, xi = model$params$dist_params$xi)
  
    if (spec@distribution == "norm")
        z = rnorm(n)

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
    
  	# remove starting series
    rownames(data) = as.character(1:n)
    if(n.start > 0)
    	  data = data[-(1:n.start),]

    # Return Values
    data <- data[, c(3,2,1)]
    colnames(data) <- c("Series", "Volatility", "Innovations")    

    # Return Value
    data
}
