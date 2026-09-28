#' Generalized Asymmetric t (GAT) Distribution
#'
#' Functions to compute the density, distribution function, quantile function,
#' and to generate random variates for the Generalized Asymmetric t (GAT)
#' distribution.
#'
#' This distribution corresponds to the \code{t3} distribution described in
#' Paolella (1997) and Mittnik and Paolella (2000). The GAT family includes
#' the Student's t, Laplace, Cauchy, and normal distributions as special cases.
#' In particular, as \eqn{\nu \to \infty} the distribution approaches normality.
#'
#' @name GAT
#' @rdname GAT
#' @aliases gat dgat pgat qgat rgat
#'
#' @param x Numeric vector of values for calculating density. 
#' @param q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations for random generation.
#' @param mean Location parameter.
#' @param sd Scale parameter (must be > 0).
#' @param nu Tail (shape) parameter (must be > 0).
#' @param d Second shape parameter (must be > 0).
#' @param xi Asymmetry parameter (must be > 0; \eqn{\xi = 1} gives symmetry).
#' @param log Logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return
#' \item{dgat}{density values}
#' \item{pgat}{distribution function values}
#' \item{qgat}{quantile function values}
#' \item{rgat}{random variates}
#'
#' @references
#' Mittnik, S., Paolella, M. S. (2000).
#' Prediction of Financial Downside-Risk with Heavy-Tailed Conditional
#' Distributions.
#'
#' Paolella, M. (1997).
#' Tail Estimation and Conditional Modeling of Heteroskedastic Time-Series.
#' PhD Thesis, Institute of Statistics and Econometrics,
#' Christian Albrechts University of Kiel.
#' 
#' Bertocchi, M., Giacometti, R., Ortobelli, S., & Rachev, S. T. (2005).
#' The impact of different distributional hypothesis on returns in asset allocation.
#' \emph{Finance Letters}, 3(1), 17-27.
#'
#' @author Thiago do Rego Sousa
#'
#' @examples
#' par(mfrow = c(2, 2))
#' set.seed(1000)
#' r <- rgat(n = 1000)
#' plot(r, type = "l", main = "GAt Random Values")
#'
#' hist(r, probability = TRUE, border = "white")
#' x <- seq(min(r), max(r), length = 201)
#' lines(x, dgat(x), lwd = 2)
#'
#' plot(sort(r), (1:1000)/1000, main = "Probability", ylab = "Probability")
#' lines(x, pgat(x), lwd = 2)
#'
#' round(qgat(pgat(q = seq(-10, 10, by = 0.5))), 6)
#'
#' @rdname GAT
#' @export
dgat <- 
  function(x, mean = 0, sd = 1, nu = 2, d = 3, xi = 1, log = FALSE)
  {   
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Scaled distance from the mode: negative side multiplied by xi, positive divided by xi
    z = (x - mean ) / sd
    arg = ifelse(z < 0, -z * xi, z / xi)
    
    # Work on the log scale so far tails do not underflow
    log_k = -log( ( xi + 1/xi ) * 1/d * nu^(1/d) * beta(1/d, nu) )
    result = log_k - (nu + 1/d) * log1p(arg^d / nu) - log(sd)
    if(!log) result = exp(result)
    
    # Return Value
    result
  }


#' @rdname GAT
#' @export
pgat <- 
  function(q, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)
  {   
    
    # Params:
    if (length(mean) == 5) {
      xi = mean[5]
      d  = mean[4]
      nu = mean[3]
      sd = mean[2]
      mean = mean[1]
    }    
    
    # Error treatment of input parameters
    if(sd <= 0  || nu <= 0 || xi <= 0 || d <= 0)
      stop("Failed to verify condition:
           sd <= 0 || nu <= 0 || xi <= 0 || d <= 0")
    
    # Both tails use w = pw / (nu + pw) with the scaled distance pw. On the negative side
    # the upper beta tail replaces pbeta(1 - w, nu, 1/d), which loses all accuracy near 0
    z = (q - mean ) / sd
    neg = z <= 0
    pw = ifelse(neg, (-z * xi)^d, (z / xi)^d)
    w = ifelse(is.infinite(pw), 1, pw / (nu + pw))
    ifelse(neg,
           1/(1 + xi^2) * pbeta(w, 1/d, nu, lower.tail = FALSE),
           1/(1 + xi^2) + 1/(1 + xi^(-2)) * pbeta(w, 1/d, nu))
  }


#' @rdname GAT
#' @export
qgat <- 
  function(p, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    # Define auxiliary functions
    Lp <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( 1 + xi^2 ) * p, nu, 1/d)
    }
    Up <- function (p = p, nu = nu, d = d, xi = xi)
    {
        qbeta( ( p - 1 / ( 1 + xi^2 ) ) * ( 1 + xi^(-2) ), 1/d, nu)
    }
    
    # Compute quantiles located at (-Inf,0] and at (0,+Inf)
    F0 = pgat(0, mean = 0, sd = 1, nu = nu, d = d, xi = xi)
    n = length(p)
    result = rep(NA,n)
    indexLessThanF0 = which (p <= F0, arr.ind = TRUE)
    sizeIndex = length(indexLessThanF0)
    if(sizeIndex == 0) {
      
        U = Up (p = p, nu = nu, d = d, xi = xi)
        result = ( U * nu / (1 - U) )^( 1/d ) * xi
        
    } else if (sizeIndex == n) {
      
         L = Lp (p = p, nu = nu, d = d, xi = xi)
        result = - ( nu/L - nu )^( 1/d ) * 1/xi
        
    } else if (TRUE) {
      
        L = Lp (p = p[indexLessThanF0], nu = nu, d = d, xi = xi)
        U = Up (p = p[-indexLessThanF0], nu = nu, d = d, xi = xi)
        result[indexLessThanF0] = - ( nu/L - nu )^( 1/d ) * 1/xi
        result[-indexLessThanF0] = ( U * nu / (1 - U) )^( 1/d ) * xi 
    }
    
    # Return Value:
    result * sd + mean
  }


#' @rdname GAT
#' @export
rgat <-  
  function(n, mean = 0, sd = 1, nu = 2, d = 3, xi = 1)  
  {   
    
    randomUnif = runif(n = n, min = 0, max = 1)
    result = qgat(p = randomUnif, mean = mean, sd = sd, nu = nu, d = d, xi = xi)
    
    # Return Value:
    result
  }


#' GAt innovations
#'
#' Generalized asymmetric t distribution of Paolella (1997) with location 0,
#' scale 1 and parameters `nu` (tail), `d` (peakedness) and `xi` (asymmetry,
#' 1 is symmetric), see [dgat()]. Moments \eqn{E|z|^\delta} exist for
#' \eqn{\delta < \nu d}; the GARCH mode power is \eqn{\delta = 2}.
#'
#' @inheritParams gs_stable
#' @return A `gs_dist` object.
#' @examples
#' d <- gs_gat()
#' d$max_power(c(nu = 2, d = 4, xi = 1))
#' @export
gs_gat <- function(lower = c(nu = 0.05, d = 0.1, xi = 0.05),
                   upper = c(nu = 100, d = 50, xi = 20),
                   start = c(nu = 2, d = 4, xi = 1))
{
  gs_dist(
    name = "gat",
    log_density = function(z, par)
      dgat(z, nu = par[["nu"]], d = par[["d"]], xi = par[["xi"]], log = TRUE),
    random = function(n, par) rgat(n, nu = par[["nu"]], d = par[["d"]], xi = par[["xi"]]),
    cdf = function(q, par) pgat(q, nu = par[["nu"]], d = par[["d"]], xi = par[["xi"]]),
    quantile = function(p, par) qgat(p, nu = par[["nu"]], d = par[["d"]], xi = par[["xi"]]),
    par_names = c("nu", "d", "xi"),
    lower = lower, upper = upper, start = start,
    default_delta = 2,
    max_power = function(par) par[["nu"]] * par[["d"]],
    mean = function(par) {
      nu <- par[["nu"]]; d <- par[["d"]]; xi <- par[["xi"]]
      if (nu * d <= 1) return(NA_real_)
      (xi - 1 / xi) * nu^(1 / d) * beta(2 / d, nu - 1 / d) / beta(1 / d, nu)
    },
    aparch_moment = function(gamma, delta, par)
      .gat_aparch_moment(par[["nu"]], par[["d"]], par[["xi"]], gamma, delta),
    check = FALSE)
}
