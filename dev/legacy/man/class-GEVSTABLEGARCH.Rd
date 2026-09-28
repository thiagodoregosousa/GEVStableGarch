\name{GEVSTABLEGARCH-class}


\Rdversion{1.1}


\docType{class}


\alias{GEVSTABLEGARCH-class}


\alias{show,GEVSTABLEGARCH-update}


\title{Class \code{"GEVSTABLEGARCH"}}


\description{

The class GEVSTABLEGARCH represents an ARMA-GARCH/APARCH 
model estimated wih function \code{\link{gsFit}}

}


\section{Objects from the Class}{

Objects can be created by calling function \code{\link{gsFit}}.
This object contain the estimated parameters of a time series process.

}


\section{Slots}{
  \describe{
  
    \item{\code{call}:}{Object of class \code{"call"}: 
        the call of the \code{\link{gsFit}} function.
        } 
        
    \item{\code{formula}:}{Object of class \code{"formula"}: 
        a formula object specifying mean and variance equation.
        }
        
    \item{\code{method}:}{Object of class \code{"character"}: 
        a string describing the optimization method used to search for the optimum value
        }
    
    \item{\code{convergence}:}{Object of class \code{"numeric"}:
        an integer code. 0 indicates successful convergence of of the estimation method used to perform the optimization of the log-likelihood function. A value different from zero indicates a failure in achieving convergence. Notice that sometimes
the optimization algorithm will return a "true" convergence, even 
when the optimized negative log-likelihood equals to 1e99. In this 
case, we set the variable \code{convergence} to 1 to indicate that 
convergence was not achieved. In the general case, the variable \code{convergence} assumes the same 
value reported by the internal solver used inside function \code{gsFit}. Hence, aditional interpretation of the \code{convergence} codes can be made by using the R help of the corresponding optimization routine: \code{solnp} (\code{"sqp"} and \code{"sqp.restriction"} algorithms) or \cr \code{nlminb} (\code{"nlminb"} and \code{"nlminb+nm"} algorithms). 
        }
        
    \item{\code{messages}:}{Object of class \code{"list"}:
        a character string giving additional informations collected
        during estimation.
        }
    
   \item{\code{data}:}{Object of class \code{"numeric"}: 
        a numeric vector containing the data of
        the time series to be estimated.
        }
    
        \item{\code{fit}:}{Object of class \code{"list"}: 
        a list with the results from the parameter estimation.
        }
    
        \item{\code{residuals}:}{Object of class \code{"numeric"}: 
        a numeric vector with the residual values.
        }
    
    \item{\code{h.t}:}{Object of class \code{"numeric"}: 
        a numeric vector with the conditional variances.
        }
        
    \item{\code{sigma.t}:}{Object of class \code{"numeric"}: 
        a numeric vector with the conditional standard deviations.
        }
        
    \item{\code{title}:}{Object of class \code{"character"}: 
        a string with the title.
        }
        
    \item{\code{description}:}{Object of class \code{"character"}: 
        a string with a description.
        }

  }
}
\section{Methods}{

  \describe{
  
      \item{show}{\code{signature(object = "GEVSTABLEGARCH")}: 
        prints an object of class 'GEVSTABLEGARCH'.
        }  
	 }
}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

}


\author{

Thiago do Rego Sousa for the latest modifications \cr
Diethelm Wuertz for the original implementation of the fGARCH-class from package fGarch

}




\keyword{classes}