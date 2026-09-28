\name{GEVSTABLEGARCHSPEC-class}


\Rdversion{1.1}


\docType{class}


\alias{GEVSTABLEGARCHSPEC-class}


\alias{show,GEVSTABLEGARCHSPEC-update}


\title{Class \code{"GEVSTABLEGARCHSPEC"}}


\description{

This class represents an ARMA-GARCH/APARCH model with innovations following an univariate distribution with zero location and unit scale.

}


\section{Objects from the Class}{

Objects can be created by calls of the function \code{gsSpec}. This object 
specifies an univariate ARMA-GARCH/APARCH model with GEV or stable conditional
distribution. Note that we only allow the user
to create an object that do not violate the model specifications. For more details about
the definition of those models see:
see Mittnik and Paolella (2000) (GAt innovations), Mittnik et al. (2002) (stable innovations), Zhao et al. (2011) (GEV innovations) Wuertz et al. (2009) (finite variance innovations).

}


\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"call"}: 
        the call of the \code{\link{gsSpec}} function.
        } 
    \item{\code{formula}:}{Object of class \code{"formula"}: 
        a list with two formula entries for the mean and variance 
        equation of the combined ARMA-GARCH/APARCH model.
        }
    \item{\code{model}:}{Object of class \code{"list"}: 
        a list with the model parameters specifying the model parameters (see the model definition described in function \code{\link{gsFit}}.
        }
    \item{\code{presample}:}{Object of class \code{"matrix"}: 
        a numeric matrix with presample values.
        }
    \item{\code{distribution}:}{Object of class \code{"character"}: 
        a character string with the name of the conditional distribution.
        }
    \item{\code{rseed}:}{Object of class \code{"numeric"}: 
        an integer with the random number generator seed.
        }
  }
}


\section{Methods}{
  \describe{
  
    \item{show}{\code{signature(object = "GEVSTABLEGARCHSPEC")}: 
    prints an object of class 'GEVSTABLEGARCH' 
    }
    
	 }
}


\references{


    Mittnik, S., Paolella, M.S. (2000).
   \emph{Prediction of Financial Downside-Risk with Heavy-Tailed Conditional Distributions}
    Available at SSRN 391261.
    
    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

    Wuertz, D., Chalabi, Y., Luksan, L. (2009).
    \emph{Parameter Estimation of ARMA Models with GARCH/ APARCH Errors: An R and SPlus SoftwareImplementation}.   
    Journal of Statistical Software, forthcoming,
    \href{http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf}{http://www-stat.wharton.upenn.edu/~steele/...WurtzEtAlGarch.pdf}.

}


\author{

Thiago do Rego Sousa for the latest modifications \cr
Diethelm Wuertz for the original implementation of the fGARCHSPEC-class from package fGarch

}


\keyword{classes}
