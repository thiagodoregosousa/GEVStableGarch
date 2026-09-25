\name{gsSelect}
\alias{gsSelect}


\title{

Selects the best model according to goodness-of-fit criteria

}


\description{

This function estimates ARMA-GARCH/APARCH models with varying order and returns the one with the smallest goodness-of-fit criteria

}


\usage{
gsSelect(data, order.max = c(1, 1, 1, 1), 
selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, 
cond.dist = c("stableS0", "stableS1", "stableS2", "gev", 
"gat", "norm", "std", "sstd", "skstd", "ged"), 
include.mean = TRUE, 
algorithm = c("sqp", "sqp.restriction", "nlminb", "nlminb+nm"), ...)
}


\arguments{

      \item{data}{
      Data with the time series to be estimated. It must be a numeric vector not contain NA, NULL or
      Inf values.
        }
        
      \item{order.max}{
      Maximum order of models to search. It must by a vector of the type \cr
      \code{c(mMax, nMax, pMax, qMax)} with each entry representing the maximum order of the model that will be fitted when searching for the best model.
        }
        
      \item{selection.criteria}{
      The goodness-of-fit criterion to be used when searching for the best model. Three different criterias are allowed: AIC, AICc or BIC. For more information see Brockwell and Davis (1996).
        }
        
      \item{is.aparch}{
      Boolean variable indicating whether to search for ARMA-GARCH or ARMA-APARCH models.
        }
        
\item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"GAt"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}.
        }
        
      \item{include.mean}{
    This is a boolean variable. It intercept is TRUE than we estimate the model with intercept.
        }
        
        \item{algorithm}{
      The algorithm to be used to search for the optimum value. The current version of the \pkg{GEVStableGarch} package implements
four different optimization procedures, namely the \code{"sqp"}, \code{"sqp.restriction"} (enforce stationarity) , \code{"nlminb"} and \code{"nlminb+nm"}. See the details section for more information about the estimation algorithms.
      }
        
      \item{\dots}{
      Additional parameters to be passed to function \code{\link{gsFit}} during estimation.
        }
}


\value{
    Returns a S4 object of class "GEVSTABLEGARCH" with the best model. See \code{\link{GEVSTABLEGARCH-class}} for details. 
}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

}


\author{

Thiago do Rego Sousa.

}


\examples{

# Best ARMA-GARCH model within the range ARMA(0,0)-GARCH(1,0) to ARMA(0,0)-GARCH(1,1)
# using the Corrected Akaike Information Criteria (AICc)
data(dem2gbp)
x = dem2gbp[,1]
model = gsSelect (data = x, order.max = c(0,0,1,1), is.aparch = FALSE, 
          algorithm = "sqp", cond.dist = "gev", selection.criteria = "AIC", 
          include.mean = FALSE)

}