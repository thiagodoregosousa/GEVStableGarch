\name{GEVStableGarch-package}


\alias{GEVStableGarch-package}
\alias{GEVStableGarch}


\docType{package}


\title{

ARMA-GARCH/APARCH modelling with GEV and stable distributions

}


\description{

This package is designed to perform maximum likelihood estimation of ARMA-GARCH/APARCH models 
with generalized extreme distribution (GEV), stable, generalized asymmetric t (GAt) and skew Student's t from Fernandez and Steel (1998). The package also allows the researcher to restrict the search within the stationarity region (see the \code{\link{gsFit}} function). Other common conditional distribution (normal, Student's t and generalized error distribution (GED)) are also allowed since they are very important for testing purposes.

}


\details{

    \tabular{ll}{
      Package: \tab GEVStableGarch \cr
      Type: \tab Package \cr
      Version: \tab 1.1 \cr
      Date: \tab 2015-07-19 \cr
      License: \tab GPL(>=2) \cr
      Depends: \tab R(>= 2.15.0), fGarch, fExtremes, stabledist, skewt, Rsolnp \cr
    }


GARCH models have proven to be highly effective for analyzing financial data over the past decades. In particular, the combination of ARMA-GARCH models with stable and GEV distributions was successfully applied for forecasting volatility and for the measurement of Value at Risk (VAR).

Choosing the normal distribution as probability distribution for the innovations was a common choice in the beginning of the development of ARCH-type models. But recent research Nolan (1999), Mittnik et al. (2002), Mittnik and Paolella (2003), Curto et al. (2006), Frain (2009), Zhao et al. (2011) has shown that other distributions should be considered, specially because normal distribution can not account for fat tails and asymmetry found in real data.

This package contains functions for simulating and estimating ARMA-GARCH or ARMA-APARCH models using the maximum likelihood technique (MLE) under different assumptions: GEV, stable, GAt (also known as t3-distribution) and skew Student's t (Fernandez and Steel (1998)).

The current version of package \pkg{GEVStableGarch} has a new algorithm that allows the user to enforce stationarity during estimation. Aditionally, it contains functions for selecting the best model according to a predetermined goodness-of-fit criteria (see \code{\link{gsSelect}}).

}


\section{Time Series Simulation}{

  Contains functions to simulate ARMA-GARCH/APARCH processes with conditional GEV or stable
  distributions. 
  Note: These routines were adapted from functions \code{garchSpec} and \code{garchSim} available in \pkg{fGarch} package to make the interfaces more similar.
  
  Functions:
  
  \tabular{ll}{
    \code{gsSpec} \tab Specifies an univariate ARMA-GARCH/APARCH model, \cr \cr
    \code{gsSim} \tab Simulates an ARMA-GARCH/APARCH process. }
  
}


\section{Parameter Estimation}{

  Contains functions to fit the parameters of ARMA-GARCH/APARCH
  time series processes.
  
  Functions:
  
  \tabular{ll}{
    \code{gsFit} Fits the parameters of an ARMA-GARCH/APARCH process. \cr This function also
    provides an algorithm to enforce stationarity during estimation, \cr \cr
    \code{gsSelect} Selects the best model according to a goodness-of-fit
    criteria. } 
}


\section{Other Conditional Distribution Functions}{

  Contains functions to compute density, distribution, quantiles and generate random values using important conditional distributions used in the garch literature.
  
  Functions:
  
  \tabular{ll}{
    \code{[dpqr]skstd}  Skew Student's t distribution function from Fernandez and Steel (1998), \cr \cr
    \code{[dpqr]gat}  Generalized Asymmetric t distribution (GAt). \cr
    The GAt distribution was also refered in the literature as t3-distribution. \cr \cr
    \code{gsMomentAparch}  Computes APARCH Moments of the form \eqn{E( |Z| - \gamma Z ) ^ \delta} for several distributions \cr 
    using analytical expressions (except for the GEV case).
  } 
  
}



\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano 
and Silvia Regina Costa Lopes.

Maintainer: Thiago do Rego Sousa <thiagoestatistico@gmail.com>

}


\keyword{ stable distribution }
\keyword{ GEV distribution }
\keyword{ GAt distribution }
\keyword{ t3-distribution }
\keyword{ ARMA-GARCH }
\keyword{ ARMA-APARCH }