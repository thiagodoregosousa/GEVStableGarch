\name{gsSim}
\alias{gsSim}
\title{
Simulation of ARMA-GARCH/APARCH process
}


\description{

This functions simulates time series following ARMA-GARCH/APARCH models with GEV and stable distributions. This function was developed through the available code from \pkg{fGarch} package available at CRAN.

}


\usage{

gsSim(spec = gsSpec(), n = 100, n.start = 100)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{spec}{
A model specified with function \code{\link{gsSpec}}.
}
  \item{n}{
The size of simulated time series.
}
  \item{n.start}{
  Length of the "burn-in" period of the simulated time series.
}

}


\details{

The initial values of the time series
are fixed and the recursion formulas of the model are used to simulate the dynamics of the
process. We do not verify the stationarity conditions of the model because the simulation of
non-stationary process could also be of interest by the researcher.

}


\value{
  The function returns an object containing the following items:
	\item{ model }{
    A string describing the estimated model.
  }
	\item{ cond.dist }{
    The conditional distribution used to fit the model.
  }
	\item{ series}{
    An array with three columns, where the first column contains
    the simulated ARMA-GARCH/APARCH process (\code{garch}), the second column the conditional
    standard deviations (\code{"sigma"}), and the last column the innovations (\code{"eps"}).
  }
}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

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
Diethelm Wuertz for the original implementation of the garchSpec function from package fGarch

}


\examples{
# Simulation of a ARMA-APARCH process with stable conditional distribution
#x <- GSgarch.Sim(N = 2500, mu = 0.1,a = c(0.2,0.3),b = c(0.2,0.5),
#omega = 0.1, alpha = c(0.1,0.2),beta = c(0.1,0.1),gm=c(0.3,-0.3),
#delta = 1,skew = 0.3,shape = 1.9, cond.dis = "stable")
}