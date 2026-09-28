\name{gsMomentAparch}
\alias{gsMomentAparch}


\title{

Computation of moments for several conditional distribution

}


\description{

Computation of the moments expression \eqn{E( |Z| - \gamma Z ) ^ \delta} where \eqn{\mathcal{D}_{\vartheta}(0,1)} is the density of the innovations with zero location and unit scale and \eqn{\vartheta} are additional distributional parameters that describe the skew and the shape of the distribution.

}


\usage{

gsMomentAparch(cond.dist = c("stableS1", "gev", "gat", "norm", 
"std", "sstd", "skstd", "ged"), shape = 1.5, skew = 0, delta = 1, gm = 0)

}


\arguments{

        \item{cond.dist}{ a character string naming conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"}, \code{"stableS1"}, \code{"stableS2"}, \code{"gev"} and \code{"gat"}. Other common distributions are also possible such as the \code{"norm"}, \code{"std"}, \code{"sstd"}, \code{"skstd"}, \code{"ged"}.
        }
        
  \item{shape, skew}{
        The shape and skew parameter for the conditional distribution. For the stable distribution the shape and skew are the index of stability (\eqn{\alpha}) and the asymmetry parameter (\eqn{\beta}). For the Generalized Asymmetric t distribution there are two shape parameters (\eqn{\nu} and \eqn{d}) and one asymmetry parameter called \eqn{\xi}. For the skew Student's t distribution the shape is the degrees of freedom (\eqn{\nu}) and \eqn{\xi} is the asymmetry parameter. Finally, for the GEV distribution the shape parameter is called \eqn{\xi}.
}

  \item{delta, gm}{
        The \eqn{\delta} and \eqn{\gamma} parameteres of the ARMA-GARCH/APARCH model.
}

See the R documentation of function \code{\link{gsFit}} for more details.

}


\details{

For each conditional distribution, the evaluation of the moment expression \eqn{E(|Z| - \gamma Z)^\delta} is crucial if one is interested in finding a stationary solution. Analytical expressions to compute these moments are available throughout the GARCH literature for many conditional distributions and we have derived an efficient formula for the asymmetric stable case.
It is true that we can use numerical integration to compute those moments but there are several drawbacks that represent limitations. The first one is that numerical integration is always slower than using a closed expression to compute those formulas. The second is that this technique fails when the distribution becomes very pick or very concentrated on a small portion of the real line. For example, if we use the R base function \code{integrate} to calculate the area total probability on the interval \eqn{(-\infty,\infty)} for the standard gat distribution with parameters \eqn{\nu = 2, d = 0.05} and \eqn{\xi = 1} we get 
\eqn{3.690067e-07}, which is not true since the total probability must be 1. Therefore, it is very important for us to have the exact expression for all those moments in order to evaluate correctly the stationarity restrictions. 

The \code{gsMomentAparch} function can be used to compute those expressions for the following distributions: Normal, Student's t, skewed Student's t from Fernandez and Steel (1998), GED, stable in its 1-parameterization and the GAt distribution. Only for the GEV distribution we use numerical integration.

Note: The \pkg{GEVStableGarch} package implements a slighty different version of the Student's skew distribution that is different from the one implemented inside the \pkg{fGarch} package. Therefore, the APARCH moment formula is only valid for the \pkg{GEVStableGarch} implementation (\code{skstd}) and not for the \pkg{fGarch} implementation \code{sstd}.
}


\value{

Returns the following expression for several conditional distributions
\eqn{E( |Z| - \gamma Z) ^ \delta}. If any of the input parameters are outside the parameter space, this function returns \eqn{\infty}.

}


\references{

    Ding, Z., Granger, C., Engle, R.F. (1993).
    \emph{A Long Memory Property of Stock Market Returns and a New Model}.
    Journal of Empirical Finance, 1, 83--106.

    Diongue, A.K. (2008).
    \emph{An investigation of Stable-Paretian Asymmetric Power GARCH Model}.
    Journal des Sciences, 8(4), 15--26.

    Lambert, P., Laurent, S. (2001).
    \emph{Modelling Financial Time Series Using GARCH-type Models with a Skewed Student Distribution for the Innovations}.
    Institut de Statistique, Universite Catholique de Louvain, Discussion Paper 0125.

    Mittnik, S., Paolella, M.S. (2000).
    \emph{Conditional Density and Value-At-Risk Prediction of Asian Currency Exchange Rates}. 
    Journal of Forecasting, 19(4), 313--333.

    Mittnik, S., Paolella, M.S., Rachev, S.T. (2002).
    \emph{Stationarity of stable power-GARCH processes}.
    Journal of Econometrics 106, 97--107.
    
}


\author{

Thiago do Rego Sousa.

}


\examples{

# Computation of the Moment E( |Z| - gamma Z) ^ delta for several distributions

gsMomentAparch(cond.dist = "stableS1", shape = 1.1, skew = 0, delta = 1.01, gm = 0.99999)

gsMomentAparch(cond.dist = "gev", shape = -4, skew = 0, delta = 1.4, gm = 0)

gsMomentAparch(cond.dist = "gat", shape = c(1.9,2.3), skew = 0.5, delta = 0.4, gm = 0)

gsMomentAparch(cond.dist = "norm", shape = c(1.9,2.3), skew =1, delta = 11.4, gm = -0.999)

gsMomentAparch(cond.dist = "std", shape = 2.001, skew = -0.5, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "sstd", shape = 2.001, skew = 0.11, delta = 2, gm = -0.99)

gsMomentAparch(cond.dist = "skstd", shape = 5.001, skew = 0.11, delta = 3, gm = -0.5)

gsMomentAparch(cond.dist = "ged", shape = 6, skew = 0.11, delta = 5.11, gm = -0.5)

}


\keyword{ aparch }
\keyword{ moments }