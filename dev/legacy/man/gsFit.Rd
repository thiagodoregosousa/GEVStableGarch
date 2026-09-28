\name{gsFit}
\alias{gsFit}


\title{

Estimation of ARMA-GARCH/APARCH models

}


\description{

This function uses Maximum Likelihood technique to estimate the parameters of
ARMA-GARCH or ARMA-APARCH model with several conditional distributions. The user can also enforce stationarity during estimation

}


\usage{

gsFit(formula = ~garch(1, 1), data, cond.dist = c("stableS0", "stableS1", 
"stableS2", "gev", "gat", "norm", "std", "sstd", "skstd", "ged"), 
include.mean = TRUE, algorithm = c("sqp", "sqp.restriction", 
"nlminb", "nlminb+nm"), control = NULL, 
tolerance = NULL, title = NULL, description = NULL)

}


\arguments{

    \item{formula}{ 
        a formula object specifying mean and variance equation. The mean equation should be specified using the names: "arma"", "garch"" or "aparch". For example: ~arma(1,0)-garch(1,0) for AR(1)-ARCH(1), or ~arma(0,1)-garch(2,2) for MA(1)-GARCH(2,2). 
        }
        
   \item{data}{ 
        a numeric vector containing the data of
        the time series to be estimated.
        }

        \item{cond.dist}{ a character string naming the conditional distribution of innovations. The package was created to accept the following distributions: \code{"stableS0"} (stable in S0-parameterization), \code{"stableS1"} (stable in S1-parameterization), \code{"stableS2"} (stable in S2-parameterization), \code{"gev"}, \code{"gat"} and \code{"skstd"} (skew Student's t from Fernandez and Steel (1998)). Other common conditional distribution (normal, Student's t, \code{"sstd"} (skew Student's t from \pkg{fGarch} package) and GED) are also allowed since they are very important for testing purposes.
        }
        
        \item{include.mean}{
          This is a boolean variable. It intercept is TRUE than we estimate the model with intercept.}
      
        \item{algorithm}{
      The algorithm to be used to search for the optimum value. The current version of the \pkg{GEVStableGarch} package implements
four different optimization procedures, namely the \code{"sqp"}, \code{"sqp.restriction"} (enforce stationarity) , \code{"nlminb"} and \code{"nlminb+nm"}. See the details section for more information about the estimation algorithms.
      }

        \item{control}{
        Control parameters passed to the optimization routine. For the \code{"sqp"} and \cr \code{"sqp.restriction"} algorithms the control parameters are passed to the \code{solnp} function from \pkg{Rsolnp} package. For the other algorithms (\code{"nlminb"} and \cr \code{"nlminb+nm"}) the control parameters are passed to the \code{nlminb} function from \pkg{stats} package.
        }

      \item{tolerance}{
        Tolerance for parameter estimation. When specified, the tolerance parameter should be passed as a list in the following format \cr \code{list( TOLG = 1e-8, TOLSTABLE = 1e-2, TOLSTATIONARITY = 1e-3 )}, where \code{TOLG} is the tolerance used to set the boundary region of the model parameters, \code{TOLSTABLE} is the boundary of parameters for stable distribution and \cr \code{TOLSTATIONARITY} is the tolerance parameter used by the \code{"sqp.restriction"} algorithm when searching for the stationary solution (see the details). 
        }

    \item{title}{ 
        a string with the title.
        }
        
    \item{description}{
        a string with a description.
        }
}


\details{

The starting values of the model parameters are key in getting any possible convergence and they were chosen to reflect this. For example, the GEV shape starting value was chosen to be 0 since in this case the support of the distribution with zero location and unit scale will encompasses both positive and negative parts of the real line (the support is at (-\eqn{\infty},+\eqn{\infty})). 

The parameters will be interpreted according to the following 
equations (see Wurtz et al. ,2009)
  \deqn{				
				X_t =  \mu + \sum_{i=1}^m a_i X_{t-i} + \sum_{j=1}^n b_j \varepsilon_{t-j} + \varepsilon_t 
		}
	\deqn{	
				  \varepsilon_t = \sigma_t z_t,\;\;\;z_t  \stackrel{iid}{\sim} D(0,1)\;,\\
		}
	\deqn{	
				\sigma_t^\delta   = \omega + \sum_{i=1}^p \alpha_i(\varepsilon_{t-i}-\gamma_i|\varepsilon_{t-i}|)^\delta + \sum_{j=1}^q \beta_j \sigma_{t-j}^\delta		
		}
where \eqn{\mathcal{D}_{\vartheta}(0,1)} is the density of the innovations with zero location and unit scale and \eqn{\vartheta} are additional distributional parameters that describe the skew and the shape of the distribution.

\subsection{Etimation Algorithms}{

Most software packages implement the estimation of GARCH models without imposing stationarity, but restricting the parameter set by appropriate bounds. This last approach was implemented in the \pkg{GEVStableGarch} package through the following algorithms: \code{"sqp"}, \code{"nlminb"} and \code{"nlminb + nm"}. The first two algorithms search for the optimum value by restricting the parameter set to appropriate lower and upper bounds. The last implements a two step optimization procedure, which consists in starting the search by using the constrained routine \code{nlminb} and then performing another search using an unconstrained method (in our case the Nelder-Mead method implemented in the \code{R} base function \code{optim}). This approach was suggested by Wuertz et al. (2009) since in many cases it leads to an improved solution (in terms of the likelihood function of the data).

Finally, the \code{"sqp.restriction"} algorithm performs a constrained search to maximize the log-likelihood function in order to obtain an stationary model. It can be used to estimate models with the following conditional distributions: \code{"stableS1"}, \code{"gev"}, \code{"gat"}, \code{"skstd"}, \code{"norm"}, \code{"std"} and \code{"ged"}. 

}

\subsection{Stationarity Solutions}{

Since we want to estimate the parameters that better adjust real data to ARMA-APARCH models, the notion of stationarity is crucial. Usually, even when the data set is non-stationary in appearance, we still are able to apply transformation techniques so that the resulting time series can be reasonably modeled as a stationary process (see Brockwell and Davis (1996)).

The APARCH(p,q) model with finite variance innovations has a unique \eqn{\delta}-order stationary solution if and only if (see Ling and McAller (2002))
\deqn{
    \sum_{i=1}^p E(|Z| - \gamma_i Z)^\delta \alpha_i +\sum_{j=1}^q\beta_j < 1,
}

The case of infinite variance has a somewhat different historical background. The first assumption made on stable distributions is that the index of stability \eqn{\alpha} must be greater than one, because in this case the innovations have finite first moment. The second assumption is that they must have a \eqn{\delta}-moment finite, which means that we must restrict our model to \eqn{1 < \delta < \alpha}. 

Diongue et al. (2008) showed that the APARCH(p,q) (all coefficients \eqn{\gamma_i = 0}) model has a strictly stationary solution if and only if
\deqn{
    \sum_{i=1}^p E(|Z| - \gamma_i Z)^\delta\alpha_i +\sum_{j=1}^q\beta_j < 1,
}

where \eqn{Z} has distribution \eqn{S(\alpha,\beta;1)} (stable in 1-parameterization). The estimation of stationary models is mainly dependent on the time taken to evaluate the expression for the asymmetric stable distribution \eqn{E(|Z| - \gamma_i Z)^\delta}.

The moment expression used to verify the stationarity condition is given by 

\deqn{
E(|Z| - \gamma Z)^\delta = \frac{ (1 - \gamma )^\delta\tilde{\sigma}^{\delta+1}\Gamma(\delta+1)\Gamma(-\frac{\delta}{\alpha}) }
    { \alpha\tilde{\sigma}\Gamma\bigg[ \bigg(\frac{1}{2}+\tilde{\beta}\frac{k(\alpha)}{2\alpha}\bigg)(-\delta)\bigg]
    \Gamma\bigg[ \frac{1}{2}- \tilde{\beta}\frac{k(\alpha)}{2\alpha} + \bigg(\frac{1}{2} + \tilde{\beta}\frac{k(\alpha)}{2\alpha}\bigg)(\delta+1)\bigg] }
}

\deqn{
+ \frac{ (1 + \gamma )^\delta\tilde{\sigma}^{\delta+1}\Gamma(\delta+1)\Gamma(-\frac{\delta}{\alpha}) }
    { \alpha\tilde{\sigma}\Gamma\bigg[ \bigg(\frac{1}{2}-\tilde{\beta}\frac{k(\alpha)}{2\alpha}\bigg)(-\delta)\bigg]
    \Gamma\bigg[ \frac{1}{2} + \tilde{\beta}\frac{k(\alpha)}{2\alpha} + \bigg(\frac{1}{2} - \tilde{\beta}\frac{k(\alpha)}{2\alpha}\bigg)(\delta+1)\bigg] },
}

where

\deqn{
\tilde{\sigma} = \bigg[1 + \beta^2 \tan^2(\frac{\alpha\pi}{2})\bigg]^{\frac{1}{2\alpha}},
}

\deqn{
    k(\alpha) = \left\{
    \begin{array}{rl}
    \alpha & , \alpha < 1, \\
    \alpha-2 & , \alpha > 1,
    \end{array} \right.
}

and

\deqn{
    \tilde{\beta} = \left\{
    \begin{array}{rl}
    \frac{2}{\pi\alpha}\arctan(\beta\tan(\frac{\alpha\pi}{2})) & , 0 < \alpha < 1, \\
    \frac{2}{\pi(\alpha-2)}\arctan(\beta\tan(\frac{\pi(\alpha-2)}{2})) & , 1 < \alpha < 2.
    \end{array} \right.
}

Notice that this expression is only valid for the stable distribution in 1-parameterization, and thus
cannot be used to the estimate models with stable distribution in 0 or 2 parameterizations.

The stationarity restriction implemented inside \pkg{GEVStableGarch} package imposes the follow restriction \cr \cr     
\eqn{
   \sum_{i=1}^p E(|Z| - \gamma_i Z)^\delta\alpha_i +\sum_{j=1}^q\beta_j < 1 -
}
\code{TOL.STATIONARITY}, \cr

where the variable TOL.STATIONARITY is set to 0.001 by default. For numerical reasons, the value of this variable cannot be greater than 0.05.

}

}


\value{
  
    \code{gsFit} \cr
    returns a S4 object of class \code{"GEVSTABLEGARCH"} with the following slots: 
  
    \item{@call}{
        the call of the \code{gsFit} function.
    }
    \item{@formula}{
        a list with two formula entries, one for the mean and the other 
        one for the variance equation.
    }
    \item{@method}{
        a string denoting the optimization method.
    }  
   \item{@data}{ 
        a numeric vector containing the data of the estimated time sereis.
        }
        
    \item{@fit}{ 
    a list with the results from the parameter estimation: \cr
      \code{par} - the estimated parameters; \cr
      \code{llh} - the estimated negative log-likelihood function; \cr
      \code{hessian} - the hessian matrix returned by the optimization algorithm; \cr
      \code{ics} - the value of the goodness-of-fit measures (AIC, AICc and BIC) (See Brockwell and Davis, 2002 for more details); \cr
      \code{order} - a list with the ARMA and GARCH/APARCH orders; \cr
      \code{cond.dis} - the conditional distribution; \cr
      \code{se.coef} - standard errors of the estimated parameters; \cr
      \code{tValue} - tValue of the estimated parameters; \cr
      \code{matcoef} - an organized matrix with the estimated parameters.
    }
    
    \item{@residuals}{ 
    a numeric vector with the residual values (\eqn{\varepsilon_t}).
    }

    \item{@h.t}{
        a numeric vector with the conditional variance
        (\eqn{h_t = \sigma_t^\delta}). 
    }
    
    \item{@sigma.t}{
        a numeric vector with the conditional standard deviation (\eqn{\sigma_t}).
    }
    
    \item{@title}{
        a string with the title.
    }
    
    \item{@description}{
        a string with a description.
    } 
  
    The entries of the @fit slot show the results from the
    optimization.

}


\section{Enhances}{

The estimation of ARMA-GARCH/APARCH models with conditional stable distribution is mainly dependent
on the time taken during the calculation of density points. Although the package can be used with the R \pkg{stabledist} package, the estimation of such models is only feasible if we use the fast implementation of stable densities using the R \pkg{stable} package designed by J.P.Nolan (see \url{http://academic2.american.edu/~jpnolan/stable/stable.html}). 

The \pkg{stable} package implements a faster computation of stable 
densities that are accurately enough to perform numerical optimization and is available at \url{http://www.robustanalysis.com}.

}


\references{

    Brockwell, P.J., Davis, R.A. (1996).
    \emph{Introduction to Time Series and Forecasting}.
    Springer, New York.

    Fernandez, C., Steel, M.F. (1998).
    \emph{On Bayesian Modeling of Fat Tails and Skewness}.
    Journal of the American Statistical Association, Taylor & Francis Group, 93(441), 359--371.

    Ling, S., McAller, M. (2002).
    \emph{Necessary and Suficient Moment Conditions for the GARCH(r,s) and Asymmetric Power GARCH(r,s) Models}. 
    Econometric Theory, 18(03), 722--729.

    Mittnik, S., Paolella, M.S., Rachev, S.T. (2002).
    \emph{Stationarity of stable power-GARCH processes}.
    Journal of Econometrics 106, 97--107.

    Nolan, J.P. (1997).
    \emph{Numerical calculations of stable densities and distribution functions}.
    Communications in Statistics - Stochastic Models 13, 759--774.

    Nolan, J.P. (1997). 
    \emph{Maximum likelihood estimation and diagnostics for stable distribution}. 
    In O. E.Barndorff-Nielsen, T. Mikosch, e S. I. Resnick (Eds.),
    Levy Processes: Theory andApplications,
    Boston, Birkhauser, 379--400.

    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

    Wuertz, D., Chalabi, Y., Luksan, L. (2009).
    \emph{Parameter Estimation of ARMA Models with GARCH/ APARCH Errors: An R and SPlus SoftwareImplementation}.   
    Journal of Statistical Software, forthcoming,
    \href{http://www-stat.wharton.upenn.edu/~steele/Courses/956/RResources/GarchAndR/WurtzEtAlGarch.pdf}{http://www-stat.wharton.upenn.edu/~steele/...WurtzEtAlGarch.pdf}.

    Zhao, X., Scarrott, C.J., Oxley, L., Reale, M. (2011).
    \emph{GARCH dependence in extreme value models with Bayesian inference}.
    Mathematics an Computers in Simulation, 81, Issue 7, 1430--1440.

    Yinyu, Ye. (1987).
    \emph{Interior Algorithms for Linear, Quadratic, and Linearly Constrained Non-Linear Programming}.
    Ph.D. Thesis, Department of EES,
    Stanford University.

}

\author{

Thiago do Rego Sousa, Cira Etheowalda Guevara Otiniano
and Silvia Regina Costa Lopes

}


\seealso{

\code{\link{gsSelect}}

}

\examples{
# This examples uses the dem2gbp dataset to estimate
# an ARMA(1,1)-GARCH(1,1) with GEV conditional distribution.
data(dem2gbp)
x = dem2gbp[, 1]
gev.model = gsFit(data = x , formula = ~garch(1,1), cond.dist = "gev")
}



