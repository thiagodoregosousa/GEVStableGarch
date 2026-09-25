\name{show-methods}


\docType{methods}


\alias{show-methods}
\alias{show,GEVSTABLEGARCH-method}
\alias{show,GEVSTABLEGARCHSPEC-method}


\title{GEVSTABLEGARCH Package Show Methods}


\description{

Methods to organize the output to the user when printing objects of class GEVSTABLEGARCH and GEVSTABLEGARCHSPEC. The show methods and the classes GEVSTABLEGARCH and GEVSTABLEGARCHSPEC were developed through the available code from \pkg{fGarch} package available at CRAN in order to make the interfaces more similar to the user.

}

\section{Methods}{
\describe{

    \item{\code{signature(object = "GEVSTABLEGARCH")}}{
        Print function for objects of class \code{"GEVSTABLEGARCH"}
    }

    \item{\code{signature(object = "GEVSTABLEGARCHSPEC")}}{
        Print function for objects of class \cr \code{"GEVSTABLEGARCHSPEC"}
    }
}}


\author{

Thiago do Rego Sousa for the latest modifications \cr
Diethelm Wuertz for the original implementation of the show-methods from package fGarch

}


\references{

    Wuertz, D., Chalabi, Y., with contribution from Miklovic, M., Boudt, C., Chausse, P., and others (2013).
    \emph{fGarch: Rmetrics - Autoregressive Conditional Heteroskedastic Modelling, R package version 3010.82},   
    \url{http://CRAN.R-project.org/package=fGarch}.

}


\keyword{models}