\name{Datasets}

\alias{Datasets}

\alias{dem2gbp}
\alias{sp500dge}



\title{Time Series Data Sets} 



\description{

    Data sets used in the examples of the \pkg{fGarch} package. 
    
}



\section{dem2gbp}{

The vector \code{dem2gbp} contains daily observations of the
  Deutschmark vs British Pound foreign exchange rate log-returns. This data
  set has been promoted as an informal benchmark for GARCH time-series software validation. See McCullough and Renfro (1999), and Brooks, Burke, and Persand (2001) for details.
  The nominal returns are expressed in percent as in Bollerslev and Ghysels (1996). The sample period is from January 3, 1984, to December 31, 1991, for a total of 1974 observations.
  
}



\section{sp500dge}{

The vector \code{sp500dge} contains daily returns from the SP500 index as used in the paper
of Ding, Granger and Engle (1993). The sample period is from January 3, 1928 to August 30, 1991 for a total of 17055 observations.
  
}



\references{

    Bollerslev T., Ghysels, E. (1996) 
    \emph{Journal of Business and Economic Statistics}.
    Periodic Autoregressive Conditional Heteroscedasticity, 14(2), 139--151

    Brooks C., Burke S.P., Persand G. (2001) 
    \emph{International Journal of Forecasting}.
    Benchmarks and the Accuracy of GARCH Model Estimation, 17(1), 45--57 

    Ding, Z., Granger, C., Engle, R.F. (1993).
    \emph{A Long Memory Property of Stock Market Returns and a New Model}.
    Journal of Empirical Finance, 1, 83--106.

    McCullough B. D., Renfro C. G. (1999)
    \emph{Journal of Economic and Social Measurement}
    Benchmarks and Software Standards: A Case Study of GARCH Procedures, 25 (2), 59--71
    
}