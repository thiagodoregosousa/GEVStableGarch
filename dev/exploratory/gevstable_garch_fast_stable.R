# Manually load files from GEVStableGarch changing 
#.GSgarch.dstable for the appropriate functions from libstable4u in 
# https://cran.r-project.org/web/packages/libstable4u/index.html

# I found that these are the same
# pars = c(alpha, beta, gamma, delta)
# stable_pdf(x, pars) AND  stabledist::dstable(x, alpha, beta, gamma, delta, pm = 0)

pacman::p_load(GEVStableGarch, numDeriv, libstable4u)

.GSgarch.dstable = function (x, alpha = 1.5, beta = 0, gamma = 1, delta = 0, param = 0) 
{
  return(stable_pdf(x = x, pars = c(alpha, beta, gamma, delta), parametrization = param))
}



.armaGarchDist = function (z, hh, shape = 1.5, skew = 0, cond.dist = c("stableS0", 
                                                      "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                      "skstd", "ged"), TOLG = 1e-08) 
{
  cond.dist = match.arg(cond.dist)
  if (length(z) != length(hh)) 
    stop("Error: Vectors 'z' and 'hh' have different length.")
  if (sum(is.na(hh)) > 0 || min(hh) == 0) {
    return(1e+99)
  }
  if (cond.dist == "norm") 
    return(-sum(log(dnorm(x = z/hh)/hh)))
  if (cond.dist == "std") {
    if (!(shape > 2)) 
      return(Inf)
    return(-sum(log(dstd(x = z/hh, nu = shape)/hh)))
  }
  if (cond.dist == "sstd") {
    if (!(shape > 2) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dsstd(x = z/hh, nu = shape, xi = skew)/hh)))
  }
  if (cond.dist == "skstd") {
    if (!(shape > 2) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dskstd(x = z/hh, nu = shape, xi = skew)/hh)))
  }
  if (cond.dist == "gat") {
    if (!(shape[1] > 0) || !(shape[2] > 0) || !(skew > 0)) 
      return(1e+99)
    return(-sum(log(dgat(x = z/hh, nu = shape[1], d = shape[2], 
                         xi = skew)/hh)))
  }
  if (cond.dist == "ged") {
    if (!(shape > 0)) 
      return(1e+99)
    return(-sum(log(dged(x = z/hh, nu = shape)/hh)))
  }
  if (cond.dist == "gev") {
    if ((shape[1] <= -0.5) || (shape[1] >= 0.5)) 
      return(1e+99)
    if (abs(shape[1]) < 1e-07) {
      llh <- sum(log(hh)) + sum(z/hh) + sum(exp(-z/hh))
      return(llh)
    }
    sig <- hh
    y <- 1 + shape * (z/hh)
    llh <- sum(log(sig)) + sum(y^(-1/shape)) + sum(log(y)) * 
      (1/shape + 1)
    return(llh)
  }
  if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) {
    if (!(shape > 1) || !(shape < 2) || !(abs(skew) < 1)) 
      return(1e+99)
    y <- z/hh
    if (sum(is.na(y)) || sum(is.nan(y)) || sum(is.infinite(y))) 
      return(1e+99)
    if (cond.dist == "stableS0") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 0)/hh))
    if (cond.dist == "stableS1") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 1)/hh))
    if (cond.dist == "stableS2") 
      result = -sum(log(.GSgarch.dstable(x = z/hh, alpha = shape, 
                                         beta = skew, param = 2)/hh))
    return(result)
  }
}



gsFit = function (formula = ~garch(1, 1), data, cond.dist = c("stableS0", 
                                                      "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                      "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
                                                                                                          "sqp.restriction", "nlminb", "nlminb+nm"), control = NULL, 
          tolerance = NULL, title = NULL, description = NULL) 
{
  DEBUG = FALSE
  cond.dist = match.arg(cond.dist)
  algorithm = match.arg(algorithm)
  if (!is.numeric(data) || any(is.na(data)) || any(is.null(data)) || 
      any(!is.finite(data))) 
    stop("Invalid 'data' input. It may be contain NA, NULL or Inf.")
  CALL = match.call()
  if (is.null(tolerance)) 
    tolerance = list(TOLG = 1e-08, TOLSTABLE = 0.01, TOLSTATIONARITY = 0.001)
  if (tolerance$TOLSTATIONARITY > 0.05) 
    stop("TOLSTATIONARITY can not be > 0.05")
  if (is.null(title)) 
    title = "ARMA-GARCH modelling"
  formula.input <- formula
  formula <- .getFormula(formula)
  m <- formula$formula.order[1]
  n <- formula$formula.order[2]
  p <- formula$formula.order[3]
  q <- formula$formula.order[4]
  APARCH <- formula$isAPARCH
  formula.mean <- ""
  formula.var <- ""
  if (m > 0 || n > 0) 
    formula.mean <- formula$formula.mean
  else formula.mean <- ""
  if (p > 0) 
    formula.var <- formula$formula.var
  else formula.var <- ""
  if (algorithm == "sqp.restriction") {
    if (formula.var == "") 
      stop("sqp.restriction should only be used with GARCH/APARCH models.")
    if (any(cond.dist == c("stableS0", "stableS2"))) 
      stop("sqp.restriction algorithm can only be used with \n stable distribution in S1 parametrization, i.e., cond.dist = 'stableS1'.")
  }
  if (formula.var == "") 
    stop("Pure ARMA model not allowed")
  printRes = TRUE
  ARMAonly <- FALSE
  AR <- FALSE
  MA <- FALSE
  GARCH <- FALSE
  if (m == 0) 
    AR <- TRUE
  if (n == 0) 
    MA <- TRUE
  if (q == 0) 
    GARCH <- TRUE
  optim.finished <- FALSE
  if (AR == TRUE) 
    m <- 1
  if (MA == TRUE) 
    n <- 1
  if ((p == 0) && (q == 0)) 
    ARMAonly = TRUE
  if (GARCH == TRUE && !ARMAonly) 
    q <- 1
  data <- data
  N <- length(data)
  out <- NULL
  messages <- NULL
  out$order <- c(m, n, p, q, include.mean, APARCH)
  TMPvector <- c(if (m != 0 || n != 0) c("arma(", m, ",", n, 
                                         ")-"), if (APARCH == FALSE) c("garch(", p, ",", q, ")"), 
                 if (APARCH == TRUE) c("aparch(", p, ",", q, ")"))
  TMPorder <- paste(TMPvector, collapse = "")
  TMPvectorintercept <- c("include.mean:", if (include.mean == 
                                               TRUE) "TRUE", if (include.mean == FALSE) "FALSE")
  TMPintercept <- paste(TMPvectorintercept, collapse = "")
  out$model <- paste(TMPorder, "##", TMPintercept, collapse = "")
  out$cond.dist <- cond.dist
  out$data <- data
  optim.finished <- FALSE
  if (cond.dist == "gat") 
    lengthShape = 2
  else lengthShape = 1
  if (DEBUG) 
    print(c("lengthShape", lengthShape))
  arCheck <- function(ar) {
    p <- max(which(c(1, -ar) != 0)) - 1
    if (!p) 
      return(TRUE)
    all(Mod(polyroot(c(1, -ar[1L:p]))) > 1)
  }
  garchLLH = function(parm) {
    if (sum(is.nan(parm)) != 0) {
      return(1e+99)
    }
    mu <- parm[1]
    a <- parm[(1 + 1):(2 + m - 1)]
    b <- parm[(1 + m + 1):(2 + m + n - 1)]
    omega <- parm[1 + m + n + 1]
    alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
    gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
    beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
                                            p + q - 1)]
    delta <- parm[2 + m + n + 2 * p + q + 1]
    skew <- parm[3 + m + n + 2 * p + q + 1]
    shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
                                                 2 * p + q + lengthShape)]
    if (!APARCH) {
      gm = rep(0, p)
      delta = 2
      if (any(cond.dist == c("stableS0", "stableS1", "stableS2"))) 
        delta = 1
    }
    if (AR == TRUE) 
      a <- 0
    if (MA == TRUE) 
      b <- 0
    if (GARCH == TRUE) 
      beta <- 0
    if (include.mean == FALSE) 
      mu <- 0
    cond.general <- FALSE
    cond.normal <- FALSE
    cond.student <- FALSE
    cond.gev <- FALSE
    cond.stable <- FALSE
    parset <- c(omega, alpha, if (!GARCH) beta, delta)
    cond.general <- any(parset < tolerance$TOLG)
    if (cond.general || cond.student || cond.gev || cond.stable) {
      return(1e+99)
    }
    if (!arCheck(a)) 
      return(1e+99)
    if (AR == TRUE && MA == TRUE) {
      filteredSeries <- .filterAparch(data = data, p = p, 
                                      q = q, mu = mu, omega = omega, alpha = alpha, 
                                      beta = beta, gamma = gm, delta = delta)
      z <- filteredSeries[, 1]
      hh <- filteredSeries[, 2]
    }
    if (AR == FALSE || MA == FALSE) {
      filteredArma <- .filterArma(data = data, size.data = N, 
                                  m = m, n = n, mu = mu, a = a, b = b)
      filteredSeries <- .filterAparch(data = filteredArma, 
                                      p = p, q = q, mu = 0, omega = omega, alpha = alpha, 
                                      beta = beta, gamma = gm, delta = delta)
      z <- filteredSeries[, 1]
      hh <- filteredSeries[, 2]
    }
    if (optim.finished) {
      out$residuals <<- as.numeric(z)
      out$sigma.t <<- as.numeric(hh)
      out$h.t <<- as.numeric(hh^delta)
    }
    llh.dens <- .armaGarchDist(z = z, hh = hh, shape = shape, 
                               skew = skew, cond.dist = cond.dist)
    llh <- llh.dens
    if (is.nan(llh) || is.infinite(llh) || is.na(llh)) {
      llh <- 1e+99
    }
    llh
  }
  start <- .getStart(data = data, m = m, n = n, p = p, q = q, 
                     AR = AR, MA = MA, cond.dist = cond.dist, TOLG = tolerance$TOLG, 
                     TOLSTABLE = tolerance$TOLSTABLE)
  if (DEBUG) {
    print("start")
    print(start)
  }
  garch.stationarity <- function(parm) {
    omega <- parm[1 + m + n + 1]
    alpha <- parm[(2 + m + n + 1):(3 + m + n + p - 1)]
    gm <- parm[(2 + m + n + p + 1):(3 + m + n + p + p - 1)]
    beta <- parm[(2 + m + n + 2 * p + 1):(3 + m + n + 2 * 
                                            p + q - 1)]
    delta <- parm[2 + m + n + 2 * p + q + 1]
    skew <- parm[3 + m + n + 2 * p + q + 1]
    shape <- parm[(4 + m + n + 2 * p + q + 1):(4 + m + n + 
                                                 2 * p + q + lengthShape)]
    model = list(omega = omega, alpha = alpha, gm = gm, beta = beta, 
                 delta = delta, skew = skew, shape = shape)
    .stationarityAparch(model = model, formula = formula, 
                        cond.dist = cond.dist)
  }
  modelLLH <- garchLLH
  if (algorithm == "nlminb") {
    fit1 <- nlminb(start[1, ], objective = modelLLH, lower = start[2, 
    ], upper = start[3, ], control = control)
    out$llh <- fit1$objective
    out$par <- fit1$par
    out$hessian <- fit1$hessian
    out$convergence <- fit1$convergence
  }
  if (algorithm == "nlminb+nm") {
    fit1.partial <- nlminb(start[1, ], objective = modelLLH, 
                           lower = start[2, ], upper = start[3, ], control = control)
    fit1 <- optim(par = fit1.partial$par, fn = modelLLH, 
                  method = "Nelder-Mead", hessian = TRUE)
    out$llh <- fit1$value
    out$par <- fit1$par
    out$hessian <- fit1$hessian
    out$convergence <- fit1$convergence
  }
  if (algorithm == "sqp") {
    fit1 <- solnp(pars = start[1, ], fun = modelLLH, LB = start[2, 
    ], UB = start[3, ], control = control)
    out$llh <- fit1$values[length(fit1$values)]
    out$par <- fit1$pars
    out$hessian <- fit1$hessian
    
    
    # # compute hessian using 
    # out$hessian <- numDeriv::hessian(
    #   func = modelLLH,
    #   x = fit1$pars,
    #   method = "Richardson",
    #   method.args = list(
    #     eps          = 1e-4,
    #     d            = 0.1,
    #     zero.tol     = sqrt(.Machine$double.eps / 7e-7),
    #     r            = 4,
    #     v            = 2,
    #     show.details = FALSE
    #   )
    # )
    # print("Two alternatives to compute hessian, one returned from optim and another from package numDeriv")
    # print(fit1$hessian)
    # print(out$hessian)
    # print(sqrt(diag(solve(fit1$hessian))))
    # print(sqrt(diag(solve(out$hessian))))
    
    
    
    
    out$convergence <- fit1$convergence
  }
  if (algorithm == "sqp.restriction") {
    fit1 <- solnp(pars = start[1, ], fun = modelLLH, ineqfun = garch.stationarity, 
                  ineqLB = 0, ineqUB = 1 - tolerance$TOLSTATIONARITY, 
                  LB = start[2, ], UB = start[3, ], control = control)
    out$llh <- fit1$values[length(fit1$values)]
    out$par <- fit1$pars
    out$hessian <- fit1$hessian
    sizeHessian = length(out$hessian[1, ]) - 1
    out$hessian = out$hessian[1:sizeHessian, 1:sizeHessian]
    out$convergence <- fit1$convergence
  }
  if ((any(c("sqp", "nlminb") == algorithm)) && !is.numeric(try(sqrt(diag(solve(out$hessian))), 
                                                                silent = TRUE))) {
    fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                          method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
    if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
      out$hessian = fit1.partial$hessian
      if (DEBUG) 
        print(abs(fit1.partial$par - fit1$par))
    }
  }
  if ((algorithm == "sqp.restriction")) {
    fit1.partial <- optim(par = fit1$par, fn = modelLLH, 
                          method = "Nelder-Mead", hessian = TRUE, control = list(maxit = 1))
    
  
    
    
    if (sum(abs(fit1.partial$par - fit1$par)) == 0) {
      out$hessian = fit1.partial$hessian
      if (DEBUG) 
        print(abs(fit1.partial$par - fit1$par))
    }
  }
  if (DEBUG) 
    print(fit1)
  if (out$llh == 1e+99) 
    out$convergence = 1
  optim.finished = TRUE
  modelLLH(out$par)
  if (!ARMAonly) {
    outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
                                                           1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
                                                                                                           m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
                                                                                                                                                                         m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                        m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
                                                                                                                                                                                                                                                                               m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                                                                                                                                              m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
                                                                                                                                                                                                                                                                                                                                                                                "stableS0", "stableS1", "stableS2", "gat") == cond.dist)) (3 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                             m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "stableS1", "stableS2", "sstd", "skstd", "ged", "gat") == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             cond.dist)) (4 + m + n + 2 * p + q + 1):(4 + m + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n + 2 * p + q + lengthShape))
  }
  else {
    outindex <- c(if (include.mean) 1, if (AR == FALSE) (1 + 
                                                           1):(2 + m - 1), if (MA == FALSE) (1 + m + 1):(2 + 
                                                                                                           m + n - 1), if (!ARMAonly) (1 + m + n + 1), if (!ARMAonly) (2 + 
                                                                                                                                                                         m + n + 1):(3 + m + n + p - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                        m + n + p + 1):(3 + m + n + p + p - 1), if (!GARCH) (2 + 
                                                                                                                                                                                                                                                                               m + n + 2 * p + 1):(3 + m + n + 2 * p + q - 1), if (APARCH) (2 + 
                                                                                                                                                                                                                                                                                                                                              m + n + 2 * p + q + 1), if (any(c("sstd", "skstd", 
                                                                                                                                                                                                                                                                                                                                                                                "stableS0", "stableS1", "stableS2", "gat") == cond.dist)) (1 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                             m + n + 2 * p + q + 1), if (any(c("std", "gev", "stableS0", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               "stableS1", "stableS2", "sstd", "skstd", "ged", "gat") == 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             cond.dist)) (2 + m + n + 2 * p + q + 1):(2 + m + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        n + 2 * p + q + lengthShape), length(out$par))
  }
  outnames <- c(if (include.mean) "mu", if (AR == FALSE) paste("ar", 
                                                               1:m, sep = ""), if (MA == FALSE) paste("ma", 1:n, sep = ""), 
                if (!ARMAonly) "omega", if (!ARMAonly) paste("alpha", 
                                                             1:p, sep = ""), if (APARCH) paste("gamma", 1:p, sep = ""), 
                if (!GARCH) paste("beta", 1:q, sep = ""), if (APARCH) "delta", 
                if (any(c("sstd", "stableS0", "stableS1", "stableS2", 
                          "gat", "skstd") == cond.dist)) "skew", if (any(c("std", 
                                                                           "gev", "stableS0", "stableS1", "stableS2", "sstd", 
                                                                           "ged", "gat", "skstd") == cond.dist)) paste("shape", 
                                                                                                                       1:lengthShape, sep = ""), if (ARMAonly) "sigma")
  if (DEBUG) {
    print(c("out", out))
    print(c("outindex", outindex))
  }
  out$par <- out$par[outindex]
  names(out$par) <- outnames
  out$hessian <- out$hessian[outindex, outindex]
  nParam <- length(out$par)
  out$aic = 2 * out$llh + 2 * nParam
  out$aicc = 2 * out$llh + 2 * (nParam + 1) * N/(N - nParam - 
                                                   2)
  out$bic = 2 * out$llh + nParam * log(N)
  out$ics = c(out$aic, out$bic, out$aicc)
  names(out$ics) <- c("AIC", "BIC", "AICc")
  if (printRes) {
    solveHessianFailed = FALSE
    out$se.coef <- 0
    out$se.coef <- try(sqrt(diag(solve(out$hessian))), silent = TRUE)
    if (!is.numeric(out$se.coef)) {
      solveHessianFailed = TRUE
      messages$solving.hessian.matrix = "Error inverting Hessian Matrix"
      out$matcoef <- cbind(out$par, rep(NA, length(out$par)), 
                           rep(NA, length(out$par)), rep(NA, length(out$par)))
      dimnames(out$matcoef) = dimnames(out$matcoef) = list(names(out$par), 
                                                           c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
    }
    else {
      out$tval <- try(out$par/out$se.coef, silent = TRUE)
      out$matcoef = cbind(out$par, if (is.numeric(out$se.coef)) 
        out$se.coef, if (is.numeric(out$tval)) 
          out$tval, if (is.numeric(out$tval)) 
            2 * (1 - pnorm(abs(out$tval))))
      dimnames(out$matcoef) = list(names(out$tval), c(" Estimate", 
                                                      " Std. Error", " t value", "Pr(>|t|)"))
    }
    cat("\nFinal Estimate of the Negative LLH:\n")
    cat("-LLH:", out$llh)
    if (out$convergence == 0) 
      messages$optimization.algorithm = "Algorithm achieved convergence"
    else messages$optimization.algorithm = "Algorithm did not achieved convergence"
    cat("\nCoefficient(s):\n")
    printCoefmat(round(out$matcoef, digits = 6), digits = 6, 
                 signif.stars = TRUE)
  }
  out$order <- c(formula$formula.order[1], formula$formula.order[2], 
                 formula$formula.order[3], formula$formula.order[4])
  names(out$order) <- c("m", "n", "p", "q")
  fit <- list(par = out$par, llh = out$llh, hessian = out$hessian, 
              ics = out$ics, order = out$order, cond.dist = cond.dist, 
              se.coef = out$se.coef, tval = out$tval, matcoef = out$matcoef)
  new("GEVSTABLEGARCH", call = as.call(match.call()), formula = formula.input, 
      method = "Max Log-Likelihood Estimation", convergence = out$convergence, 
      messages = messages, data = data, fit = fit, residuals = out$residuals, 
      h.t = out$h.t, sigma.t = as.vector(out$sigma.t), title = as.character(title), 
      description = as.character(description))
}


gsSelect = function (data, order.max = c(1, 1, 1, 1), selection.criteria = c("AIC", "AICc", "BIC"), is.aparch = FALSE, cond.dist = c("stableS0", 
                                                                                                                   "stableS1", "stableS2", "gev", "gat", "norm", "std", "sstd", 
                                                                                                                   "skstd", "ged"), include.mean = TRUE, algorithm = c("sqp", 
                                                                                                                                                                       "sqp.restriction", "nlminb", "nlminb+nm"), ...) 
{
  cond.dist = match.arg(cond.dist)
  algorithm = match.arg(algorithm)
  selection.criteria = match.arg(selection.criteria)
  if (!is.numeric(data) || !is.vector(data)) 
    stop("data set must be a numerical one dimensional vector")
  goodness.of.fit.min <- 1e+99
  goodness.of.fit.current <- 1e+99
  fit.min <- list()
  fit <- list()
  order.list <- .getOrder(order.max = order.max)
  order.list.size <- length(order.list[, 1])
  for (i in 1:order.list.size) {
    m = order.list[i, 1]
    n = order.list[i, 2]
    p = order.list[i, 3]
    q = order.list[i, 4]
    if (is.aparch == TRUE) {
      formula = as.formula(paste("~ arma(", m, ", ", n, 
                                 ") + aparch(", p, ", ", q, ")", sep = "", collapse = NULL))
    }
    else {
      formula = as.formula(paste("~ arma(", m, ", ", n, 
                                 ") + garch(", p, ", ", q, ")", sep = "", collapse = NULL))
    }
    cat("\n------------------------------------------------------------------------------------------\n")
    cat(paste("Model: ", paste(formula)[2], " with '", cond.dist, 
              "' conditional distribution", sep = ""))
    cat("\n------------------------------------------------------------------------------------------\n")
    fit = gsFit(data = data, formula = formula, cond.dist = cond.dist, 
                algorithm = algorithm, include.mean = include.mean, 
                ...)
    goodness.of.fit.current = fit@fit$ics[selection.criteria]
    if (goodness.of.fit.current < goodness.of.fit.min) {
      fit.min <- fit
      goodness.of.fit.min <- fit@fit$ics[selection.criteria]
    }
  }
  cat("\n------------------------------------------------------------------------------------------\n")
  cat(paste("Best Model: ", paste(fit.min@formula)[2]))
  cat("\n------------------------------------------------------------------------------------------\n")
  return(fit.min)
}
