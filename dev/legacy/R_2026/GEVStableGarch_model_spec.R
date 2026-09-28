#' GEVStableGarch_model_spec
#' 
#' Specifices ARMA-GARCH/APARCH model 
#' 
#' It sets valid orders for ARMA, GARCH and if APARCH allowing for the user to 
#' set which conditional distribution to use and if the model will include a mean
#' 
#' @param m,n,p,q integers specifying model order ARMA(m,n)-GARCH/APARCH(p,q)
#' @param aparch logical indicating if model is APARCH (TRUE) or GARCH(FALSE)
#' @param include_mean logical indicating whether the model includes a mean
#' @param cond_dist Character. Innovation distribution; one of
#'   \code{"stable"}, \code{"GEV"}, or \code{"GAT"}.
#'
#' @return A list containing the model specification with components:
#' \describe{
#'   \item{spec}{A list with:}
#'   \itemize{
#'     \item{order:}{ A vector with components \code{m}, \code{n}, \code{p}, \code{q}}
#'     \item{aparch}{ Logical; APARCH included or not}
#'     \item{include_mean}{ Logical; whether mean term is included}
#'     \item{cond_dist}{ Character; innovation distribution}
#'   }
#' }
GEVStableGarch_model_spec = function(m = 0, n = 0, p = 1, q = 0, 
                                     aparch = FALSE,
                                     include_mean = FALSE, 
                                     cond_dist = "stable")
{
  # check model is valid
  if (!cond_dist %in% c("stable", "norm", "gev")) {
    stop("`cond_dist` must be one of: 'stable', 'norm', 'gev'.")
  }
  
  if (!cond_dist %in% c("stable", "norm", "gev")) {
    stop("`cond_dist` must be one of: 'stable', 'norm', 'gev'.")
  }
  
  # value constraints
  if (m < 0) stop("`m` must be >= 0.")
  if (n < 0) stop("`n` must be >= 0.")
  if (q < 0) stop("`q` must be >= 0.")
  if (p < 1) stop("`p` must be >= 1.")
  
  if (!is.logical(aparch)) stop("`aparch` must be TRUE/FALSE.")
  if (!is.logical(include_mean)) stop("`include_mean` must be TRUE/FALSE.")
  
  # set model as a list with empty parameters
  spec = list( 
    order = list(m = m, n = n, p = p, q = q), 
    aparch = aparch, 
    include_mean = include_mean, 
    cond_dist = cond_dist)
  class(spec) <- "GEVStableGarch_spec"
  
  # return 
  spec
}


#' Set model parameters for a GEV/Stable/GAT-GARCH specification
#'
#' This function attaches valid parameter values to an existing
#' model specification created by \code{GEVStableGarch_model_spec()}.
#' It performs consistency checks between the specification orders
#' and the supplied parameter vectors.
#'
#' @param spec A model specification list produced by
#'   \code{GEVStableGarch_model_spec()}.
#' @param ar Numeric vector of AR coefficients of length \code{m}.
#'   Must be supplied if \code{m > 0}, otherwise must be \code{NULL}.
#' @param ma Numeric vector of MA coefficients of length \code{n}.
#'   Must be supplied if \code{n > 0}, otherwise must be \code{NULL}.
#' @param omega Single positive numeric value for the GARCH constant term.
#' @param alpha Numeric vector of ARCH coefficients of length \code{p}. Should be at least one. 
#' @param beta Numeric vector of GARCH coefficients of length \code{q}.
#'   Must be supplied if \code{q > 0}, otherwise must be \code{NULL}.
#' @param aparch Single positive numeric APARCH parameter. Required only when
#'   \code{spec$aparch = TRUE}; must be omitted otherwise.
#' @param mean Single numeric conditional mean parameter. Required only when
#'   \code{spec$include_mean = TRUE}; must be omitted otherwise.
#' @param dist_params List of numeric distributional parameters associated with
#'   the chosen conditional distribution in \code{spec$cond_dist}. 
#'   Parameters name should match the choosen distribution names.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item dimension compatibility between \code{(m,n,p,q)} and supplied vectors ar, ma, alpha and beta. 
#'   \item positivity constraints for \code{omega} and \code{aparch}
#'   \item conditional inclusion of \code{mean} and \code{aparch}
#'   \item numeric type of all distribution parameters
#' }
#'
#' @return
#' An object of class \code{"GEVStableGarch_model"} containing:
#' \itemize{
#'   \item \code{spec} - the model specification
#'   \item \code{params} - validated model parameters
#' }
#'
#' @seealso
#' \code{\link{GEVStableGarch_model_spec}} for creating the specification.
#'
#' @examples
#' spec <- GEVStableGarch_model_spec(m = 1, n = 1, p = 1, q = 1)
#' model <- GEVStableGarch_set_params(
#'   spec,
#'   ar = 0.2,
#'   ma = -0.1,
#'   omega = 0.1,
#'   alpha = 0.05,
#'   beta = 0.9,
#'   dist_params = list(alpha = 1.7, beta = 0)
#' )

GEVStableGarch_set_params <- function(
    spec,
    ar = NULL,
    ma = NULL,
    omega = NULL,
    alpha = NULL,
    beta = NULL,
    aparch = NULL,
    mean = NULL,
    dist_params = NULL
) {
  if (class(spec) != "GEVStableGarch_model") stop("`spec` must be of class GEVStableGarch_spec.")
  
  ord <- spec$order
  m <- ord$m; n <- ord$n; p <- ord$p; q <- ord$q
  
  include_mean <- spec$include_mean
  aparch_flag  <- spec$aparch
  
  # ----- flags -----
  ar_flag <- m > 0
  ma_flag <- n > 0
  p_flag  <- p > 0
  q_flag  <- q > 0
  
  ## ---------- MEAN ----------
  if (include_mean) {
    if (is.null(mean) || !is.numeric(mean) || length(mean) != 1)
      stop("`mean` must be a single numeric value when include_mean = TRUE.")
  } else {
    if (!is.null(mean))
      stop("`mean` must not be supplied when include_mean = FALSE.")
  }
  
  ## ---------- APARCH ----------
  if (aparch_flag) {
    if (is.null(aparch) || !is.numeric(aparch) || length(aparch) != 1 || aparch <= 0)
      stop("`aparch` must be a single numeric value > 0 when aparch = TRUE.")
  } else {
    if (!is.null(aparch))
      stop("`aparch` must not be supplied when aparch = FALSE.")
  }
  
  ## ---------- AR ----------
  if (ar_flag) {
    if (is.null(ar) || !is.numeric(ar) || length(ar) != m)
      stop("`ar` must be numeric of length m when m > 0.")
  } else {
    if (!is.null(ar))
      stop("`ar` must not be supplied when m = 0.")
  }
  
  ## ---------- MA ----------
  if (ma_flag) {
    if (is.null(ma) || !is.numeric(ma) || length(ma) != n)
      stop("`ma` must be numeric of length n when n > 0.")
  } else {
    if (!is.null(ma))
      stop("`ma` must not be supplied when n = 0.")
  }
  
  ## ---------- omega ----------
  if (is.null(omega) || !is.numeric(omega) || length(omega) != 1 || omega <= 0)
    stop("`omega` must be a single positive numeric value.")
  
  ## ---------- alpha (p) ----------
  if (p_flag) {
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != p)
      stop("`alpha` must be numeric of length p when p > 0.")
  } else {
    if (!is.null(alpha))
      stop("`alpha` must not be supplied when p = 0.")
  }
  
  ## ---------- beta (q) ----------
  if (q_flag) {
    if (is.null(beta) || !is.numeric(beta) || length(beta) != q)
      stop("`beta` must be numeric of length q when q > 0.")
  } else {
    if (!is.null(beta))
      stop("`beta` must not be supplied when q = 0.")
  }
  
  ## ---------- distribution parameters ----------
  if (is.null(dist_params) || !is.list(dist_params))
    stop("`dist_params` must be a list.")
  if (!all(vapply(dist_params, is.numeric, logical(1))))
    stop("All `dist_params` must be numeric.")
  
  # ---------- assemble params ----------
  params <- list(
    arma = list(ar = ar, ma = ma),
    garch = list(omega = omega, alpha = alpha, beta = beta),
    aparch = aparch,
    mean = mean,
    dist_params = dist_params
  )
  
  # ---------- build model ----------
  model <- list(
    spec = spec,
    params = params
  )
  
  class(model) <- "GEVStableGarch_model"
  return(model)
}



