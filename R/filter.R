# Vectorized filters in the style of Wuertz et al. (2006): the lagged terms are
# built as shifted vectors and the recursions run through stats::filter.

# ARMA residuals e_t = x_t - mu - sum ar_i x_{t-i} - sum ma_j e_{t-j},
# conditional on e_t = 0 for the first max(m, n) observations (as fGarch)
.filter_arma <- function(x, u)
{
  N <- length(x); m <- length(u$ar); n <- length(u$ma); k <- max(m, n)
  if (k == 0) return(x - u$mu)
  idx <- (k + 1):N
  w <- x[idx] - u$mu
  for (i in seq_len(m)) w <- w - u$ar[i] * x[idx - i]
  if (n > 0) w <- as.numeric(stats::filter(w, -u$ma, method = "recursive", init = rep(0, n)))
  c(rep(0, k), w)
}

# APARCH recursion for h_t = sigma_t^delta. The first max(p, q) values are set to
# h_init, by default mean(|e|^delta)
.filter_aparch <- function(e, u, h_init = NULL)
{
  N <- length(e); p <- length(u$alpha); q <- length(u$beta); pq <- max(p, q)
  h0 <- if (is.null(h_init)) mean(abs(e)^u$delta) else h_init
  arch <- rep(0, N)
  for (i in seq_len(p)) {
    lagged <- (abs(e) - u$gamma[i] * e)^u$delta
    arch[(i + 1):N] <- arch[(i + 1):N] + u$alpha[i] * lagged[1:(N - i)]
  }
  h <- rep(h0, N)
  idx <- (pq + 1):N
  h[idx] <- if (q > 0)
    as.numeric(stats::filter(u$omega + arch[idx], u$beta, method = "recursive", init = rep(h0, q)))
  else u$omega + arch[idx]
  h
}

# Residuals, h and conditional scale of the full model
.filter_model <- function(x, u, h_init = NULL)
{
  e <- .filter_arma(x, u)
  h <- .filter_aparch(e, u, h_init)
  list(e = e, h = h, sigma = h^(1 / u$delta))
}
