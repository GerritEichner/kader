#2345678901234567890123456789012345678901234567890123456789012345678901234567890
### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/kader/R/
###  bias_var_mse.R
### Functions for computation of the estimators of bias, scaled variance and
### mse for "Kernel adjusted density estimation" of Srihera & Stute (2011)
### and for "Rank Transformations in Kernel Density Estimation" of Eichner &
### Stute (2013).
### R 3.4.1, 13./14./17.2./21./22./23.8.2017 (6./7./10.2.2015 /
###  21./24./26./28./31.10./2./4./9.11./5.12.2016)
###*****************************************************************************

#' Estimators of Bias and Scaled Variance
#'
#' ``Workhorse'' function for vectorized (in \eqn{\sigma}) computation of both
#' the bias estimator and the scaled variance estimator of eq. (2.3) in Srihera
#' & Stute (2011), and for the analogous computation of the bias and scaled
#' variance estimator for the rank-transformation method in the paragraph
#' after eq. (6) in Eichner & Stute (2013).
#'
#' @param sigma Numeric vector \eqn{(\sigma_1, \ldots, \sigma_s)} with
#'              \eqn{s \ge 1}.
#' @param Ai Numeric vector expecting \eqn{(x_0 - X_1, \ldots, x_0 - X_n) / h},
#'           where (usually) \eqn{x_0} is the point at which the density is to
#'           be estimated for the data \eqn{X_1, \ldots, X_n} with
#'           \eqn{h = n^(-1/5)}.
#' @param Bj Numeric vector expecting \eqn{(-J(1/n), \ldots, -J(n/n))} in case
#'           of the rank-transformation method, but \eqn{(\hat{\theta} - X_1,
#'           \ldots, \hat{\theta} - X_n)} in case of the non-robust
#'           Srihera-Stute-method. (Note the missing denominator \eqn{h} in
#'           comparison to argument \code{Bj} of \code{\link{adaptive_fnhat}}!)
#' @param h Numeric scalar, where (usually) \eqn{h = n^(-1/5)}.
#' @param K Kernel function with vectorized in- & output.
#' @param fnx \eqn{f_n(x_0) =} \code{mean(K(Ai))/h}, where typically
#'            \eqn{h = n^{(-1/5)}}.
#' @param ticker Logical; determines if a 'ticker' documents the iteration
#'               progress through \code{sigma}. Defaults to FALSE.
#'
#' @return A list with components \code{BiasHat} and \code{VarHat.scaled}, both
#'         numeric vectors of same length as \code{sigma}.
#'
#' @references Srihera & Stute (2011) and Eichner & Stute (2013): see
#'             \link{kader}.
bias_AND_scaledvar <- function(sigma, Ai, Bj, h, K, fnx, ticker = FALSE) {
  nr <- length(Ai)
  nc <- length(Bj)
  Ai.times.nc <- rep(Ai, times = nc)

  BV <- sapply(seq_along(sigma),
    function(j) {
      if(ticker) {  # Iteration "ticker".
        if(j == 1L) {
          message("sigma:1", appendLF = FALSE)
         } else message(",", j, appendLF = FALSE)
        }

      sig <- sigma[j]
      Y <- rep(Bj/sig, each = nr)   # Equivalent to AiBj <-
      AiBj <- Ai.times.nc + Y       # outer(Ai, Bj/sig, "+"),
      dim(AiBj) <- c(nr, nc)        # but approx. twice as fast.

      h_sig <- h / sig
      KInt <- try(stats::integrate(f = kfn_vectorized,  # Integral
            lower = -Inf, upper = Inf,                  # in bias
            K = K, xixj = AiBj, h_sig = h_sig),         # estimator.
            silent = TRUE)

      proto.biashat <- if(!inherits(KInt, "try-error"))
        KInt$value else
        { cat("\n");  print(KInt); NA }

      AiBj <- K(AiBj / h_sig)                        # Quantities whose
      varhat <- stats::var(.rowMeans(AiBj, nr, nc) + # sample var. is the
                           .colMeans(AiBj, nr, nc))  # var. estimator.
                     # Z_i on p. 431 in ES2013 or before (2.3) in SS2011.

      c(B = proto.biashat, V = varhat)
    });   if(ticker) message(".", appendLF = FALSE)

  list(BiasHat = unname(BV["B",]) / h - fnx,
       VarHat.scaled = sigma*sigma / (nr * h*h*h*h) * unname(BV["V",]))
}


#' MSE Estimator
#'
#' Vectorized (in \eqn{\sigma}) function of the MSE estimator in eq. (2.3) of
#' Srihera & Stute (2011), and of the analogous estimator in the paragraph after
#' eq. (6) in Eichner & Stute (2013).
#'
#' @inheritParams bias_AND_scaledvar
#'
#' @return A vector with corresponding MSE values for the values in
#'         \code{sigma}.
#'
mse_hat <- function(sigma, Ai, Bj, h, K, fnx, ticker = FALSE) {
  BV <- bias_AND_scaledvar(sigma = sigma, Ai = Ai, Bj = Bj, h = h, K = K,
                           fnx = fnx, ticker = ticker)
  BV$BiasHat * BV$BiasHat + BV$VarHat.scaled
}
