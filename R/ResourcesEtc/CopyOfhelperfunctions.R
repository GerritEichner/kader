#2345678901234567890123456789012345678901234567890123456789012345678901234567890
### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/kader/R/
###  helperfunctions.R
### Helper functions for "Kernel adjusted density estimation" of Srihera &
### Stute (2011) and for "Rank Transformations in Kernel Density Estimation"
### of Eichner & Stute (2013)
### R 3.3.2, 21./24./26./27./28./31.10./4./9.11.2016
###*****************************************************************************

#' Epanechnikov kernel
#'
#' Vectorized evaluation of the Epanechnikov kernel.
#'
#' @param x Vector.
#' @return A vector of the Epanechnikov kernel evaluated at the values in
#'         \code{x}.
epanechnikov <- function(x)
  (abs(x) < sqrt(5)) * 3 / (4 * sqrt(5)) * (1 - x * x / 5)


#' Rectangular kernel
#'
#' Vectorized evaluation of the rectangular kernel.
#'
#' @param x Vector.
#' @param a Scalar. Lower bound of kernel support; defaults to -0.5.
#' @param b Scalar. Upper bound of kernel support; defaults to 0.5.
#' @return A vector of the rectangular kernel evaluated at the values in
#'         \code{x}.
rectangular <- function(x, a = -0.5, b = 0.5) (a < x & x < b) / (b - a)


#' Convolution of a kernel function K with fn
#'
#' Vectorized evaluation of the convolution of the kernel function K with fn.
#'
#' Vectorized (in u) evaluation of - a more explicit representation of - the
#' integrand \eqn{K(u) * f_n( \ldots - h^2/\sigma * u)} that is used in the
#' computation of the bias estimator before eqn. (2.3) in S. & S. (2011).
#' Also used for the analogous computation of the respective bias estimator
#' in the paragraph after eqn. (6) in E. & S. (2013).
#'
#' @param u Vector.
#' @param K Kernel function with vectorized in- & output.
#' @param xixj Matrix.
#' @param h_sig Scalar.
#'
#' @return A vector of \eqn{(K * f_n)(u)} evaluated at the values in \code{u}.
#'
#' @section Note:
#'  An alternative implementation could be
#'  \code{K(u) * sapply(h_sig * u, function(v) mean(K(xixj - v)))}
kfn_vectorized <- function( u, K, xixj, h_sig) {
  XU <- outer( xixj, h_sig * u, "-")
  K( u) * colMeans( K( XU), dims = 2)
}



#' Minimization of the estimated MSE
#'
#' Minimization of the estimated MSE as function of \eqn{\sigma} in four steps.
#'
#' @inheritParams bias_AND_scaledvar
#' @param VarHat.scaled Vector of estimates of the scaled variance
#'                      (for values of \eqn{\sigma} in \code{sigma}).
#' @param BiasHat.squared Vector of estimates of the squared bias
#'                        (for values of \eqn{\sigma} in \code{sigma}).
#'
#' @return A list with components \code{sigma.adap}, \code{msehat.min} and
#'         \code{discr.min.smaller} whose meanings are as follows:
#'   \tabular{ll}{
#'    \code{sigma.adap} \tab Minimizer of MSE-estimator. \cr
#'    \code{msehat.min} \tab Minimum of MSE-estimator. \cr
#'    \code{discr.min.smaller} \tab TRUE iff the numerically found minimum was
#'                             smaller than the discrete one.
#'    }
minimize_MSEHat <- function( VarHat.scaled, BiasHat.squared,
                             sigma, Ai, Bj, h, K, fnx, plot = FALSE, ...) {

# Step 1: determine first (= smallest) maximizer of VarHat.scaled (!!!)
#         on the grid in sigma.
message( "Step 1: Smallest maximizer of VarHat.scaled on sigma-grid.")
varhatscaled.maxidx <- which.max( VarHat.scaled) # 1st sigma-index, where
                                                 # VarHat.scaled is maximal.

if( plot) points( sigma[ varhatscaled.maxidx], # mark maximum in graph.
  VarHat.scaled[ varhatscaled.maxidx],
  pch = 6, cex = 2, col = "violet", lwd = 2)


# Step 2: determine first (= smallest) minimizer of MSEHat on the sigma-grid
#         LEFT OF the first maximizer of VarHat.scaled.
#***************************************************************************
message( "Step 2: Smallest minimizer of MSEHat on sigma-grid LEFT OF ",
         "smallest maximizer of VarHat.scaled.")
MSEHat <- BiasHat.squared + VarHat.scaled
msehat.minidx <- which.min( MSEHat[ 1:varhatscaled.maxidx])
 # 1. sigma-index where MSEHat is minimal left of 1st
 # maximizer of VarHat.scaled.

sigadap <- sigma[ msehat.minidx]
 # sigma-value for which MSEHat is minimal
 # left of 1st maximizer of VarHat.scaled.
msehat.min <- MSEHat[ msehat.minidx]  # Pertaining minimal MSEHat-value


# Step 3: determine a range around the yet-found (discrete) minimizer of
#         MSEHat within which a finer search for the "true" minimum con-
#         tinues using numerical minimization
#***********************************************************************
message( "Step 3: Finer search for 'true' minimum of MSEHat ",
         "using numerical minimization.")
sigrange <- c(      sigma[ max(       1,        msehat.minidx - 3)],
  min( sigma[ min( length( sigma), msehat.minidx + 3)],
    sigma[ varhatscaled.maxidx]) )
if( diff( sigrange) < .Machine$double.eps^0.25) sigrange <- range( sigma)

if( plot) segments( x0 = sigrange[ 1], y0 = par( "usr")[ 3],     # mark range
  x1 = sigrange[ 2], col = "cyan1", lwd = 2)   # in graph

# Numerical minimization:
suppressMessages(
msehat.opt <- optimize( mse_hat, interval = sigrange, Ai = Ai, Bj = Bj, h = h,
  K = K, fnx = fnx)
)


# Step 4: check if the numerically determined minimum is indeed better,
#         i.e., smaller than the discrete one! And if not keep the first.
#************************************************************************
message( "Step 4: Check if numerically determined minimum is smaller ",
         "than discrete one.")
discrete.minimum.smaller.than.numerical <- FALSE
if( msehat.opt$objective < msehat.min) {
  sigadap <- msehat.opt$minimum
  msehat.min <- msehat.opt$objective
} else {
  cat( "\nGrid search 'better than' optimize().\n")
  discrete.minimum.smaller.than.numerical <- TRUE
  if( plot) legend( "top", legend = "Minimum adjusted!", bty = "n", cex = 1.3)
}


if( plot) {
  points( c( msehat.opt$minimum,   sigadap), # mark "both mini-
    c( msehat.opt$objective, msehat.min),    # ma" in graph
    pch = c( 3, 4), cex = 2, col = c( "green", "red"), lwd = 2)

  legend( "topright", pch = c( 6, NA, 3, 4), lty = c( NA, "solid", NA, NA),
    col = c( "violet", "cyan1", "green", "red"), lwd = 2, bty = "n",
    legend = c( expression( paste( "Max. of ", sigma^2/(n*h^4) *
        widehat(Var)(sigma))),
      expression( paste( "Search-range for num. min. of ",
        widehat(MSE)(sigma))),
      substitute( "Num. min." == group( "(", list( x, y), ")"),
        list( x = round( msehat.opt$minimum, 4),
          y = round( msehat.opt$objective, 4))),
      substitute( "Discr. min." ==
          group( "(", list( x, y), ")"),
        list( x = round( sigadap, 4),
          y = round( msehat.min, 4))) ) )
}

list( sigma.adap = sigadap, msehat.min = msehat.min,
      discr.min.smaller = discrete.minimum.smaller.than.numerical)
}
