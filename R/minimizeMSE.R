#2345678901234567890123456789012345678901234567890123456789012345678901234567890
### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/kader/R/
###  minimizeMSE.R
### Minimization of MSE-estimator for "Kernel adjusted density estimation" of
### Srihera & Stute (2011) and for "Rank Transformations in Kernel Density
### Estimation" of Eichner & Stute (2013).
### R 3.4.2, 13./14./15./16./17./24.2./22./23.8./28.9./2.10.2017
###  (21./24./26./27./28./31.10./4./9.11./5.12.2016)
###*****************************************************************************


#' Minimization of Estimated MSE
#'
#' Minimization of the estimated MSE as function of \eqn{\sigma} in four steps.
#'
#' Step 1: determine first (= smallest) maximizer of \code{VarHat.scaled} (!)
#' on the grid in \code{sigma}. Step 2: determine first (= smallest) minimizer
#' of estimated MSE on the \eqn{\sigma}-grid LEFT OF the first maximizer of
#' \code{VarHat.scaled}. Step 3: determine a range around the yet-found
#' (discrete) minimizer of estimated MSE within which a finer search for the
#' ``true'' minimum is continued using numerical minimization. Step 4: check if
#' the numerically determined minimum is indeed better, i.e., smaller than the
#' discrete one; if not keep the first.
#'
#' @inheritParams bias_AND_scaledvar
#' @param VarHat.scaled Vector of estimates of the scaled variance
#'                      (for values of \eqn{\sigma} in \code{sigma}).
#' @param BiasHat.squared Vector of estimates of the squared bias
#'                        (for values of \eqn{\sigma} in \code{sigma}).
#' @param plot Should graphical output be produced? Defaults to \code{FALSE}.
#' @param \ldots Currently ignored.
#'
#' @return A list with components \code{sigma.adap}, \code{msehat.min} and
#'         \code{discr.min.smaller} whose meanings are as follows:
#'   \tabular{ll}{
#'    \code{sigma.adap} \tab Found minimizer of MSE estimator. \cr
#'    \code{msehat.min} \tab Found minimum of MSE estimator. \cr
#'    \code{discr.min.smaller} \tab TRUE iff the numerically found minimum was
#'                             smaller than the discrete one.
#'    }
#'
#' @examples
#' require(stats)
#'
#' set.seed(2017);     n <- 100;     Xdata <- sort(rnorm(n))
#' x0 <- 1;      Sigma <- seq(0.01, 10, length = 11)
#'
#' h <- n^(-1/5)
#' Ai <- (x0 - Xdata)/h
#' fnx0 <- mean(dnorm(Ai)) / h   # Parzen-Rosenblatt estimator at x0.
#'
#'  # For non-robust method:
#' Bj <- mean(Xdata) - Xdata
#'#  # For rank transformation-based method (requires sorted data):
#'# Bj <- -J_admissible(1:n / n)   # rank trafo
#'
#' BV <- kader:::bias_AND_scaledvar(sigma = Sigma, Ai = Ai, Bj = Bj,
#'   h = h, K = dnorm, fnx = fnx0, ticker = TRUE)
#'
#' kader:::minimize_MSEHat(VarHat.scaled = BV$VarHat.scaled,
#'   BiasHat.squared = (BV$BiasHat)^2, sigma = Sigma, Ai = Ai, Bj = Bj,
#'   h = h, K = dnorm, fnx = fnx0, ticker = TRUE, plot = FALSE)
#'
minimize_MSEHat <- function(VarHat.scaled, BiasHat.squared, sigma, Ai, Bj, h,
  K, fnx, ticker = FALSE, plot = FALSE, ...) {

 message("Minimizing MSEHat:")

 # Step 1: determine first (= smallest) maximizer of VarHat.scaled (!!!)
 #         on the grid in sigma.
 #**********************************************************************
 message("Step 1: Search smallest maximizer of VarHat.scaled on sigma-grid.")
 varhatscaled.maxidx <- which.max(VarHat.scaled) # 1st sigma-index, where
                                                 # VarHat.scaled is maximal.

 if(plot) {
   graphics::points(sigma[varhatscaled.maxidx], # mark maximum in graph.
                    VarHat.scaled[varhatscaled.maxidx],
                    pch = 6, cex = 2, col = "violet", lwd = 2)
 }

 # Step 2: determine first (= smallest) minimizer of MSEHat on the sigma-grid
 #         LEFT OF the first maximizer of VarHat.scaled.
 #***************************************************************************
 message("Step 2: Search smallest minimizer of MSEHat on sigma-grid to the\n",
         "        LEFT of just found smallest maximizer of VarHat.scaled.")
 MSEHat <- BiasHat.squared + VarHat.scaled
 msehat.minidx <- which.min(MSEHat[ 1:varhatscaled.maxidx])
  # First sigma-index where MSEHat is minimal left of first
  # maximizer of VarHat.scaled.

 sigadap <- sigma[ msehat.minidx]
  # sigma-value for which MSEHat is minimal
  # left of first maximizer of VarHat.scaled.
 msehat.min <- MSEHat[ msehat.minidx]  # Pertaining minimal MSEHat-value


 # Step 3: determine a range around the yet-found (discrete) minimizer of
 #         MSEHat within which a finer search for the "true" minimum con-
 #         tinues using numerical minimization
 #***********************************************************************
 message("Step 3: Finer search for 'true' minimum of MSEHat using\n",
         "        numerical minimization. (May take a while.)")
 sigrange <- c(    sigma[max(      1,       msehat.minidx - 3)],
               min(sigma[min(length(sigma), msehat.minidx + 3)],
                   sigma[varhatscaled.maxidx]) )
 if(diff(sigrange) < .Machine$double.eps^0.25) sigrange <- range(sigma)

 if(plot) {
   graphics::segments(
     x0 = sigrange[1], y0 = graphics::par("usr")[3] / 2,   # mark range
     x1 = sigrange[2], col = "cyan1", lwd = 2)             # in graph
 }

 # Numerical minimization:
 # suppressMessages( {
 msehat.opt <- stats::optimize(mse_hat, interval = sigrange, Ai = Ai, Bj = Bj,
   h = h, K = K, fnx = fnx, ticker = ticker);  if(ticker) message("")
 # } )


 # Step 4: check if the numerically determined minimum is indeed better,
 #         i.e., smaller than the discrete one! And if not keep the first.
 #************************************************************************
 message("Step 4: Check if numerically determined minimum is smaller\n",
         "        than discrete one.")
 discrete.minimum.smaller.than.numerical <- FALSE
 if(msehat.opt$objective < msehat.min) {
   sigadap <- msehat.opt$minimum
   msehat.min <- msehat.opt$objective
   message("        Yes, optimize() was 'better' than grid search.")
 } else {
   message("        No, grid search was 'at least as good' as optimize().")
   discrete.minimum.smaller.than.numerical <- TRUE
   if(plot) {
     graphics::legend("bottom", legend = "Minimum adjusted!", bty = "n",
       cex = 1.3)
   }
 }


 if(plot) {
   graphics::points(c(msehat.opt$minimum, sigadap),       # mark "both min-
                    c(msehat.opt$objective, msehat.min),  # nima" in graph
                    pch = c(3, 4), cex = 2, col = c("green", "red"), lwd = 2)

   graphics::legend("topright", pch = c(6, NA, 3, 4), # cex = 0.8,
     lty = c(NA, "solid", NA, NA), lwd = 2, bty = "n", inset = c(0.005, 0),
     col = c("violet", "cyan1", "green", "red"),
     legend = c(
       expression("Found max. of"~sigma^2/(n*h^4) * widehat(Var)(sigma)),
       expression("Search-range for num. min. of"~widehat(MSE)(sigma)),
       substitute("Num. min." == group("(", list(x, y), ")"),
         list(x = round(msehat.opt$minimum, 4),
              y = round(msehat.opt$objective, 4))),
       substitute("Discr. min." == group("(", list(x, y), ")"),
         list(x = round(sigadap, 4),
              y = round(msehat.min, 4))) ) )
   }

 list(sigma.adap = sigadap, msehat.min = msehat.min,
      discr.min.smaller = discrete.minimum.smaller.than.numerical)
}
