# #2345678901234567890123456789012345678901234567890123456789012345678901234567890
# ### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/karde/R/
# ###  kde_L2adaptive.R
# ### Functions for computation of the estimators of L2-bias, scaled L2-variance
# ### and adaptive density for "Rank Transformations in Kernel Density Estimation
# ### --- The L2-approach" by Eichner & Stute for Wroclaw 2015.
# ### R 3.3.1, 21./24./26.10.2016 (6./7./10.2.2015)
# ###*****************************************************************************
#
# #' Estimators of "L2-bias" and scaled "L2-variance"
# #'
# #' Vectorized (in sigma) computation of estimators of "L2-bias" and
# #' scaled "L2-variance".
# #'
# #' @param sigma Vector \eqn{(\sigma_1, \ldots, \sigma_s)} with \eqn{s \geq 1}.
# #' @param X Vector of order statistics, i.e., increasing sorted (!) data:
# #'          \eqn{(X_1:n, ..., X_n:n)}.
# #' @param h Scalar, where (usually) h = hadap.
# #' @param K Kernel fct. with vectorized in- & output.
# #' @param rotK Kernel fct. for "rule-of-thumb"-estimator with vectorized
# #'             in- & output.
# #' @param negJ Vector \eqn{(-J(1/n), \ldots, -J(n/n))} in case of
# #'        rank-trafo-method.
# #' @return A list with components \code{L2BiasHat} and \code{L2VarHat.scaled}.
# L2bias_hatANDvar_hat_scaled <- function( sigma, X, h, K, rotK, negJ) {
#  X.h <- X / h
#  BVXlist <- lapply( X.h,
#                     function( X.hj) { # cat( "X") # For possible process
#                                                   # feedback.
#                       Aj <- X.hj - X.h
#                       fnaXj <- mean( rotK( Aj)) / h
#                       BV <- bias_hatANDvar_hat_scaled( sigma = sigma,
#                               Aj = Aj, Bi = negJ, h = h, K = K,
#                               fnax0 = fnaXj)
#                     } )
#  BiasX <- sapply( BVXlist, "[[", "BiasHat")
#  VscaX <- sapply( BVXlist, "[[", "VarHat.scaled")
#  if( is.array( BiasX)) {
#    BiasX <- rowMeans( BiasX)
#    VscaX <- rowMeans( VscaX)
#   } else {
#    BiasX <- mean( BiasX)
#    VscaX <- mean( VscaX)
#   }
#  list( L2BiasHat = BiasX, L2VarHat.scaled = VscaX)
#  }
#
#
# #' "L2-MSE"-estimator
# #'
# #' Vectorized version (vectorized in sigma) of the computation of the
# #' "L2-MSE"-estimator.
# #'
# #' @inheritParams L2bias_hatANDvar_hat_scaled
# #' @return A vector with corresponding L2-MSE-values for the values in
# #'         \code{sigma}.
# L2mse_hat <- function( sigma, X, h, K, rotK, negJ) {
#  # print( sigma)
#  BV <- L2bias_hatANDvar_hat_scaled( sigma, X, h, K, rotK, negJ)
#  BV$L2BiasHat * BV$L2BiasHat + BV$L2VarHat.scaled
#  }
#
#
# #' Minimization of "L2-MSE" (as function of sigma).
# #'
# #' @inheritParams L2bias_hatANDvar_hat_scaled
# #' @param i Counter (currently: number of Monte-Carlo sample).
# #' @param plot Logical, indicating if graphical output should be produced.
# #' @return A list with ...
# argmin_of_L2MSE <- function( X, sigma, h, K, rotK, negJ, i = NA, plot = FALSE) {
#   # Computation of L2MSEHat on the discrete grid in sigma
#   #******************************************************
#  cat( "a") # Process feedback.
#  sigma.readjusted <- 0
#  repeat { # check if select sigma-region is too flat.
#   BV <- L2bias_hatANDvar_hat_scaled( sigma = sigma, X = X, h = h, rotK, K = K,
#                                      negJ = negJ)
#   L2VarHat.scaled <- BV$L2VarHat.scaled
#   L2BiasHat.squared <- BV$L2BiasHat * BV$L2BiasHat
#   L2MSEHat <- L2BiasHat.squared + L2VarHat.scaled
#   if( plot) {
#    mainstr <- NA # paste0( "MC sample ", i, ": Rang-trafo-L2-method")
#    matplot( sigma, cbind( L2BiasHat.squared, L2VarHat.scaled, L2MSEHat),
#             type = "l", lty = c( "twodash", "longdash", "solid"),
#             lwd = rep( 3, 3), col = "black", #c( "blue", "violet", "magenta3"),
#             main = mainstr, ylim = c( 0, max( L2MSEHat, na.rm = TRUE)),
#             xlab = expression( sigma),
#             ylab = NA) # expression( widehat(MSE)(sigma))
#
#    # The legend for (Bias_L2(\sigma))^2 and Var_L2(\sigma):
#     legend( "top", lty = c( "twodash", "longdash", "solid"), lwd = 3,
# #            col = c( "blue", "violet", "magenta3"),
#             cex = 1.2, bty = "n", y.intersp = 1.5, # horiz = TRUE,
#             legend = c( expression( bar(widehat(plain(Bias)))^2~(sigma)),
#                         expression( sigma^2*n^{-1}*h^{-4}~
#                                      bar(widehat(plain(Var)))(sigma)),
#                         expression( bar(widehat(plain(MISE)))(sigma))) )
#
#    if( sigma.readjusted > 0)
#     legend( "topleft", legend = "'Sigma' extended!", bty = "n", cex = 1.3)
#    }
#
#   if( max( diff( range( L2BiasHat.squared)), diff( range( L2VarHat.scaled)))
#       > .Machine$double.eps^0.25) break
#
#   if( plot) dev.set( dev.cur())
#   warning( "L2Biashat.squared or L2VarHat.scaled were - almost - constant.")
#   sigma <- 10 * sigma;   sigma.readjusted <- sigma.readjusted + 1
#   cat( "\nNew sigma:", sigma[ 1], "-", sigma[ length( sigma)], "\n")
#   }    # end of repeat{...}
#
#
#   # Minimization of the estimated MSE in sigma (more precise: on the
#   # sigma-grid) in four steps:
#   # Step 1: Location of (first) maximum of L2VarHat.scaled (!) on the grid
#   #***********************************************************************
#  cat( "b") # Process feedback.
#  varhatscaled.maxidx <- which.max( L2VarHat.scaled) # First sigma-index, where
#                                                     # L2VarHat.scaled maximal
# # if( plot) points( sigma[ varhatscaled.maxidx],
# #                   L2VarHat.scaled[ varhatscaled.maxidx],
# #                   pch = 6, cex = 2, col = "violet", lwd = 2)
#
#   # Step 2: Location of (first) minimum of MSEHat on the grid LEFT
#   #         of the location of (first) maximum of VarHat.scaled
#   #***************************************************************
#  msehat.minidx <- which.min( L2MSEHat[ 1:varhatscaled.maxidx])
#                    # First sigma-index, where L2MSEHat left of first maximum
#                    # of L2VarHat.scaled is minimal.
#  sigadap <- sigma[ msehat.minidx] # sigma-value, for which L2MSEHat left of
#                                   # first maximum of L2VarHat.scaled is minimal.
#  L2msehat.min <- L2MSEHat[ msehat.minidx] # Pertaining minimal L2MSEHat-value.
#
#   # Neighborhood of location of minimum, within which the search continues:
#  sigrange <- c(      sigma[ max(       1,        msehat.minidx - 3)],
#                 min( sigma[ min( length( sigma), msehat.minidx + 3)],
#                      sigma[ varhatscaled.maxidx]) )
#  if( diff( sigrange) < .Machine$double.eps^0.25) sigrange <- range( sigma)
# # if( plot) segments( x0 = sigrange[ 1], y0 = par( "usr")[ 3],
# #                     x1 = sigrange[ 2], col = "cyan1", lwd = 2)
#
#   # Step 3: Numerical minimization in neighborhood of discrete
#   # minimum of L2MSEHat
#   #***********************************************************
#  cat( "c") # Process feedback.
#  L2msehat.opt <- optimize( L2mse.hat, interval = sigrange, X = X,
#                            h = h, K = K, negJ = negJ)
#
#  if( plot) {
#   points( c( L2msehat.opt$minimum,   sigadap)[1],
#           c( L2msehat.opt$objective, L2msehat.min)[1],
#           pch = c( 3, 4)[2], cex = 2, # col = c( "green", "red"),
#           lwd = 2)
# #  legend( "top", pch = c( 6, NA, 3, 4), lty = c( NA, "solid", NA, NA),
# #          col = c( "violet", "cyan1", "green", "red"), lwd = 2, bty = "n",
# #          legend = c( expression( paste( "Max. of ", sigma^2/(n*h^4) *
# #                                           widehat(Var)(sigma))),
# #                      expression( paste( "Search-range for num. min. of ",
# #                                         widehat(MSE)(sigma))),
# #                      substitute( "Num. min." == group( "(", list( x, y), ")"),
# #                            list( x = round( L2msehat.opt$minimum, 4),
# #                                  y = round( L2msehat.opt$objective, 4))),
# #                      substitute( "Discr. min." ==
# #                                    group( "(", list( x, y), ")"),
# #                                  list( x = round( sigadap, 4),
# #                                        y = round( L2msehat.min, 4))) ) )
#  }
#
#   # Step 4: Control, if numerically found minimum is in fact
#   #         "better" than the discrete one.
#   #*********************************************************
#  discrete.minimum.less.than.numerical <- FALSE
#  if( L2msehat.opt$objective < L2msehat.min) {
#    sigadap <- L2msehat.opt$minimum
#    L2msehat.min <- L2msehat.opt$objective
#   } else {
# #   cat( "\ngrid search '<=' optim()\n")
#    discrete.minimum.less.than.numerical <- TRUE
# #   if( plot) legend( "center", legend = "Minimum adjusted!",
# #                     bty = "n", cex = 1.3)
#   }
#
#  list( sigma = sigadap, L2msehat = L2msehat.min,
#        sigma.adjusted = sigma.readjusted,
#        discmin.less = discrete.minimum.less.than.numerical)
# }
#
#
#
# #' Kernel adaptive density estimator for the L2-adapted value of sigma.
# #'
# #' @param sigmaL2 Scalar with L2-adapted value of sigma.
# #' @param Aj Vector (x0 - X_1, ..., x0 - X_n) / h, where (usually) h = hadap.
# #' @param Bi Vector (-J(1/n), ..., -J(n/n)) in case of rank-trafo-method,
# #'        and (hat theta - X_1, ..., hat theta - X_n) for S.-S.-method.
# #' @param h Scalar, where (usually) h = hadap.
# #' @param K Kernel function with vectorized in- & output, (usually) K = adapK.
# #' @return A list with components \code{fnx0}, \code{sigma} and \code{h}.
# fn_L2adaptive <- function( sigmaL2, Aj, Bi, h, K) {
#  K.arg <- outer( sigmaL2 / h * Aj, Bi / h, "+")
#  K.arg <- K.arg[ c( FALSE, rep( TRUE, nrow( K.arg)))]  # Without diagonal.
#
#  list( fnx0 = sigmaL2 / (h * h) * mean( K( K.arg)), sigma = sigmaL2, h = h)
#  }
