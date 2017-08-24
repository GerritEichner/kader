### Estimation and minimization of the estimated AMSE function in
### Eichner & Stute (2012): "Kernel adjusted nonparametric regression", Journal
### of Statistical Planning and Inference.
### R 3.1.1, 11.8.2014
###*****************************************************************************

 ## Helper functions
 ##*****************

 # Kernel function
 #****************
Kernel <- function( x) dnorm( x)   # Gauss kernel 

 # Regression function
 #********************
 # Here: a polynomial of degree 4 with one maximum (or minimum), one point of
 # inflection, and one saddle point.
 # Memo: for p(x) = a * (x - x1) * (x - x2)^3 + b the max. (or min.) is at
 # x = (3 * x1 + x2) / 4, the pt. of inflection at x = (x1 + x2) / 2, and the
 # saddle pt. at x = x2.
 # Note: for poly4()'s default values a minimum is at x = 2, a pt. of
 # inflection at x = 4, and a saddle pt. an x = 8; an "arbitrary" pt. would,
 # e.g., be at x = 5.
poly4 <- function( x, x1 = 0, x2 = 8, a = 0.01, b = 0)
 a * (x - x1) * (x - x2)^3 + b

 # Nadaraya-Watson estimator (m_n in eq. (1.1) in the paper)
 #**********************************************************
NadWat <- function( x, X, Y, h, K = Kernel) {
  # x:      vector (x_1, ..., x_r) or scalar
  # X, Y:   vectors (X_1, ..., X_n), (Y_1, ..., Y_n)
  # h:      scalar
  # K:      fct. with vectorial in- & output
 M <- K( outer( x/h, X/h, "-"))   # r x n
 drop( (M %*% Y) / rowSums( M))   # r x 1
 }

 # Weights (W_{ni} below eq. (2.1) in the paper)
 #**********************************************
W_n <- function( sigma, h, xXh, thetaXh, K = Kernel) {
  # sigma:     vector (sigma_1, ..., sigma_s) or scalar
  # h:         scalar
  # xXh:       vector ((x - X_1)/h, ..., (x - X_n)/h)
  # thetaXh:   vector ((\theta - X_1)/h, ..., (\theta - X_n)/h)
  # K:         Fct. with vectorial in- & output
 A <- outer( sigma/h, xXh)       # s x n
 A <- outer( A, thetaXh, "+")    # s x n x n
 A <- rowSums( K( A), dims = 2)  # s x n
 A / rowSums( A)   # s x n: (W_{ni}( x, \sigma_r))_{1<=r<=s, 1<=i<=n}
 }

 # Bias estimator (Bias_n(\sigma) in 2nd paragraph below eq. (2.2) in paper)
 #**************************************************************************
Bias_n <- function( sigma, h, xXh, thetaXh, K = Kernel, mmDiff) {
  # sigma, X, h, xXh, thetaXh, K:   see W_n()!
  # mmDiff:   vector (m_n(X_1) - m_n(x), ..., m_n(X_n) - m_n(x))
 drop( W_n( sigma, h, xXh, thetaXh, K) %*% mmDiff)  # s x 1
 }

 # Variance estimator (Var_n(\sigma) in 2nd paragraph below eq. (2.2) in paper)
 #*****************************************************************************
Var_n <- function( sigma, h, xXh, thetaXh, K = Kernel, YmDiff2) {
  # sigma, h, xXh, thetaXh, K:   see W_n()!
  # YmDiff2:   vector ( (Y_1 - m_n(x))^2, ..., (Y_n - m_n(x))^2 )
 wx <- W_n( sigma, h, xXh, thetaXh, K)   # s x n
 drop( (wx * wx) %*% YmDiff2)  # s x 1
 }


 ## An exemplary simulation
 ##************************

 # Simulation parameters and data generation
 #******************************************
m  <- poly4   # Choice of the regression function (here: a degree-4-polynomial)
x0 <- 1   # The point x_0 at which the AMSE-optimal kernel adjusted nonpara-
          # metric estimation of m should take place. (Recall: for poly4()'s
          # default values a minimum is at 2, a pt. of inflection at 4, and a
          # saddle pt. an 8; an "arbitrary" pt. would, e.g., be at 5.)

n    <- 50   # Sample size
sdeps <- 1   # Standard deviation of the \epsilon_i: \sqrt(Var(Y|X=x))
             # (here: constant in x)

design.ctr <- x0 + 0.5   # "centre" and "scale" of the design, i.e. in the
design.scl <- 1          # normal scenario below, expected value and std. dev., 
                         # resp., of the distribution of the x_i's (i >= 1)

set.seed( 42)   # starting value of the pseudo-random numbers generator
                # to allow reproducibility of (x_i, Y_i)'s, i = 1, ..., n
x <- rnorm( n, mean = design.ctr, sd = design.scl)   # the design x_1, ... x_n
Y <- m( x) + rnorm( length( x), sd = sdeps)          # the Y_i's


xgrid <- design.ctr + seq( -3.3 * design.scl, 3.3 * design.scl, length = 101)
   # arbitrarily chosen x-grid, on which the regression function m will be
   # drawn.



 ## Computation of the AMSE-optimal kernel adjusted
 ## nonparametric regression estimator of m(x0)
 ##************************************************

 # Range within which and grid on which the scale parameter \sigma of the
 # adaptive method will be considered:
sig.range <- c( 0.01, 1.5 * diff( range( xgrid)))
sig.grid <- seq( sig.range[ 1], sig.range[ 2], length = 501)

 # Choice of the window width h for the adaptive method:
h <- n^(-1/5)

 # Using the Nadaraya-Watson estimators m_n(x_0) and m_n(x_1), ..., m_n(x_n)
 # with the window width h from above to derive quantities for the estimation
 # of Var(\sigma) and Bias(\sigma) at x_0:
mnx0     <- NadWat( x0, x, Y, h)   # m_n(x_0)
mnx      <- NadWat( x, x, Y, h)    # m_n(x_i) for i = 1, ..., n
mnx_mnx0 <- mnx - mnx0             # m_n(x_i) - m_n(x_0) for i = 1, ..., n
Y_mnx2   <- (Y - mnx)^2            # (Y_i - m_n(x_i))^2 for i = 1, ..., n

 # Quantities that do not depend on \sigma, but will repeatedly be needed in
 # what follows (makes the algorithmn slightly more efficient):
x0xh <- (x0 - x) / h;     thetaxh <- (mean( x) - x) / h

 # Estimators of Var_x0(\sigma) and Bias_x0(\sigma) on the \sigma-grid (to
 # also be able to draw them for visualisation purposes later):
Bn  <- Bias_n( sigma = sig.grid, h = h, xXh = x0xh, thetaXh = thetaxh,
               mmDiff = mnx_mnx0)
Vn2 <- Var_n( sigma = sig.grid, h = h, xXh = x0xh, thetaXh = thetaxh,
              YmDiff2 = Y_mnx2)

 # Composing the estimator of AMSE_x0(\sigma) on the \sigma-grid using the
 # estimators of Var_x0(\sigma) and Bias_x0(\sigma) on the \sigma-grid (as
 # in 2nd paragraph below eq. (2.2) in paper):
An <- Vn2 + Bn * Bn

 # Minimizer and minimum of the estimator of AMSE_x0(\sigma) using "grid
 # search" on the \sigma-grid:
min_idx <- which.min( An)
sig_opt_x0 <- sig.grid[ min_idx];   amse_min <- An[ min_idx]

 # The AMSE-optimal kernel adjusted nonparametric regression estimator
 # at x_0, i.e., \hat m_n(x_0):
W <- W_n( sigma = sig_opt_x0, h = h, xXh = x0xh, thetaXh = thetaxh)

mnhatx0 <- W %*% Y   # = \hat m_n(x_0)





 ## Grafical display for the current data set
 ##******************************************
 # preparing the layout of the grafics device
op <- par( mfrow = c( 3,1), mar = c( 3,3,2,0.1)+0.1, mgp = c( 1.5,0.5,0),
           tcl = -0.3, cex.main = 2)

 # The scatter plot of the "raw data":
plot( Y ~ x, xlim = range( x, xgrid), ylim = range( Y, mnhatx0, na.rm = TRUE),
      main = bquote( n == .(n)), xlab = "x", ylab = "y")

 # The "true" regression function m:
lines( xgrid, m( xgrid), lty = 2) 

 # The AMSE-optimal kernel adjusted nonparametric regression estimator
 # at x_0, i.e., the point (x_0, \hat m_n(x_0)):
points( x0, mnhatx0, col = "red", pch = 4, cex = 2)

 # The legend for the "true" regression function m and for the point
 # (x_0, \hat m_n(x_0)):
legend( "top", lty = c( 2, NA), pch = c( NA, 4), col = c( "black", "red"),
        bty = "n", cex = 1.5,
        legend = c( as.expression( bquote( paste( "m  with  ",
                                                  sigma( paste( Y, "|", X == x))
                                                   == .(sdeps)) ) ),
                    as.expression( bquote( paste( hat(m)[n](x[0]), "  at  ",
                                                  x[0] == .(x0)) ) ) ) )

 # Visualisation of the estimators of (Bias_x0(\sigma))^2 and of Var_x0(\sigma)
 # on the \sigma-grid:
matplot( sig.grid, cbind( Bn*Bn, Vn2), type = "l", lty = 1:2,
         col = c( "black", "red"), xlab = expression( sigma), ylab = "")

 # The legend for (Bias_x0(\sigma))^2 and Var_x0(\sigma):
legend( "top", lty = 1:2, col = c( "black", "red"), bty = "n", cex = 1.5,
        legend = c( expression( paste( widehat(plain(Bias))[x[0]]^2, (sigma))),
                    expression( widehat(plain(Var))[x[0]](sigma))))

 # Visualisation of the estimator of AMSE_x0(\sigma) on the \sigma-grid
 # together with the point indicating the detected minimum, and a legend:
plot( sig.grid, An, type = "l", xlab = expression( sigma), ylab = "")
points( sig_opt_x0, amse_min, pch = 4, col = "red", cex = 2)
legend( "top", lty = c( 1, NA), pch = c( NA, 4), col = c( "black", "red"),
        bty = "n", cex = 1.5,
        legend = c( expression( widehat( plain( AMSE))[x[0]](sigma)),
                    substitute( group( "(", list( plain( Minimizer),
                                                  plain( Minimum)), ")")
                                   == group( "(", list( x, y), ")"),
                                list( x = signif( sig_opt_x0, 4),
                                      y = signif( amse_min, 4))) ) )

par( op)   # Restoring the previous settings of the graphics device.
