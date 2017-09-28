#2345678901234567890123456789012345678901234567890123456789012345678901234567890
### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/kader/R/
###  kre_adaptive.R
### Functions for computing the adaptive regression estimator for "Kernel
### adjusted nonparametric regression" of Eichner & Stute (2012).
### R 3.4.1, 16./17./24.2./21./22./23./24.8./28.9.2017 (11.8.2014 / 25.11.2016)
###*****************************************************************************

#' The Classical Nadaraya-Watson Regression Estimator
#'
#' In its arguments \code{x} and \code{dataX} vectorized function to compute
#' the classical Nadaraya-Watson estimator (as it is \eqn{m_n} in eq. (1.1)
#' in Eichner & Stute (2012)).
#'
#' Implementation of the classical Nadaraya-Watson estimator as in eq. (1.1) in
#' Eichner & Stute (2012) at given location(s) in \code{x} for data \eqn{(X_1,
#' Y_1), \ldots, (X_n, Y_n)}, a kernel function \eqn{K} and a bandwidth \eqn{h}.
#'
#' @export
#'
#' @param x Numeric vector with the location(s) at which the Nadaraya-Watson
#'          regression estimator is to be computed.
#' @param dataX Numeric vector \eqn{(X_1, \ldots, X_n)} of the x-values from
#'              which (together with the pertaining y-values) the estimate is
#'              to be computed.
#' @param dataY Numeric vector \eqn{(Y_1, \ldots, Y_n)} of the y-values from
#'              which (together with the pertaining x-values) the estimate is
#'              to be computed.
#' @param K A kernel function (with vectorized in- & output) to be used for the
#'          estimator.
#' @param h Numeric scalar for bandwidth \eqn{h}.
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @examples
#' require(stats)
#'
#'  # Simulation parameters and data generation
#'  #******************************************
#'  # Regression function: a polynomial of degree 4 with one maximum (or
#'  # minimum), one point of inflection, and one saddle point.
#'  # Memo: for p(x) = a * (x - x1) * (x - x2)^3 + b the max. (or min.)
#'  # is at x = (3*x1 + x2)/4, the point of inflection is at x =
#'  # (x1 + x2)/2, and the saddle point at x = x2.
#' m <- function(x, x1 = 0, x2 = 8, a = 0.01, b = 0) {
#'  a * (x - x1) * (x - x2)^3 + b
#' }
#'  # Note: for m()'s default values a minimum is at x = 2, a point of
#'  # inflection at x = 4, and a saddle point an x = 8; an "arbitrary"
#'  # point would, e.g., be at x = 5.
#'
#' n <- 100      # Sample size.
#' set.seed(42)   # to guanrantee reproducibility.
#' X <- runif(n, min = -3, max = 15)        # x_1, ..., x_n
#' Y <- m(X) + rnorm(length(X), sd = 5)   # Y_1, ..., Y_n
#'
#' x <- seq(-3, 15, by = 0.5) # Locations at which the Nadaraya-Watson kernel
#'                           # estimator of m shall be computed.
#'
#' fnhat <- nadwat(x = x, dataX = X, dataY = Y, K = dnorm, h = n^(-1/5))
#'
#' plot(x = X, y = Y)
#' lines(x = x, y = fnhat, col = "blue")
#' curve(m, add = TRUE, col = "red")
#'
nadwat <- function(x, dataX, dataY, K, h) {
 M <- K(outer(x/h, dataX/h, "-"))   # length(x) x length(Y)
 drop((M %*% dataY) / rowSums(M))   # length(x) x 1
 }



#' Weights \eqn{W_{ni}} of Eichner & Stute (2012)
#'
#' Function, vectorized in its first argument \code{sigma}, to compute the
#' ``updated'' weights \eqn{W_{ni}} in eq. (2.1) of Eichner & Stute (2012) for
#' the kernel adjusted regression estimator.
#'
#' Note that it is not immediately obvious that \eqn{W_{ni}} in eq. (2.1) of
#' Eichner & Stute (2012) is a function of \eqn{\sigma}. In fact, \eqn{W_{ni}
#' = W_{ni}(x; h, \theta, \sigma)} as can be seen on p. 2542 ibid. The
#' computational version implemented here, however, is given in (15.19) of
#' Eichner (2017). Pre-computed \eqn{(x - X_i)/h} and \eqn{(\theta - X_i)/h},
#' \eqn{i = 1, \ldots, n} are expected for efficiency reasons (and are
#' currently prepared in function \code{kare}).
#'
#' @export
#'
#' @param sigma Numeric vector \eqn{(\sigma_1, \ldots, \sigma_s)} with
#'              \eqn{s \ge 1} with values of the scale parameter \eqn{\sigma}.
#' @param xXh Numeric vector expecting the pre-computed h-scaled differences
#'            \eqn{(x - X_1)/h, \ldots, (x - X_n)/h} where \eqn{x} is the
#'            single (!) location for which the weights are to be computed,
#'            the \eqn{X_i}'s are the data values, and \eqn{h} is the numeric
#'            bandwidth scalar.
#' @param thetaXh Numeric vector expecting the pre-computed h-scaled differences
#'                \eqn{(\theta - X_1)/h, \ldots, (\theta - X_n)/h} where
#'                \eqn{\theta} is the numeric scalar location parameter, and the
#'                \eqn{X_i}'s and \eqn{h} are as in \code{xXh}.
#' @param K A kernel function (with vectorized in- & output) to be used for the
#'          estimator.
#' @param h Numeric scalar for bandwidth \eqn{h} (as ``contained'' in
#'          \code{thetaXh} and \code{xXh}).
#'
#' @return If \code{length(sigma)} > 1 a numeric matrix of the dimension
#'         \code{length(sigma)} by \code{length(xXh)} with elements
#'         \eqn{(W_{ni}(x; h, \theta, \sigma_r))} for \eqn{r = 1, \ldots,}
#'         \code{length(sigma)} and \eqn{i = 1, \ldots,} \code{length(xXh)};
#'         otherwise a numeric vector of the same length as \code{xXh}.
#'
#' @references Eichner & Stute (2012) and Eichner (2017): see \link{kader}.
#'
#' @seealso \code{\link{bias_ES2012}} and \code{\link{var_ES2012}} which both
#'          call this function, and \code{\link{kare}} which currently does
#'          the pre-computing.
#'
weights_ES2012 <- function(sigma, xXh, thetaXh, K, h) {
  if (length(xXh) != length(thetaXh))
    stop("Lengths of xXh and thetaXh must be equal!")
  A <- outer(sigma/h, xXh)       # length(sigma) x n
  A <- outer(A, thetaXh, "+")    # length(sigma) x n x n
  A <- rowSums(K(A), dims = 2)   # length(sigma) x n
  drop(A / rowSums(A)) # (W_{ni}(x; \sigma_r))_{1<=r<=length(sigma), 1<=i<=n}
  }



#' Bias Estimator of Eichner & Stute (2012)
#'
#' Bias estimator \eqn{Bias_n(\sigma)}, vectorized in \eqn{\sigma}, on p. 2540
#' of Eichner & Stute (2012).
#'
#' The formula can also be found in eq. (15.21) of Eichner (2017).
#' Pre-computed \eqn{(x - X_i)/h}, \eqn{(\theta - X_i)/h}, and
#' \eqn{m_n(X_i) - m_n(x)} are expected for efficiency reasons (and are
#' currently prepared in function \code{kare}).
#'
#' @export
#'
#' @inheritParams weights_ES2012
#' @param mmDiff Numeric vector expecting the pre-computed differences
#'               \eqn{m_n(X_1) - m_n(x), \ldots, m_n(X_n) - m_n(x)}.
#'
#' @return A numeric vector of the length of \code{sigma}.
#'
#' @references Eichner & Stute (2012) and Eichner (2017): see \link{kader}.
#'
#' @seealso \code{\link{kare}} which currently does the pre-computing.
#'
bias_ES2012 <- function(sigma, h, xXh, thetaXh, K, mmDiff) {
  drop(weights_ES2012(sigma, xXh, thetaXh, K, h) %*% mmDiff) # length(sigma) x 1
 }



#' Variance Estimator of Eichner & Stute (2012)
#'
#' Variance estimator \eqn{Var_n(\sigma)}, vectorized in \eqn{\sigma}, on p.
#' 2540 of Eichner & Stute (2012).
#'
#' The formula can also be found in eq. (15.22) of Eichner (2017).
#' Pre-computed \eqn{(x - X_i)/h}, \eqn{(\theta - X_i)/h}, and
#' \eqn{(Y_i - m_n(x))^2} are expected for efficiency reasons (and are
#' currently prepared in function \code{kare}).
#'
#' @export
#'
#' @inheritParams weights_ES2012
#' @param YmDiff2 Numeric vector of the pre-computed squared differences
#'                \eqn{(Y_1 - m_n(x))^2, \ldots, (Y_n - m_n(x))^2}.
#'
#' @return A numeric vector of the length of \code{sigma}.
#'
#' @references Eichner & Stute (2012) and Eichner (2017): see \link{kader}.
#'
#' @seealso \code{\link{kare}} which currently does the pre-computing.
#'
var_ES2012 <- function(sigma, h, xXh, thetaXh, K, YmDiff2) {
 wx <- weights_ES2012(sigma, xXh, thetaXh, K, h)   # length(sigma) x n
 drop((wx * wx) %*% YmDiff2)   # length(sigma) x 1
 }



#' Kernel Adaptive Regression Estimator
#'
#' Wrapper function which does some preparatory calculations and then calls
#' the actual ``workhorse'' functions which do the main computations for
#' kernel adaptive regression estimation of Eichner & Stute (2012). Finally,
#' it structures and returns the obtained results. Summarizing information
#' and technical details can be found in Eichner (2017).
#'
#' @export
#'
#' @param x.points Vector of location(s) at which the regression estimate is
#'                 to be computed.
#' @param data Data frame or list with one component named \code{x} which
#'             contains the vector of regressor values \eqn{x_1, \ldots, x_n}
#'             and one named \code{y} which holds the vector of pertaining
#'             response values \eqn{y_1, \ldots, y_n} (in the corresponding
#'             order) of the data from which the estimate is to be computed
#'             at the values given in \code{x.points}. Pairs \eqn{(x_i, y_i)}
#'             with \code{NA} or an infinite value in a least one of their
#'             elements are removed (and a warning is issued).
#' @param kernel A character string naming the kernel to be used for the
#'               adaptive estimator. This must partially match one of
#'               "gaussian", "rectangular" or "epanechnikov", with default
#'               "gaussian", and may be abbreviated to a unique prefix.
#'               (Currently, this kernel is also used for the initial,
#'               non-adaptive Nadaraya-Watson regression estimator which enters
#'               into the estimators of bias and variance as described in the
#'               references.)
#' @param Sigma Vector of value(s) of the scale parameter \eqn{\sigma}.
#'              If of length 1 no adaptation is performed. Otherwise
#'              considered as the grid over which the optimization of the
#'              adaptive method will be performed. Defaults to
#'              \code{seq(0.01, 10, length = 51)}.
#' @param h Numeric scalar for bandwidth \eqn{h}. Defaults to NULL and is then
#'          internally set to \eqn{n^{-1/5}}.
#' @param theta Numeric scalar for value of location parameter \eqn{\theta}.
#'              Defaults to NULL and is then internally set to the arithmetic
#'              mean of \eqn{x_1, \ldots, x_n}.
#'
#' @return If \code{length(x.points)} = 1, a list of eight components with the
#'  following names and meanings:
#'  \tabular{ll}{
#'   \code{x} \tab Scalar \eqn{x}-value in \code{x.points} at which the
#'                 regression estimator was computed. \cr
#'   \code{y} \tab Estimated scalar value of \eqn{m(x)} at point in
#'                 \code{x.points}. \cr
#'   \code{sigma.adap} \tab The found scalar minimizer of the MSE-estimator,
#'                          i.e., the adaptive smoothing parameter value. \cr
#'   \code{msehat.min} \tab The found scalar minimum of the MSE-estimator. \cr
#'   \code{Sigma} \tab Vector with the \eqn{\sigma}-grid on which the
#'                     minimization process was performed. \cr
#'   \code{Bn} \tab Vector with the estimator of bias on that
#'                  \eqn{\sigma}-grid. \cr
#'   \code{Vn2} \tab Ditto for the variance. \cr
#'   \code{MSE} \tab Ditto for the MSE. \cr
#'  }
#'
#'  If \code{length(x.points)} > 1, a list with the same component names as
#'  above, but then
#'  \tabular{ll}{
#'   \code{x} \tab Vector \code{x.points} with x-values at which the regression
#'                 estimator was computed. \cr
#'   \code{y} \tab Vector of estimated values of \eqn{m(x)} at the x-values in
#'                 \code{x.points}. \cr
#'   \code{sigma.adap} \tab Vector of the found minimizers of the MSE-estimator,
#'                          i.e., the adaptive smoothing parameter values. \cr
#'   \code{msehat.min} \tab Vector of the found minima of the MSE-estimator. \cr
#'   \code{Sigma} \tab Vector with the \eqn{\sigma}-grid on which the
#'                     minimization process was performed. \cr
#'   \code{Bn} \tab (\code{length(Sigma)} by \code{length(x.points)})-matrix
#'                  with the estimated values of the bias on the
#'                  \eqn{\sigma}-grid in their columns (which correspond to the
#'                  x-values). \cr
#'   \code{Vn2} \tab Ditto for the variance. \cr
#'   \code{MSE} \tab Ditto for the MSE. \cr
#'   }
#'
#' @references Eichner & Stute (2012) and Eichner (2017): see \link{kader}.
#'
#' @examples
#' require(stats)
#'
#'  # Simulation parameters and data generation
#'  #******************************************
#'  # Regression function: a polynomial of degree 4 with one maximum (or
#'  # minimum), one point of inflection, and one saddle point.
#'  # Memo: for p(x) = a * (x - x1) * (x - x2)^3 + b the max. (or min.)
#'  # is at x = (3*x1 + x2)/4, the point of inflection is at x =
#'  # (x1 + x2)/2, and the saddle point at x = x2.
#' m <- function(x, x1 = 0, x2 = 8, a = 0.01, b = 0) {
#'  a * (x - x1) * (x - x2)^3 + b
#' }
#'  # Note: for m()'s default values a minimum is at x = 2, a point of
#'  # inflection at x = 4, and a saddle point an x = 8; an "arbitrary"
#'  # point would, e.g., be at x = 5.
#'
#' x0 <- 5   # The point x_0 at which the MSE-optimal kernel adjusted
#'  # nonparametric estimation of m should take place. (Recall: for m's
#'  # default values a minimum is at 2, a point of inflection at 4, and
#'  # a saddle point an 8; an "arbitrary" point would, e.g., be at 5.)
#'
#' n <- 100   # Sample size
#' sdeps <- 1   # Std. dev. of the \epsilon_i: \sqrt(Var(Y|X=x))
#'              # (here: constant in x)
#' design.ctr <- x0 + 0.5   # "centre" and "scale" of the design, i.e.,
#' design.scl <- 1  # in the normal scenario below, expected value and
#'                  # std. dev. of the distribution of the x_i's.
#'
#' set.seed(42)   # to guanrantee reproducibility.
#'
#' x <- rnorm(n, mean = design.ctr, sd = design.scl)   # x_1, ..., x_n
#' Y <- m(x) + rnorm(length(x), sd = sdeps)            # Y_1, ..., Y_n
#' data <- data.frame(x = x, y = Y)
#'
#'  # Computing the kernel adaptive regression estimator values
#'  #**********************************************************
#' x.points <- seq(-3.3 * design.scl, 3.3 * design.scl, length = 101) +
#'    design.ctr  # x-grid on which to draw and estimate the regr. fct. m.
#'
#' Sigma <- seq(0.01, 10, length = 51)   # \sigma-grid for minimization.
#' fit <- kare(x.points = x0, data = data, Sigma = Sigma)
#'
#' \dontrun{
#'  # Grafical display for the current data set
#'  #******************************************
#'  # Storing the curent settings of the graphics device
#'  # and changing its layout for the three plots to come:
#' op <- par(mfrow = c(3, 1), mar = c(3, 3, 2, 0.1)+0.1,
#'    mgp = c(1.5, 0.5, 0), tcl = -0.3, cex.main = 2)
#'
#'  # The scatter plot of the "raw data":
#' plot(y ~ x, data = data, xlim = range(data$x, x.points),
#'    ylim = range(data$y, fit$y, na.rm = TRUE),
#'    main = bquote(n == .(n)), xlab = "x", ylab = "y")
#'
#'  # The "true" regression function m:
#' lines(x.points, m(x.points), lty = 2)
#'
#'  # The MSE-optimal kernel adjusted nonparametric regression estimator
#'  # at x_0, i.e., the point (x_0, \hat m_n(x_0)):
#' points(fit$x, fit$y, col = "red", pch = 4, cex = 2)
#'
#'  # The legend for the "true" regression function m and for the point
#'  # (x_0, \hat m_n(x_0)):
#' legend("topleft", lty = c(2, NA), pch = c(NA, 4),
#'  col = c("black", "red"), bty = "n", cex = 1.2,
#'  legend = c(as.expression(bquote(paste("m  with  ",
#'                                        sigma(paste(Y, "|", X == x))
#'                                  == .(sdeps)))),
#'             as.expression(bquote(paste(hat(m)[n](x[0]), "  at  ",
#'                                        x[0] == .(x0))))))
#'
#'  # Visualizing the estimators of (Bias_x0(sigma))^2 and
#'  # Var_x0(sigma) on the sigma-grid:
#' with(fit,
#'   matplot(Sigma, cbind(Bn*Bn, Vn2), type = "l", lty = 1:2,
#'    col = c("black", "red"), xlab = expression(sigma), ylab = ""))
#'
#'  # The legend for (Bias_x0(sigma))^2 and Var_x0(sigma):
#' legend("topleft", lty = 1:2, col = c("black", "red"), bty = "n",
#'   legend = c(expression(paste(widehat(plain(Bias))[x[0]]^2, (sigma))),
#'              expression(widehat(plain(Var))[x[0]](sigma))),
#'   cex = 1.2)
#'
#'  # Visualizing the estimator of MSE_x0(sigma) on the sigma-grid together
#'  # with the point indicating the detected minimum, and a legend:
#' plot(fit$Sigma, fit$MSE, type = "l",
#'  xlab = expression(sigma), ylab = "")
#' points(fit$sigma.adap, fit$msehat.min, pch = 4, col = "red", cex = 2)
#' legend("topleft", lty = c(1, NA), pch = c(NA, 4),
#'  col = c("black", "red"), bty = "n", cex = 1.2,
#'  legend = c(expression(widehat(plain(MSE))[x[0]](sigma)),
#'             substitute(group("(", list(plain(Minimizer),
#'                                        plain(Minimum)), ")")
#'                          == group("(", list(x, y), ")"),
#'                        list(x = signif(fit$sigma.adap, 4),
#'                             y = signif(fit$msehat.min, 4)))))
#'
#' par(op) # Restoring the previous settings of the graphics device.
#' }
#'
kare <- function(x.points, data, # Someday to be adapted to args. of ksmooth()?
  kernel = c("gaussian", "epanechnikov", "rectangular"),
  Sigma = seq(0.01, 10, length = 51), h = NULL, theta = NULL) {

  data.name <- deparse(substitute(data))
  data <- as.data.frame(data[c("x", "y")])

  if(any(is.na(data))) {
    data <- stats::na.omit(data)
    warning(length(attr(data, "na.action")), " NA(s) in ", data.name,
      " removed.")
  }
  isinf <- rowSums(sapply(data, is.infinite)) > 0
  if(any(isinf)) {
    data <- data[!isinf,]
    warning(sum(isinf), " rows with infinite value(s) in ", data.name,
            " removed.")
  }

  # Select kernel for adaptive estimator:
  kernel <- match.arg(kernel)
  Kadap <- switch(kernel,
    gaussian     = stats::dnorm,
    epanechnikov = epanechnikov,   # kader:::epanechnikov,
    rectangular  = rectangular)    # kader:::rectangular)

  # Select kernel for (non-adaptive) Nadaraya-Watson estimator.
  # Independent choice not yet implemented.
  Knawa <- Kadap

  n <- nrow(data)
  if(is.null(h)) {
    h <- n^(-1/5)
    message("h set to n^(-1/5) with n = ", n, ".")
  }

  if(is.null(theta)) {
    theta <- mean(data$x)
    message("theta set to arithmetic mean of component x in ", data.name, ".")
  }

  x <- data$x
  Y <- data$y

  Res <- lapply(x.points, function(x0) {
   # Computation of the AMSE-optimal kernel adjusted
   # nonparametric regression estimator of m(x0)
   #************************************************

   # Using the Nadaraya-Watson estimators m_n(x_0) and m_n(x_1), ..., m_n(x_n)
   # with the window width h from above to derive quantities for the estimation
   # of Var(sigma) and Bias(sigma) at x_0:
  mnx0     <- nadwat(x0, x, Y, K = Knawa, h)   # m_n(x_0)
  mnx      <- nadwat(x, x, Y, K = Knawa, h)    # m_n(x_i) for i = 1, ..., n
  mnx_mnx0 <- mnx - mnx0            # m_n(x_i) - m_n(x_0) for i = 1, ..., n
  Y_mnx2   <- Y - mnx               # Y_i - m_n(x_i) for i = 1, ..., n
  Y_mnx2   <- Y_mnx2 * Y_mnx2       # (Y_i - m_n(x_i))^2 for i = 1, ..., n

   # Quantities which do not depend on sigma, but will repeatedly be needed
   # in what follows (makes the algorithmn slightly more efficient):
  x0xh <- (x0 - x) / h;     thetaxh <- (theta - x) / h

   # Estimators of Var_x0(sigma) and Bias_x0(sigma) on the sigma-grid (to
   # also be able to draw them for visualisation purposes later):
  Bn <- bias_ES2012(sigma = Sigma, h = h, xXh = x0xh, thetaXh = thetaxh,
                    K = Kadap, mmDiff = mnx_mnx0)
  Vn2 <- var_ES2012(sigma = Sigma, h = h, xXh = x0xh, thetaXh = thetaxh,
                    K = Kadap, YmDiff2 = Y_mnx2)

   # Composing the estimator of MSE_x0(sigma) on the sigma-grid using the
   # estimators of Var_x0(sigma) and Bias_x0(sigma) on the sigma-grid (as
   # on p. 2540 of Eichner & Stute (2012) (or below eq. (15.22) of Eichner
   # (2017)):
  MSE <- Vn2 + Bn * Bn

   # Minimizer and minimum of the estimator of AMSE_x0(sigma) using "grid
   # search" on the sigma-grid:
  min_idx <- which.min(MSE)
  sig_opt_x0 <- Sigma[min_idx]
  minMSE <- MSE[min_idx]

   # The weights of the AMSE-optimal kernel adjusted nonparametric regression
   # estimator at x_0:
  W <- weights_ES2012(sigma = sig_opt_x0, h = h, xXh = x0xh, thetaXh = thetaxh,
    K = Kadap)

   # And finally, the AMSE-optimal kernel adjusted nonparametric regression
   # estimator at x_0, i.e., \hat m_n(x_0):
  y <- drop(W %*% Y)

  return(list(y = y, x = x0, Bn = Bn, Vn2 = Vn2, MSE = MSE,
              sigma.adap = sig_opt_x0, msehat.min = minMSE))
  } )

  Res <- if (is.list(Res) && length(Res) == 1) {
    Res[[1]]
  } else {
    res1 <- sapply(Res, function(x)
                         unlist(x[c("x", "y", "sigma.adap", "msehat.min")]))
    res1 <- as.list(as.data.frame(t(res1)))
    estimators <- c("Bn", "Vn2", "MSE")
    res2 <- lapply(estimators, function(v) sapply(Res, "[[", v))
    names(res2) <- estimators
    c(res1, res2)
  }
  Res <- c(Res, list(Sigma = Sigma))
  return(Res[c("x", "y", "sigma.adap", "msehat.min", "Sigma", "Bn", "Vn2",
    "MSE")])
}
