#2345678901234567890123456789012345678901234567890123456789012345678901234567890
### ~/PapersBooksReportsReviewsTheses/PapersSelbst/FestschriftWS70/kader/R/
###  kde_adaptive.R
### Functions for computing the adaptive density estimator for "Kernel
### adjusted density estimation" of Srihera & Stute (2011) and for "Rank
### Transformations in Kernel Density Estimation" of Eichner & Stute (2013).
### R 3.4.2, 8./13./14./15./16./17./24.2./21./22./23./24./25.8./27./28.9./4.10.2017
###  (6./7./10.2.2015 / 21./24./26./28./31.10./2./4./8./9./18.11./5.12.2016)
###*****************************************************************************

#' ``Unified'' Function for Kernel Adaptive Density Estimators
#'
#' ``Unified'' function to compute the kernel density estimator both of Srihera
#' & Stute (2011) and of Eichner & Stute (2013).
#'
#' Implementation of both eq. (1.6) in Srihera & Stute (2011) for given and
#' fixed scalars \eqn{\sigma} and \eqn{\theta}, and eq. (4) in Eichner & Stute
#' (2013) for a given and fixed scalar \eqn{\sigma} and for a given and fixed
#' rank transformation (and, of course, for fixed and given location(s) in
#' \eqn{x}, data \eqn{(X_1, \ldots, X_n)}, a kernel function \eqn{K} and a
#' bandwidth \eqn{h}). The formulas that the computational version implemented
#' here is based upon are given in eq. (15.3) and eq. (15.9), respectively, of
#' Eichner (2017). This function rests on preparatory computations done in
#' \code{\link{fnhat_SS2011}} or \code{\link{fnhat_ES2013}}.
#'
#' @export
#'
#' @param x Numeric vector with the location(s) at which the density estimate
#'          is to be computed.
#' @param data Numeric vector \eqn{(X_1, \ldots, X_n)} of the data from which
#'             the estimate is to be computed.
#' @param K A kernel function (with vectorized in- & output) to be used for
#'          the estimator.
#' @param h Numeric scalar for bandwidth \eqn{h}.
#' @param Bj Numeric vector expecting \eqn{(-J(1/n), \ldots, -J(n/n))} as
#'           produced in \code{\link{fnhat_SS2011}} in case of the rank
#'           transformation method (using an admissible rank transformation
#'           as implemented by \code{\link{J_admissible}}), but
#'           \eqn{(\hat \theta - X_1}, \ldots, \eqn{\hat \theta - X_n)} as produced
#'           in \code{\link{fnhat_ES2013}} in case of the non-robust method.
#' @param sigma Numeric scalar for value of scale parameter \eqn{\sigma}.
#'
#' @return A numeric vector of the same length as \code{x} with the estimated
#'         density values from eq. (1.6) of Srihera & Stute (2011) or eq. (4)
#'         of Eichner & Stute (2013).
#'
#' @note In case of the rank transformation method the data are expected to
#'       be sorted in increasing order.
#'
#' @references Srihera & Stute (2011), Eichner and Stute (2013), and Eichner
#'             (2017): see \code{\link{kader}}.
#'
#' @examples
#' require(stats)
#'
#'  # The kernel density estimators for simulated N(0,1)-data and a single
#'  # sigma-value evaluated on a grid using the rank transformation and
#'  # the non-robust method:
#' set.seed(2017);     n <- 100;     Xdata <- rnorm(n)
#' xgrid <- seq(-4, 4, by = 0.1)
#' negJ <- -J_admissible(1:n / n)                 # The rank trafo requires
#' compute_fnhat(x = xgrid, data = sort(Xdata),   # sorted data!
#'   K = dnorm, h = n^(-1/5), Bj = negJ, sigma = 1)
#'
#' theta.X <- mean(Xdata) - Xdata    # non-robust method
#' compute_fnhat(x = xgrid, data = Xdata, K = dnorm, h = n^(-1/5),
#'   Bj = theta.X, sigma = 1)
#'
compute_fnhat <- function(x, data, K, h, Bj, sigma) {
  sig_hh <- sigma / (h * h)
  Bj_h <- Bj/h
  sig_hh * sapply(x, function(x0) {
    sig..x.X_hh <- sig_hh * (x0 - data)   # = A_i, length(data) x 1
    Karg <- outer(sig..x.X_hh, Bj_h, "+") # length(data) x length(data)
    # Delete diagonal to average only over pairs of different indices
    # (cf. eq. (1.6) in S. & S. (2011) or eq. (4) in E. & S. (2013)):
    Karg <- Karg[ c(FALSE, rep(TRUE, nrow(Karg)))]
    mean(K(Karg))
  })   # length(x) x 1
}



#' (Non-robust) Kernel Density Estimator of Srihera & Stute (2011)
#'
#' Implementation of eq. (1.6) in Srihera & Stute (2011) for given and fixed
#' scalars \eqn{\sigma} and \eqn{\theta} (and, of course, for fixed and given
#' location(s) in \eqn{x}, data \eqn{(X_1, \ldots, X_n)}, a kernel function
#' \eqn{K} and a bandwidth \eqn{h}).
#'
#' The formula upon which the computational version implemented here is based
#' is given in eq. (15.3) of Eichner (2017). This function does mainly only a
#' simple preparatory computation and then calls \code{\link{compute_fnhat}}
#' which does the actual work.
#'
#' @export
#'
#' @param x Numeric vector with the location(s) at which the density estimate
#'          is to be computed.
#' @param data Numeric vector \eqn{(X_1, \ldots, X_n)} of the data from which
#'             the estimate is to be computed. Missing or infinite values are
#'             not allowed and entail an error.
#' @param K A kernel function to be used for the estimator.
#' @param h Numeric scalar for bandwidth \eqn{h}.
#' @param theta Numeric scalar for value of location parameter \eqn{\theta}.
#' @param sigma Numeric scalar for value of scale parameter \eqn{\sigma}.
#'
#' @return An object with class "density" whose underlying structure is
#'         a list containing the following components (as described in
#'         \code{\link[stats]{density}}), so that the \code{print} and
#'         \code{plot} methods for \code{density}-objects are
#'         immediately available):
#' \tabular{ll}{
#'  \code{x}  \tab the n coordinates of the points where the density is
#'                 estimated. \cr
#'  \code{y}  \tab the estimated density values from eq. (1.6) in Srihera &
#'                 Stute (2011). \cr
#'  \code{bw} \tab the bandwidth used. \cr
#'  \code{n}  \tab the sample size. (Recall: missing or infinite values are
#'                 not allowed here.) \cr
#'  \code{call}      \tab the call which produced the result. \cr
#'  \code{data.name} \tab the deparsed name of the x argument. \cr
#'  \code{has.na}    \tab logical, for compatibility (always FALSE). \cr\cr
#'  Additionally: \tab \cr
#'  \code{theta} \tab as in Arguments. \cr
#'  \code{sigma} \tab as in Arguments. \cr
#'  }
#'
#' @seealso \code{\link{fnhat_ES2013}}.
#'
#' @references Srihera & Stute (2011) and Eichner (2017): see \link{kader}.
#'
#' @examples
#' require(stats);   require(grDevices);    require(datasets)
#'
#'  # Simulated N(0,1)-data and one sigma-value
#' set.seed(2017);     n <- 100;     d <- rnorm(n)
#' xgrid <- seq(-4, 4, by = 0.1)
#' (fit <- fnhat_SS2011(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'   theta = mean(d), sigma = 1))
#' \donttest{
#' plot(fit, ylim = range(0, dnorm(0), fit$y), col = "blue")
#' curve(dnorm, add = TRUE);   rug(d, col = "red")
#' legend("topleft", lty = 1, col = c("blue", "black", "red"),
#'   legend = expression(tilde(f)[n], phi, "data")) }
#' \donttest{
#'  # The same data, but several sigma-values
#' sigmas <- seq(1, 4, length = 4)
#' (fit <- lapply(sigmas, function(sig)
#'   fnhat_SS2011(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'     theta = mean(d), sigma = sig)))
#'
#' ymat <- sapply(fit, "[[", "y")
#' matplot(x = xgrid, y = ymat, type = "l", lty = 1, col = 3:6,
#'   ylim = range(0, dnorm(0), ymat), main = "", xlab = "", ylab = "Density")
#' curve(dnorm, add = TRUE);   rug(d, col = "red")
#' legend("topleft", lty = 1, col = c("black", "red", NA), bty = "n",
#'   legend = expression(phi, "data", tilde(f)[n]~"in other colors")) }
#' \donttest{
#'  # Old-Faithful-eruptions-data and several sigma-values
#' d <- faithful$eruptions;     n <- length(d);     er <- extendrange(d)
#' xgrid <- seq(er[1], er[2], by = 0.1);    sigmas <- seq(1, 4, length = 4)
#' (fit <- lapply(sigmas, function(sig)
#'    fnhat_SS2011(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'      theta = mean(d), sigma = sig)))
#'
#' ymat <- sapply(fit, "[[", "y");     dfit <- density(d, bw = "sj")
#' plot(dfit, ylim = range(0, dfit$y, ymat), main = "", xlab = "")
#' rug(d, col = "red")
#' matlines(x = xgrid, y = ymat, lty = 1, col = 3:6)
#' legend("top", lty = 1, col = c("black", "red", NA), bty = "n",
#'   legend = expression("R's est.", "data", tilde(f)[n]~"in other colors")) }
fnhat_SS2011 <- function(x, data, K, h, theta, sigma) {
  if(any(is.na(data) | is.infinite(data)))
    stop("Missing or infinite values in data not allowed!")
  if(length(sigma) != 1) stop("sigma must have length 1!")

  theta.X <- theta - data   # = B_j, length(data) x 1

  y <- compute_fnhat(x = x, data = data, K = K, h = h, Bj = theta.X,
    sigma = sigma)   # length(x) x 1

  res <- list(x = x, y = y, bw = h, n = length(data),
    call = match.call(), data.name = deparse(substitute(data)),
    has.na = FALSE, theta = theta, sigma = sigma)
  class(res) <- "density"
  return(res)
  }




#' Robust Kernel Density Estimator of Eichner & Stute (2013)
#'
#' Implementation of eq. (4) in Eichner & Stute (2013) for a given and fixed
#' scalar \eqn{\sigma}, for rank transformation function \eqn{J} (and, of
#' course, for fixed and given location(s) in \eqn{x}, data \eqn{(X_1, \ldots,
#' X_n)}, a kernel function \eqn{K}, and a bandwidth \eqn{h}).
#'
#' The formula upon which the computational version implemented here is based
#' is given in eq. (15.9) of Eichner (2017). This function does mainly only a
#' simple preparatory computation and then calls \code{\link{compute_fnhat}}
#' which does the actual work.
#'
#' @export
#'
#' @param x Numeric vector with the location(s) at which the density
#'          estimate is to be computed.
#' @param data Numeric vector \eqn{(X_1, \ldots, X_n)} of the data from
#'             which the estimate is to be computed. Missing values are not
#'             allowed and entail an error.
#' @param K A kernel function to be used for the estimator.
#' @param h Numeric scalar for bandwidth \eqn{h}.
#' @param ranktrafo A function used for the rank transformation.
#' @param sigma Numeric scalar for value of scale parameter \eqn{\sigma}.
#'
#' @return An object with class "density" whose underlying structure is
#'         a list containing the following components (as described in
#'         \code{\link[stats]{density}}), so that the \code{print} and
#'         \code{plot} methods for \code{density}-objects are
#'         immediately available):
#' \tabular{ll}{
#'  \code{x}  \tab the n coordinates of the points where the density is
#'                 estimated. \cr
#'  \code{y}  \tab the estimated density values from eq. (4) in Eichner & Stute
#'                 (2013). \cr
#'  \code{bw} \tab the bandwidth used. \cr
#'  \code{n}  \tab the sample size. (Recall: missing or infinite values are not
#'                 allowed here.) \cr
#'  \code{call}      \tab the call which produced the result. \cr
#'  \code{data.name} \tab the deparsed name of the x argument. \cr
#'  \code{has.na}    \tab logical, for compatibility (always FALSE). \cr\cr
#'  Additionally: \tab \cr
#'  \code{ranktrafo} \tab as in Arguments. \cr
#'  \code{sigma} \tab as in Arguments. \cr
#'  }
#'
#' @seealso \code{\link{fnhat_SS2011}}.
#'
#' @references Eichner & Stute (2013) and Eichner (2017): see \link{kader}.
#'
#' @examples
#' require(stats);   require(grDevices);   require(datasets)
#'
#'  # Simulated N(0,1)-data and one sigma-value
#' set.seed(2016);     n <- 100;     d <- rnorm(n)
#' xgrid <- seq(-4, 4, by = 0.1)
#' (fit <- fnhat_ES2013(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'   ranktrafo = J2, sigma = 1) )
#' \donttest{
#' plot(fit, ylim = range(0, dnorm(0), fit$y), col = "blue")
#' curve(dnorm, add = TRUE);   rug(d, col = "red")
#' legend("topleft", lty = 1, col = c("blue", "black", "red"),
#'   legend = expression(hat(f)[n], phi, "data"))
#' }
#' \donttest{
#'  # The same data, but several sigma-values
#' sigmas <- seq(1, 4, length = 4)
#' (fit <- lapply(sigmas, function(sig)
#'   fnhat_ES2013(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'     ranktrafo = J2, sigma = sig)) )
#'
#' ymat <- sapply(fit, "[[", "y")
#' matplot(x = xgrid, y = ymat, type = "l", lty = 1, col = 2 + seq(sigmas),
#'   ylim = range(0, dnorm(0), ymat), main = "", xlab = "", ylab = "Density")
#' curve(dnorm, add = TRUE);   rug(d, col = "red")
#' legend("topleft", lty = 1, col = c("black", "red", NA), bty = "n",
#'   legend = expression(phi, "data", hat(f)[n]~"in other colors"))
#' }
#' \donttest{
#'  # Old-Faithful-eruptions-data and several sigma-values
#' d <- faithful$eruptions;     n <- length(d);     er <- extendrange(d)
#' xgrid <- seq(er[1], er[2], by = 0.1);    sigmas <- seq(1, 4, length = 4)
#' (fit <- lapply(sigmas, function(sig)
#'    fnhat_ES2013(x = xgrid, data = d, K = dnorm, h = n^(-1/5),
#'      ranktrafo = J2, sigma = sig)) )
#'
#' ymat <- sapply(fit, "[[", "y");     dfit <- density(d, bw = "sj")
#' plot(dfit, ylim = range(0, dfit$y, ymat), main = "", xlab = "")
#' rug(d, col = "red")
#' matlines(x = xgrid, y = ymat, lty = 1, col = 2 + seq(sigmas))
#' legend("top", lty = 1, col = c("black", "red", NA), bty = "n",
#'   legend = expression("R's est.", "data", hat(f)[n]~"in other colors"))
#' }
fnhat_ES2013 <- function(x, data, K, h, ranktrafo, sigma) {
  if(any(is.na(data) | is.infinite(data)))
    stop("Missing or infinite values in data not allowed!")
  if(length(sigma) != 1) stop("sigma must have length 1!")

  data.name <- deparse(substitute(data))

  if (is.unsorted(data)) data <- sort(data)
  n <- length(data)
  negRT <- -ranktrafo(1:n / n)   # = B_j, n x 1. This is why
                                 # the data *must* be sorted!
  y <- compute_fnhat(x = x, data = data, K = K, h = h, Bj = negRT,
    sigma = sigma)   # length(x) x 1

  res <- list(x = x, y = y, bw = h, n = n, call = match.call(),
    data.name = data.name, has.na = FALSE, ranktrafo = ranktrafo,
    sigma = sigma)
  class(res) <- "density"
  return(res)
  }



#' Specialized ``Workhorse'' Function for Kernel Adaptive Density Estimators
#'
#' Common specialized computational ``workhorse'' function to compute the kernel
#' adaptive density estimators both in eq. (1.6) of Srihera & Stute (2011) and
#' in eq. (4) of Eichner & Stute (2013) (together with several related
#' quantities) with a \eqn{\sigma} that minimizes the estimated MSE using an
#' estimated \eqn{\theta}. This function is ``specialized'' in that it expects
#' some pre-computed quantities (in addition to the point(s) at which the
#' density is to be estimated, the data, etc.). In particular, the estimator of
#' \eqn{\theta} (which is typically the arithmetic mean of the data) is
#' expected to be already ``contained'' in those pre-computed quantities, which
#' increases the computational efficiency.
#'
#' The computational procedure in this function can be highly iterative because
#' for each point in \code{x} (and hence for each row of matrix \code{Ai}) the
#' MSE estimator is computed as a function of \eqn{\sigma} on a (usually fine)
#' \eqn{\sigma}-grid provided through \code{sigma}. This happens by repeated
#' calls to \code{\link{bias_AND_scaledvar}()}. The minimization in \eqn{\sigma}
#' is then performed by \code{\link{minimize_MSEHat}()} using both a discrete
#' grid-search and the numerical optimization routine implemented in base R's
#' \code{optimize()}. Finally, \code{\link{compute_fnhat}()} yields the actual
#' value of the density estimator for the adapted \eqn{\sigma}, i.e., for the
#' MSE-estimator-minimizing \eqn{\sigma}.
#' (If necessary the computation over the \eqn{\sigma}-grid is repeated after
#' extending the range of the grid until the estimator functions for both bias
#' and variance are \emph{not constant} across the \eqn{\sigma}-grid.)
#'
#' @export
#'
#' @param x Numeric vector \eqn{(x_1, \ldots, x_k)} of location(s) at which the
#'          density estimate is to be computed.
#' @param data Numeric vector \eqn{(X_1, \ldots, X_n)} of the data from which
#'             the estimate is to be computed.
#' @param K Kernel function with vectorized in- & output.
#' @param h Numeric scalar, where (usually) \eqn{h = n^{-1/5}}.
#' @param sigma Numeric vector \eqn{(\sigma_1, \ldots, \sigma_s)} with
#'              \eqn{s \ge 1}.
#' @param Ai Numeric matrix expecting in its i-th row \eqn{(x_i - X_1, \ldots,
#'           x_i - X_n)/h}, where (usually) \eqn{x_1, \ldots, x_k} with
#'           \eqn{k =} \code{length(x)} are the points at which the density is
#'           to be estimated for the data \eqn{X_1, \ldots, X_n} with
#'           \eqn{h = n^{-1/5}}.
#' @param Bj Numeric vector expecting \eqn{(-J(1/n), \ldots, -J(n/n))} in
#'           case of the rank transformation method, but \eqn{(\hat \theta -
#'           X_1, \ldots, \hat \theta - X_n)} in case of the non-robust
#'           Srihera-Stute-method.
#' @param fnx Numeric vector expecting \eqn{(f_n(x_1), \ldots, f_n(x_k))} with
#'            \eqn{f_n(x_i)} the Parzen-Rosenblatt estimator at \eqn{x_i}, i.e.,
#'            \eqn{f_n(x_i) =} \code{mean(K(Ai[i,]))/h} where here typically
#'            \code{h} \eqn{= n^{-1/5}}.
#' @param ticker Logical; determines if a 'ticker' documents the iteration
#'               progress through \code{sigma}. Defaults to FALSE.
#' @param plot Logical or character or numeric and indicates if graphical
#'             output should be produced. Defaults to \code{FALSE} (i.e., no
#'             graphical output is produced). If it is a character string or
#'             a numeric value, graphical output will be written to numbered
#'             pdf-files (one for each element of \code{x}, in the current
#'             working directory) whose names start with the provided
#'             ``value'' after converting it into a character string
#'             followed by the index number of the pertaining
#'             \code{x}-element. (Parts of the graphical output are
#'             generated by \code{\link{minimize_MSEHat}}.)
#' @param parlist A list of graphical parameters; affects only the pdf-files
#'               (if any are created at all). Default: \code{NULL}.
#' @param \ldots Possible further arguments passed to \code{minimize_MSEHat()}
#'               (where they are currently ignored).
#'
#' @return A list of as many lists as elements in \code{x}, each with components
#'   \code{x}, \code{y}, \code{sigma.adap}, \code{msehat.min},
#'   \code{discr.min.smaller}, and \code{sig.range.adj} whose meanings are as
#'   follows:
#'   \tabular{ll}{
#'    \code{x} \tab the n coordinates of the points where the density is
#'             estimated. \cr
#'    \code{y} \tab the estimate of the density value \eqn{f(x)}. \cr
#'    \code{sigma.adap} \tab Minimizer of MSE-estimator (from function
#'                      \code{\link{minimize_MSEHat}}). \cr
#'    \code{msehat.min} \tab Minimum of MSE-estimator (from function
#'                      \code{\link{minimize_MSEHat}}). \cr
#'    \code{discr.min.smaller} \tab TRUE iff the numerically found minimum was
#'                             smaller than the discrete one (from function
#'                             \code{\link{minimize_MSEHat}}). \cr
#'    \code{sig.range.adj} \tab Number of adjustments of sigma-range. \cr
#'    }
#'
#' @references Srihera & Stute (2011) and Eichner & Stute (2013): see
#'             \link{kader}.
#'
#' @examples
#' \dontrun{
#' require(stats)
#'
#'  # Kernel adaptive density estimators for simulated N(0,1)-data
#'  # computed on an x-grid using the rank transformation and the
#'  # non-robust method:
#' set.seed(2017);     n <- 100;     Xdata <- sort(rnorm(n))
#' x <- seq(-4, 4, by = 0.5);     Sigma <- seq(0.01, 10, length = 51)
#' h <- n^(-1/5)
#'
#' x.X_h <- outer(x/h, Xdata/h, "-")
#' fnx <- rowMeans(dnorm(x.X_h)) / h   # Parzen-Rosenblatt estim. at
#'                                     # x_j, j = 1, ..., length(x).
#'  # non-robust method:
#' theta.X <- mean(Xdata) - Xdata
#' adaptive_fnhat(x = x, data = Xdata, K = dnorm, h = h, sigma = Sigma,
#'   Ai = x.X_h, Bj = theta.X, fnx = fnx, ticker = TRUE, plot = TRUE)
#'
#'  # rank transformation-based method (requires sorted data):
#' negJ <- -J_admissible(1:n / n)   # rank trafo
#' adaptive_fnhat(x = x, data = Xdata, K = dnorm, h = h, sigma = Sigma,
#'   Ai = x.X_h, Bj = negJ, fnx = fnx, ticker = TRUE, plot = TRUE)
#' }
#'
adaptive_fnhat <- function(x, data, K, h, sigma, Ai, Bj, fnx, ticker = FALSE,
  plot = FALSE, parlist = NULL, ...) {
  # For the following quantities on the lhs of each "=" see kare():
  # Ai = x.X_h.matrix in both methods,
  # Bj = theta.X for the non-robust method of S. & S. (2011), and
  # Bj = negRT for the robust method of E. & S. (2013).

  message("\nFor each element in x: Computing estimated values of\n",
          "bias and scaled variance on the sigma-grid.")
  message("Note: x has ", length(x), " element(s) and the sigma-grid ",
          length(sigma), ".\n")
  nxns <- length(x) * length(sigma)
  if(nxns > 500) {
    message("This may need more time than enough to get yourself\n",
      "a tea or a coffee instead of watching the 'ticker'\n", appendLF = FALSE)
   } else if(nxns > 50) {
    message("Please, have a little patience; this may need some time.\n",
      "You may want to watch the 'ticker' as time goes by\n", appendLF = FALSE)
   } else message("As a little distraction, the 'ticker' documents the\n",
      "computational progress", appendLF = FALSE)
  message(" (if you have set ticker = TRUE).")

  if(is.character(plot) || is.numeric(plot)) {
    prefix <- as.character(plot)
    idxformatwidth <- ceiling(log10(length(x)+0.1))
    plot <- TRUE
  } else prefix <- NULL

  lapply(seq_along(x),
    function(i) {   # For ith element in x, i.e., for ith row of matrix Ai:
       # Computation of the MSE estimator (MSEHat) as a function of \sigma for
       # values on a - usually fine - \sigma-grid provided by sigma.
      if(ticker) message("x[", i, "]:", appendLF = FALSE)
      sigma.range.adjusted <- 0
      repeat {   # If necessary the computation over the \sigma-grid is
         # repeated after extending the range of the grid until both
         # estimator functions are *not constant* across the grid.

        BV <- bias_AND_scaledvar(sigma = sigma, Ai = Ai[i,], Bj = Bj,
          h = h, K = K, fnx = fnx[i], ticker = ticker);   message("")

        VarHat.scaled <- BV$VarHat.scaled
        BiasHat.squared <- BV$BiasHat * BV$BiasHat
        MSEHat <- BiasHat.squared + VarHat.scaled
        if (max(diff(range(BiasHat.squared)),
                diff(range(VarHat.scaled))) > .Machine$double.eps^0.25) break

        warning("Biashat.squared or VarHat.scaled was - almost - constant!")
        sigma <- 10 * sigma
        sigma.range.adjusted <- sigma.range.adjusted + 1
        message("New range for sigma: ", sigma[1], " - ", sigma[length(sigma)])
        }   # end of repeat {....}

      if(plot) { # Note that plot = FALSE by default.
         # Draws the squared bias estimator (Bias_x0(\sigma))^2, the
         # scaled variance estimator \sigma^2/(h^4 n) Var_x0(\sigma),
         # and the MSE estimator on the \sigma-grid.

        if(!is.null(prefix)) {
          grDevices::pdf(paste0(prefix, formatC(i, width = idxformatwidth,
            flag=0), ".pdf"))
          graphics::par(parlist)
        }

        graphics::matplot(sigma, cbind(BiasHat.squared, VarHat.scaled, MSEHat),
          type = "l", lty = c("twodash", "longdash", "solid"), lwd = 2,
          col = c("blue", "violet", "magenta3"), main = paste("x =", x[i]),
          ylim = c(0, max(MSEHat, na.rm = TRUE)),
          xlab = expression(sigma), ylab = NA)

         # Legend for (Bias_x0(\sigma))^2 and Var_x0(\sigma):
        graphics::legend("topleft", lty = c("twodash", "longdash", "solid"),
          lwd = 2, col = c("blue", "violet", "magenta3"), bty = "n",
          legend = c(expression(widehat(plain(Bias))^2~(sigma)),
                     expression(sigma^2/(n*h^4)~widehat(plain(Var))(sigma)),
                     expression(widehat(plain(MSE))(sigma))),
          inset = c(0.05, 0)) # cex = 0.8

        if(sigma.range.adjusted > 0)
          graphics::legend("center", legend = "Range of sigma was extended.",
                            bty = "n", cex = 1.3)
        }   # end of if(plot) {....}

       # Minimizing MSEHat = BiasHat.squared + VarHat.scaled as a function of
       # \sigma using two different methods; see minimize_MSEHat()!
      minMSE <- minimize_MSEHat(VarHat.scaled, BiasHat.squared, sigma = sigma,
        Ai = Ai[i,], Bj = Bj, h = h, K = K, fnx = fnx[i], ticker = ticker,
        plot = plot, ...);  if(ticker) message("")

      if(!is.null(prefix)) grDevices::dev.off()

       # Computation of kernel density estimator with adapted \sigma
      y <- compute_fnhat(x = x[i], data = data, K = K, h = h, Bj = Bj,
        sigma = minMSE$sigma.adap)

      c(list(x = x[i], y = y), minMSE,
        list(sig.range.adj = sigma.range.adjusted))
    } # end of function(i) {....}
  ) # end of sapply(seq_along(x), ....)
  }




#' Kernel Adaptive Density Estimator
#'
#' Wrapper function which does some preparatory calculations and then calls
#' the actual ``workhorse'' functions which do the main computations for
#' kernel adaptive density estimation of Srihera & Stute (2011) or Eichner
#' & Stute (2013). Finally, it structures and returns the obtained results.
#' Summarizing information and technical details can be found in Eichner
#' (2017).
#
#' @export
#'
#' @param x Vector of location(s) at which the density estimate is
#'          to be computed.
#' @param data Vector \eqn{(X_1, \ldots, X_n)} of the data from which the
#'             estimate is to be computed. \code{NA}s or infinite values are
#'             removed (and a warning is issued).
#' @param kernel A character string naming the kernel to be used for the
#'               adaptive estimator. This must partially match one of
#'               "gaussian", "rectangular" or "epanechnikov", with default
#'               "gaussian", and may be abbreviated to a unique prefix.
#'               (Currently, this kernel is also used for the initial,
#'               non-adaptive Parzen-Rosenblatt estimator which enters
#'               into the estimators of bias and variance as described in the
#'               references.)
#' @param method A character string naming the method to be used for the
#'               adaptive estimator. This must partially match one of
#'               "both", "ranktrafo" or "nonrobust", with default "both",
#'               and may be abbreviated to a unique prefix.
#' @param Sigma Vector of value(s) of the scale parameter \eqn{\sigma}.
#'              If of length 1 no adaptation is performed. Otherwise
#'              considered as the initial grid over which the optimization
#'              of the adaptive method will be performed. Defaults to
#'              \code{seq(0.01, 10, length = 51)}.
#' @param h Numeric scalar for bandwidth \eqn{h}. Defaults to NULL and is then
#'          internally set to \eqn{n^{-1/5}}.
#' @param theta Numeric scalar for value of location parameter \eqn{\theta}.
#'              Defaults to NULL and is then internally set to the arithmetic
#'              mean of \eqn{x_1, \ldots, x_n}.
#' @param ranktrafo Function used for the rank transformation. Defaults to
#'                  \code{\link{J2}} (with its default \code{cc = sqrt(5))}.
#' @param ticker Logical; determines if a 'ticker' documents the iteration
#'               progress through \code{Sigma}. Defaults to FALSE.
#' @param plot Logical or character or numeric and indicates if graphical
#'             output should be produced. Defaults to FALSE (i.e., no
#'             graphical output is produced) and is passed to
#'             \code{\link{adaptive_fnhat}()} which does the actual work. For
#'             details on how it is processed see there.
#' @param parlist A list of graphical parameters that is passed to
#'                \code{\link{adaptive_fnhat}()}; see there. Default:
#'                \code{NULL}.
#' @param \ldots Further arguments possibly passed down. Currently ignored.
#'
#' @return In the case of only one method a data frame whose components
#'  have the following names and meanings:
#'  \tabular{ll}{
#'   \code{x} \tab {x_0}. \cr
#'   \code{y} \tab Estimate of f(x_0). \cr
#'   \code{sigma.adap} \tab The found minimizer of the MSE-estimator, i.e.,
#'                          the adaptive smoothing parameter value. \cr
#'   \code{msehat.min} \tab The found minimum of the MSE-estimator. \cr
#'   \code{discr.min.smaller} \tab TRUE iff the numerically found minimum
#'                                 was smaller than the discrete one. \cr
#'   \code{sig.range.adj} \tab Number of adjustments of sigma-range. \cr
#'   }
#'
#'   In the case of both methods a list of two data frames of the
#'   just described structure.
#'
#' @references Srihera & Stute (2011), Eichner & Stute (2013), and Eichner
#'             (2017): see \code{\link{kader}}.
#'
#' @examples
#' require(stats)
#'
#'  # Generating N(0,1)-data
#' set.seed(2017);     n <- 80;     d <- rnorm(n)
#'
#'  # Estimating f(x0) for one sigma-value
#' x0 <- 1
#' (fit <- kade(x = x0, data = d, method = "nonrobust", Sigma = 1))
#' \donttest{
#'  # Estimating f(x0) for sigma-grid
#' x0 <- 1
#' (fit <- kade(x = x0, data = d, method = "nonrobust",
#'   Sigma = seq(0.01, 10, length = 10), ticker = TRUE))
#' }
#' \dontrun{
#'  # Estimating f(x0) for sigma-grid and Old-Faithful-eruptions-data
#' x0 <- 2
#' (fit <- kade(x = x0, data = faithful$eruptions, method = "nonrobust",
#'   Sigma = seq(0.01, 10, length = 51), ticker = TRUE, plot = TRUE))
#' }
kade <- function(x, data,  # Someday to be adapted to args. of fct. density?
  kernel = c("gaussian", "epanechnikov", "rectangular"),
  method = c("both", "ranktrafo", "nonrobust"),
  Sigma = seq(0.01, 10, length = 51), h = NULL, theta = NULL,
  ranktrafo = J2, ticker = FALSE, plot = FALSE, parlist = NULL,
  # parlist = list(mfrow = c(3, 1), mar = c(2.5, 2, 2, 0.1) + 0.1,
  #                mgp = c(1.5, 0.5, 0), tcl = -0.3, cex = 0.8),
  ...) {

  data.name <- deparse(substitute(data))

  if(any(is.na(data))) {
    data <- stats::na.omit(data)
    warning(length(attr(data, "na.action")), " NA(s) in ", data.name,
            " removed.")
    }
  if(any(is.infinite(data))) {
    data <- data[!is.infinite(data)]
    warning(sum(is.infinite(data)), " infinite value(s) in ", data.name,
      " removed.")
  }

   # Select kernel for adaptive estimator:
  kernel <- match.arg(kernel)
  Kadap <- switch(kernel,
    gaussian     = stats::dnorm,
    epanechnikov = epanechnikov,   # kader:::epanechnikov,
    rectangular  = rectangular)    # kader:::rectangular)

   # Select kernel for (non-adaptive) Parzen-Rosenblatt estimator.
   # Independent choice not yet implemented.
  Kparo <- Kadap

  n <- length(data)
  if(is.null(h)) {
    h <- n^(-1/5)
    message("h set to n^(-1/5) with n = ", n, ".")
    }

  if(is.null(theta)) {
    theta <- mean(data)
    message("theta set to arithmetic mean of data in ", data.name, ".")
    }

  ## Uncomment to open an enlarged graphics screen device:
##  x11(width = 10, height = 10, xpos = 100, ypos = -50)
  # if (plot) {
  #   op <- parlist
  #   graphics::stripchart(data, ...)   # Strip chart of the data.
  #   on.exit(graphics::par(op))
  # }

   # Select method for adaptive estimator:
  method <- match.arg(method)

  message("Using the ", appendLF = FALSE)
  if(length(Sigma) == 1) { # Only one sigma-value, whence no minimization!

    y1 <- if(method %in% c("both", "nonrobust")) {
      message("adaptive method of Srihera & Stute (2011)", appendLF = FALSE)
      fnhat_SS2011(x = x, data = data, K = Kadap, h = h, theta = theta,
                   sigma = Sigma)
      }
    if(method == "both") message("and the")

    y2 <- if(method %in% c("both", "robust")) {
      message("rank transformation-based robust adaptive method of ",
              "Eichner & Stute (2013)", appendLF = FALSE)
      # (Note: sorting the data will happen in fnhat_ES2013().)
      fnhat_ES2013(x = x, data = data, K = Kadap, h = h,
                   ranktrafo = ranktrafo, sigma = Sigma)
      }
    message(".")

    switch(method,
           nonrobust = y1, robust = y2, list(SS2011 = y1, ES2013 = y2))

   } else { # Minimization in Sigma

    # Computation of quantities that do not depend on the scale parameter
    # \sigma and will repeatedly be needed in what follows. Their computation
    # in advance makes the algorithm slightly more efficient.

    # Sorting the data is *absolutely necessary* because the implementation
    # in adaptive_fnhat() uses this fact for efficiency reasons in the use
    # of negRT.
    if(is.unsorted(data)) data <- sort(data)

    x.X_h.matrix <- outer(x/h, data/h, "-")    # rows of A_i's, length(x) x n.
    fnx <- rowMeans(Kparo(x.X_h.matrix)) / h   # Parzen-Rosenblatt estim. at
                                               # x_j, j = 1, ..., length(x).

    y1 <- if(method %in% c("both", "nonrobust")) {
      message("adaptive method of Srihera & Stute (2011)", appendLF = FALSE)
      theta.X <- theta - data   # B_j (theta = arith. mean of data).
      adaptive_fnhat(x = x, data = data, K = Kadap, h = h, sigma = Sigma,
        Ai = x.X_h.matrix, Bj = theta.X, fnx = fnx, ticker = ticker,
        plot = plot, parlist = parlist)
      }
    if(method == "both") message("and the")

    y2 <- if(method %in% c("both", "robust")) {
      message("rank transformation-based robust adaptive method of ",
              "Eichner & Stute (2013)", appendLF = FALSE)
      # Rank transformation-values: This is why the data *must* be sorted!
      negRT <- -ranktrafo(1:n / n)  # B_j
      adaptive_fnhat(x = x, data = data, K = Kadap, h = h, sigma = Sigma,
        Ai = x.X_h.matrix, Bj = negRT, fnx = fnx, ticker = ticker,
        plot = plot, parlist = parlist)
      }
    message(".")

    Y <- switch(method,
                nonrobust = y1, robust = y2, list(SS2011 = y1, ES2013 = y2))
    if(method == "both") {
      lapply(Y, function(y) do.call("rbind", lapply(y, data.frame)))
      } else do.call("rbind", lapply(Y, data.frame))
    } # end of minimization in Sigma
} # end of kade()
