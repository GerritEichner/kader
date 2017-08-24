### ~/Beratung/StuteW/KernelAdjustedNonparRegression/AMSE_Minimierung.R
### Schaetzung und Minimierung der AMSE-Funktion fuer die Simulationen in "Ker-
### nel Adjusted Nonparametric Regression" von Eichner & Stute
### R 2.12.2, 17./18./22.3./5./6./7./8./14./18./26./27./29.4./27.5./3.6.2011
###********************************************************************************

## Hilfsfunktionen
##*****************
 # R-Paket fuer "Kernel Regression Smoothing with Local Plug-in Bandwidth"
 #*************************************************************************
library( lokern)

 # Grafikfunktion
 #****************
TitleAndTopAxis <- function( at, digits = 4, ...) {
 op <- par( ...);      at0 = par( "usr")[ 1];   lh <- 1.1 * par( "cex.axis")
 B <- signif( Biases, digits);   V <- signif( Variances, digits)
 M <- signif( MSEs, digits)
 axis( 3, at = c( at0, at), tick = FALSE, hadj = 0, line = 0,
       labels = c( expression( paste( widehat(plain(MSE)), ":", sep = "")),
                   bquote( .(M[ 1])), bquote( .(M[ 2])), bquote( .(M[ 3])) ) )
 axis( 3, at = c( at0, at), tick = FALSE, hadj = 0, line = 1 * lh,
       labels = c( expression( paste( widehat(plain(Var)), ":", sep = "")),
                   bquote( .(V[ 1])), bquote( .(V[ 2])), bquote( .(V[ 3])) ) )
 axis( 3, at = c( at0, at), tick = FALSE, hadj = 0, line = 2 * lh,
       labels = c( expression( paste( widehat(plain(Bias)), ":", sep = "")),
                   bquote( .(B[ 1])), bquote( .(B[ 2])), bquote( .(B[ 3])) ) )
 axis( 3, at = at, tick = FALSE, hadj = 0, line = 2.9 * lh,
       labels = c( expression( hat(m)[n](x[0])), expression( m[n](x[0])),
                   expression( plain(lokerns)(x[0])) ) )

 title( bquote( list( h[adj] == {n^{-1/5} == .(signif( h, digits))},
                      ~~h[ideal] == .(signif( h_opt, digits)),
                      ~~bar(h[lokerns]) == .(signif( lox0.meanh, digits)),
                      ~~~x[0] == .(x0) )),
        line = 4.2 * lh, outer = TRUE)
 title( bquote( list( "AMSE-optimal,   K = Gauss kernel", ~~S == .(MCsamplesize),
                      ~~n == .(n), ~~sigma(paste(Y,"|",X==x)) == .(sd.eps) )),
        line = 5.4 * lh, outer = TRUE)

 invisible( par( op))
 }


 # Potenzielle Kern- und damit assoziierte Funktionen
 #****************************************************
Kern <- function( x) dnorm( x)   # Gausskern 

Integral_Kquadrat  <- integrate( function( u) { k <- Kern( u); k*k },
                                 lower = -Inf, upper = Inf)$value
                      # = 1 / (2 * sqrt( pi)), falls Kern Gausskern
Integral_xquadratK <- integrate( function( u) u*u * Kern( u),
                                 lower = -Inf, upper = Inf)$value
                      # = 1, falls Kern Gausskern 


 # Potenzielle Regressionsfunktionen
 #***********************************

#m.poly3 <- function( x, a0 = 1, a1 = -3, a2 = 1, a3 = 1)
# a0 + a1 * x + a2 * x*x + a3 * x*x*x   # m'(x)  = a1 + 2 a2 x + 3 a3 x^2
#                                       # m''(x) = 2 a2 + 6 a3 x 
#m.poly4 <- function( x, a0 = 1, a1 = -3, a2 = 1, a3 = 1, a4 = -1)
# a0 + a1 * x + a2 * x*x + a3 * x*x*x + a4 * x*x*x*x
#   # m'(x)  = a1 + 2 a2 x + 3 a3 x^2 + 4 a4 x^3
#   # m''(x) = 2 a2 + 6 a3 x + 12 a4 x^2
# # Effizienter: a0 + (a1 + (a2 + (a3 + a4 * x) * x) * x) * x

 # Gewuenscht ist ein Polynom 4. Grades m mit einem
 # Tief-/Hochpunkt, wo also m'(x) = 0 und m''(x) <> 0, einem
 # Wendepunkt,      wo also m''(x) = 0 und m'(x) <> 0 sowie einem
 # Sattelpunkt,     wo also m'(x) = 0 = m''(x).
 # Memo: Fuer f(x) = a * (x - x1) * (x - x2)^3 + b ist der
 # TP/HP an (3 * x1 + x2) / 4, der WP an (x1 + x2) / 2 und der SP an x2.
m.poly4 <- function( x, x1 = 0, x2 = 8, a = 0.01, b = 0)
 a * (x - x1) * (x - x2)^3 + b


#m.sin <- function( x, a0 = 0, a1 = 1, a2 = 1)   # m'(x)  =  a1 a2 cos( a2 x)
# a0 + a1 * sin( a2 * x)                         # m''(x) = -a1 a2^2 sin( a2 x)


 # Nadaraya-Watson-Schaetzer ("NW-Schaetzer", m_n im Paper)
 #**********************************************************
NadWat <- function( x, X, Y, h, K = Kern) {
  # x:    Vektor (x_1, ..., x_r) oder Skalar
  # X, Y: Vektoren (X_1, ..., X_n), (Y_1, ..., Y_n)
  # h:    Skalar
  # K:    Fkt. mit vektorieller Ein- & Ausgabe
 M <- K( outer( x/h, X/h, "-"))   # r x n
 drop( (M %*% Y) / rowSums( M))   # r x 1
 }


 # Gewichte
 #**********
Wn <- function( sigma, h, xXh, thetaXh, K = Kern) {
  # sigma:   Vektor (sigma_1, ..., sigma_s) oder Skalar
  # h:       Skalar
  # xXh:     Vektor ((x - X_1)/h, ..., (x - X_n)/h)
  # thetaXh: Vektor ((\theta - X_1)/h, ..., (\theta - X_n)/h)
  # K:       Fkt. mit vektorieller Ein- & Ausgabe

 A <- outer( sigma/h, xXh)       # s x n
 A <- outer( A, thetaXh, "+")    # s x n x n
 A <- rowSums( K( A), dims = 2)  # s x n
 A / rowSums( A)   # s x n: (W_{ni}( x, \sigma_r))_{1<=r<=s, 1<=i<=n}
 }


 # Zwei Versionen des Bias-Schaetzers: BiasA() mit m_n, BiasB() mit \hat m_n
 #****************************************************************************
BiasA <- function( sigma, h, xXh, thetaXh, K = Kern, mmDiff) {
  # sigma, X, h, xXh, thetaXh, K: Siehe Wn()!
  # mmDiff: Vektor (m_n(X_1) - m_n(x), ..., m_n(X_n) - m_n(x))

 drop( Wn( sigma, h, xXh, thetaXh, K) %*% mmDiff)  # s x 1
 }


BiasB <- function( sigma, Y, h, xXh, XXh, thetaXh, K = Kern) {
  # sigma, h, xXh, thetaXh, K: Siehe Wn()!
  # Y:   Vektor (Y_1, ..., Y_n)
  # XXh: Matrix ((X_j - X_i)/h)_{1<=j<=n, 1<=i<=n}

 wx <- Wn( sigma, h, xXh, thetaXh, K)   # s x n

 wX <- apply( XXh, 2, function( XXih) Wn( sigma, h, XXih, thetaXh, K))
       # (1<=r<=s)*(1<=j<=n) x (1<=i<=n)
 wX <- array( wX, dim = c( length( sigma), dim( XXh)))
       # (1<=r<=s) x (1<=j<=n) x (1<=i<=n)
 wX <- aperm( wX - as.vector( wx), c( 2, 1, 3))
       # (1<=j<=n) x (1<=r<=s) x (1<=i<=n)
 rowSums( wx * colSums( wX * Y))   # s x 1
               # (1<=r<=s) x (1<=i<=n)
 }


 # Zwei Versionen des Varianzschaetzer: VarA() mit m_n, VarB() mit \hat m_n
 #***************************************************************************
VarA <- function( sigma, h, xXh, thetaXh, K = Kern, YmDiff2) {
  # sigma, h, xXh, thetaXh, K: Siehe Wn()!
  # YmDiff2: Vektor ( (Y_1 - m_n(x))^2, ..., (Y_n - m_n(x))^2 )

 wx <- Wn( sigma, h, xXh, thetaXh, K)   # s x n
 drop( (wx * wx) %*% YmDiff2)  # s x 1
 }


VarB <- function( sigma, Y, h, xXh, XXh, thetaXh, K = Kern) {
  # sigma, h, xXh, thetaXh, K: Siehe Wn()!
  # Y:   Vektor (Y_1, ..., Y_n)
  # XXh: Matrix ((X_j - X_i)/h)_{1<=j<=n, 1<=i<=n}

 wx <- Wn( sigma, h, xXh, thetaXh, K)   # s x n

 wX <- apply( XXh, 2, function( XXih) Wn( sigma, h, XXih, thetaXh, K))
       # (1<=r<=s)*(1<=j<=n) x (1<=i<=n)
 wX <- aperm( array( wX, dim = c( length( sigma), dim( XXh))), c( 2, 1, 3))
       # (1<=j<=n) x (1<=r<=s) x (1<=i<=n)
 Y_twXY <- Y - t( colSums( wX * Y))   # (1<=i<=n) x (1<=r<=s)
 
 rowSums( (wx * wx) * t( Y_twXY * Y_twXY))  # s x 1
 }


 # Zwei Versionen des geschaetzten "Asymptotic Mean Squared Errors"
 # (\widehat AMSE): AMSE_A() mit m_n, AMSE_B() mit \hat m_n
 #*****************************************************************
AMSE_A <- function( sigma, h, xXh, thetaXh, K = Kern, mmDiff, YmDiff2) {
  # sigma, h, xXh, thetaXh, K: Siehe Wn()!
  # mmDiff:  Vektor (m_n(X_1) - m_n(x), ..., m_n(X_n) - m_n(x))
  # YmDiff2: Vektor ( (Y_1 - m_n(x))^2, ..., (Y_n - m_n(x))^2 )
 W <- Wn( sigma, h, xXh, thetaXh, K)  # s x n
 Bn  <- W %*% mmDiff          # = BiasA( sigma, h, xXh, thetaXh, K, mmDiff)
 Vn2 <- (W * W) %*% YmDiff2   # = VarA( sigma, h, xXh, thetaXh, K, YmDiff2)
 drop( Vn2 + Bn * Bn)         # s x 1
 }


AMSE_B <- function( sigma, Y, h, xXh, XXh, thetaXh, K = Kern) {
  # sigma, h, xXh, thetaXh, K: Siehe Wn()!
  # Y:   Vektor (Y_1, ..., Y_n)
  # XXh: Matrix ((X_j - X_i)/h)_{1<=j<=n, 1<=i<=n}

 wx <- Wn( sigma, h, xXh, thetaXh, K)   # s x n

 wX <- apply( XXh, 2, function( XXih) Wn( sigma, h, XXih, thetaXh, K))
       # (1<=r<=s)*(1<=j<=n) x (1<=i<=n)
 wX <- array( wX, dim = c( length( sigma), dim( XXh)))
       # (1<=r<=s) x (1<=j<=n) x (1<=i<=n)

 wX_wx <- aperm( wX - as.vector( wx), c( 2, 1, 3))
          # (1<=j<=n) x (1<=r<=s) x (1<=i<=n)
 Bn <- rowSums( wx * colSums( wX_wx * Y))
       # = BiasB( sigma, Y, h, xXh, XXh, thetaXh, K),   s x 1

 Y_twXY <- Y - t( colSums( aperm( wX, c( 2, 1, 3)) * Y)) # (1<=i<=n) x (1<=r<=s)
 Vn2 <- rowSums( (wx * wx) * t( Y_twXY * Y_twXY))
        # = VarB( sigma, Y, h, xXh, XXh, thetaXh, K),   s x 1

 drop( Vn2 + Bn * Bn)         # s x 1
 }



## Die Simulations-Szenarien
##***************************

##Version <- "A"   # Varianz- und Bias-Schaetzung mit NW-Schaetzer m_n fuer m
Version <- "B"   # Varianz- und Bias-Schaetzung mit \hat m_n fuer m

 # Simulationsparameter
x00 <- c( 2, 4, 5, 8)   # Punkte, an denen die AMSE-optim. "kernel adjusted
                        # nonparametric regressions" von m stattfinden sollen.
                        # Fuer m.poly4()'s Voreinstellung ist TP an 2, WP an 4
                        # und SP an 8; ein "beliebiger" Punkt liegt z. B. bei 5.
nn  <- c( 25, 50, 100)  # Stichprobenumfaenge.
sdd <- c( 0.5, 1, 2)    # Std.-abwn. der Fehler \epsilon_i: \sqrt(Var(Y|X=x)).

Resultate <- array( NA, dim = c( length( x00), length( nn), length( sdd), 3, 3),
                    dimnames = list( x0 = paste( "x0", x00, sep = "="),
                                     n  = paste( "n",  nn,  sep = "="),
                                     sd = paste( "sd", sdd, sep = "="),
                                     Schaetzer = c( "Bias", "Var", "MSE"),
                                     Methode   = c( "adjusted", "idealistic",
                                                    "classic")) )

mm <- m.poly4  # Regressionsfunktion und ihre erste sowie zweite Ableitung
mm_prime  <- function() {};   body( mm_prime)  <- D( body( mm), "x")
mm_second <- function() {};   body( mm_second) <- D( body( mm_prime), "x")
formals( mm_second) <- formals( mm_prime) <- formals( mm)

 # Dichte der Normalverteilung (der X_i) und ihre erste Ableitung
f <- dnorm;     f_prime <- function( x) -x * dnorm( x)


for( x0 in x00) {   # Punkt x_0

 design.mean <- x0 + 0.5   # "Zentrum" und "Skala" des Designs = Erwartungswert
 design.sd <- 1            # bzw. Standardabweichung der Verteilung der X_i.
 xx <- design.mean + seq( -3.3 * design.sd, 3.3 * design.sd, length = 101)
    # x-Gitter, worauf m und der NW-Schaetzer m_n ausgewertet & gezeichnet werden.

 sig.range <- c( 0.01, 1.5 * diff( range( xx))) # Zu betrachtender Wertebereich
 sig.grid <- seq( sig.range[ 1],              # samt Gitter fuer den Skalenparame-
                  sig.range[ 2], length = 501)# ter \sigma des adaptiv. Verfahrens.

 for( n in nn) { # Stichprobenumfang
  for( sd.eps in sdd) { # Standardabweichung der \epsilon_i

   set.seed( 20110426)   # Startwert des Pseudo-Zufallszahlengenerators

    # Berechnung der fuer x_0 optimalen, "idealistischen" Fensterbreite h_opt
    # zur Berechnung des NW-Schaetzers an x_0
   numerator   <- sd.eps^2  * Integral_Kquadrat * f( x0)
   denominator <- n * Integral_xquadratK^2 *
                  (f( x0) * mm_second( x0) + 2 * f_prime( x0) * mm_prime( x0))^2
   h_opt <- (numerator / denominator)^(1/5)

    # Wahl von h fuer das adaptive Verfahren
   h <- n^(-1/5)

    # Vorbereitung des Grafik-Ausgabe
   now <- format( Sys.time(), "%Y%m%d_%H%M%S")
   pdf( paste( "Plots/mBVA_x", x0, "_n", formatC( n, width = 3, flag = 0), "_s",
               formatC( 10*sd.eps, width = 2, flag = 0), "_D", now, ".pdf",
               sep = ""),
        paper = "a4", width = 0, height = 0)
   op <- par( mfrow = c( 3,1), mar = c( 3,3,3,0.1)+0.1, mgp = c( 1.5,0.5,0),
              tcl = -0.3, cex.main = 2)

    ## Ein paar beispielhafte Simulationslaeufe mit Grafikgenerierung
    ##****************************************************************
   for( bsp in 1:4) {
     # Datengenerierung Y_i = m( X_i) + \epsilon_i, i = 1, ..., n
    X <- rnorm( n, mean = design.mean, sd = design.sd)
    Y <- mm( X) + rnorm( length( X), sd = sd.eps)

     # NW-Schaetzer auf einem x-Gitter, um ihn zur Anschauung zu zeichnen
    mnxx <- NadWat( xx, X, Y, h_opt)

     # Kern-Schaetzer mit lokal-adaptiver "plug-in"-Bandbreite (aus dem R-Paket
     # "lokern"), auf selbem x-Gitter, um auch ihn zur Anschauung zu zeichnen
    loxx <- lokerns( x = X, y = Y, x.out = xx)$est

     # Vorratsberechnung von Groessen, die nicht von \sigma abhaengen, aber
     # i. F. wiederholt gebraucht werden; Algor. wird so effizienter.
    x0Xh <- (x0 - X) / h;   thetaXh <- (mean( X) - X) / h

    if( Version == "A") {   # Var.- & Bias-Schaetzung mit m_n

      # NW-Schaetzer m_n(x0) und m_n(X_1), ..., m_n(X_n) mit Fensterbreite
      # h = n^(-1/5) sowie daraus gebildete Groessen zur Schaetzung von
      # Var_x0(\sigma) und Bias_x0(\sigma)
     mnx0 <- NadWat( x0, X, Y, h)   # m_n(x0)
     mnX  <- NadWat( X, X, Y, h)    # m_n(X_i) fuer i = 1, ..., n
     mnX_mnx0 <- mnX - mnx0         # m_n(X_i) - m_n(x0) fuer i = 1, ..., n
     Y_mnX2 <- (Y - mnX)^2          # (Y_i - m_n(X_i))^2 fuer i = 1, ..., n

      # Schaetzer von Var_x0(\sigma) und Bias_x0(\sigma) auf dem \sigma-Gitter,
      # um sie zur Anschauung zu zeichnen
     Bn <- BiasA( sigma = sig.grid, h = h, xXh = x0Xh, thetaXh = thetaXh,
                  mmDiff = mnX_mnx0)
     Vn2 <- VarA( sigma = sig.grid, h = h, xXh = x0Xh, thetaXh = thetaXh,
                  YmDiff2 = Y_mnX2)

#      # Versuch 1: Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma):
#     opt <- optimize( f = AMSE_A, interval = sig.range, h = h, xXh = x0Xh,
#                      thetaXh = thetaXh, mmDiff = mnX_mnx0, YmDiff2 = Y_mnX2)
#     sig_opt_x0 <- opt$minimum;   amse_min <- opt$objective
#
#      # Versuch 2: Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma):
#     opt <- optimize( f = AMSE_A, interval = sig.range/c(1,3), h = h, xXh = x0Xh,
#                      thetaXh = thetaXh, mmDiff = mnX_mnx0, YmDiff2 = Y_mnX2)
#     sig_opt_x0_2 <- opt$minimum;   amse_min_2 <- opt$objective

     } else if( Version == "B") { # Var.- & Bias-Schaetzung mit \hat m_n

      # Weitere Vorratsgroesse, die nicht von \sigma abhaengt, aber i. F.
      # wiederholt gebraucht wird
     XjXih <- outer( X/h, X/h, "-")

      # Schaetzer von Var_x0(\sigma) und Bias_x0(\sigma) auf dem \sigma-Gitter,
      # um sie zur Anschauung zu zeichnen
     Bn <- BiasB( sigma = sig.grid, Y = Y, h = h, xXh = x0Xh, XXh = XjXih,
                  thetaXh = thetaXh)
     Vn2 <- VarB( sigma = sig.grid, Y = Y, h = h, xXh = x0Xh, XXh = XjXih,
                  thetaXh = thetaXh)

#      # Versuch 1: Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma):
#     opt <- optimize( f = AMSE_B, interval = sig.range, Y = Y, h = h,
#                      xXh = x0Xh, XXh = XjXih, thetaXh = thetaXh)
#     sig_opt_x0 <- opt$minimum;   amse_min <- opt$objective
#
#      # Versuch 2: Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma):
#     opt <- optimize( f = AMSE_B, interval = sig.range/c(1,3), Y = Y, h = h,
#                      xXh = x0Xh, XXh = XjXih, thetaXh = thetaXh)
#     sig_opt_x0_2 <- opt$minimum;   amse_min_2 <- opt$objective

    } else stop( "Version '", Version, "' gibt es nicht!")

#     # Versuch 2 "besser" als Versuch 1?
#    if( amse_min_2 < amse_min) {
#     sig_opt_x0 <- sig_opt_x0_2; amse_min <- amse_min_2
#     }

     # "Komposition" von AMSE_x0(\sigma) aus den Schaetzern von Var_x0(\sigma)
     # und Bias_x0(\sigma) auf dem \sigma-Gitter
    An <- Vn2 + Bn * Bn

     # Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma)
     # durch "grid search":
    min_idx <- which.min( An)
    sig_opt_x0 <- sig.grid[ min_idx];   amse_min <- An[ min_idx]

     # Der AMSE-optimale "kernel adjusted nonparametric regression"-Schaetzwert
     # an x_0, also \hat m_n( x_0)
    W <- Wn( sigma = sig_opt_x0, h = h, xXh = x0Xh, thetaXh = thetaXh)
    mnhatx0 <- W %*% Y

     # Grafische Darstellung fuer den aktuellen Datensatz
    plot( Y ~ X, xlim = range( X, xx),
          ylim = range( Y, mnxx, mnhatx0, loxx, na.rm = TRUE),# mm(xx)
          main = bquote( n == .(n)), xlab = "x", ylab = "y") # Das Streudiagramm
    lines( xx, mm( xx), lty = 2)    # Die "wahre" Regressionsfunktion m
    lines( xx, mnxx)                # Der NW-Schaetzer m_n von m
    lines( xx, loxx, col = "blue")  # Der lokal-adaptive Schaetzer von m
    points( x0, mnhatx0, col = "red", pch = 8, cex = 2)  # (x_0, \hat m_n(x_0))
    legend( "bottomleft", lty = c( 2, 1, 1), col = c( "black", "black", "blue"),
            bty = "n", cex = 1.5,
            legend = c( bquote( paste( "m, ", sigma(paste( Y,"|", X==x))
                                              == .(sd.eps)) ),
                        expression( m[n]), "lokern()"))
    legend( "top", pch = 8, col = "red", bty = "n", cex = 1.5,
            legend = bquote( paste( hat(m)[n](x[0]), ",   ", x[0] == .(x0))))

    matplot( sig.grid, cbind( Bn*Bn, Vn2), type = "l", lty = 1:2,
             col = c( "black", "red"), xlab = expression( sigma), ylab = "")
    legend( "top", lty = 1:2, col = c( "black", "red"), bty = "n", cex = 1.5,
            legend = c( expression( paste( widehat(plain(Bias))[x[0]]^2, (sigma))),
                        expression( widehat(plain(Var))[x[0]](sigma))))

    plot( sig.grid, An, type = "l", xlab = expression( sigma), ylab = "",
          main = expression( widehat(plain(AMSE))[x[0]](sigma)))
    points( sig_opt_x0, amse_min, pch = 13, col = "red", cex = 2)
    legend( "top", pch = 13, col = "red", bty = "n", cex = 1.5,
            legend = substitute( group( "(", list( plain( Minimizer),
                                                   plain( Minimum)), ")")
                                 == group( "(", list( x, y), ")"),
                                 list( x = signif( amse_min, 4),
                                       y = signif( sig_opt_x0, 4))))
    } # bsp-loop
   par( op)
   dev.off() ## Ende der beispielhaften Simulationslaeufe


    ## Die eigentliche Simulation   (Kommentierung siehe oben!)
    ##****************************
   MCsamplesize <- 5000

   now <- format( Sys.time(), "%Y%m%d_%H%M%S")
   set.seed( 20110426) # Derselbe Startwert des Pseudo-RNGs wie oben.
   Sim <- replicate( MCsamplesize, {   # Schleife ueber Monte-Carlo-Samples
    X <- rnorm( n, mean = design.mean, sd = design.sd)
    Y <- mm( X) + rnorm( length( X), sd = sd.eps)
 
if( n > 25 || sd.eps == 1) {
  c( mnhatx0 = NA, mnx0_hopt = NA, lox0 = NA, lox0.h = NA)
 } else {
    mnx0_hopt <- NadWat( x0, X, Y, h_opt)
    lokerns.x0 <- lokerns( x = X, y = Y, x.out = x0)
    lox0 <- lokerns.x0$est;     lox0.h <- lokerns.x0$bandwidth

    x0Xh <- (x0 - X) / h;   thetaXh <- (mean( X) - X) / h
   
    if( Version == "A") {
     mnx0 <- NadWat( x0, X, Y, h);     mnX <- NadWat( X, X, Y, h)
     mnX_mnx0 <- mnX - mnx0;           Y_mnX2 <- (Y - mnX)^2
#      # Minimierungsversuch 1:
#     opt <- optimize( f = AMSE_A, interval = sig.range, h = h, xXh = x0Xh,
#                      thetaXh = thetaXh, mmDiff = mnX_mnx0, YmDiff2 = Y_mnX2)
#     sig_opt_x0 <- opt$minimum;   amse_min <- opt$objective
#      # Minimierungsversuch 2:
#     opt <- optimize( f = AMSE_A, interval = sig.range/c(1,3), h = h, xXh = x0Xh,
#                      thetaXh = thetaXh, mmDiff = mnX_mnx0, YmDiff2 = Y_mnX2)
#     sig_opt_x0_2 <- opt$minimum;   amse_min_2 <- opt$objective

     An <- AMSE_A( sig.grid, h = h, xXh = x0Xh, thetaXh = thetaXh,
                   mmDiff = mnX_mnx0, YmDiff2 = Y_mnX2)

    } else if( Version == "B") {
     XjXih <- outer( X/h, X/h, "-")
#      # Minimierungsversuch 1:
#     opt <- optimize( f = AMSE_B, interval = sig.range, Y = Y, h = h,
#                      xXh = x0Xh, XXh = XjXih, thetaXh = thetaXh)
#     sig_opt_x0 <- opt$minimum;   amse_min <- opt$objective
#      # Minimierungsversuch 2:
#     opt <- optimize( f = AMSE_B, interval = sig.range/c(1, 3), Y = Y, h = h,
#                      xXh = x0Xh, XXh = XjXih, thetaXh = thetaXh)
#     sig_opt_x0_2 <- opt$minimum;   amse_min_2 <- opt$objective

     An <- AMSE_B( sig.grid, Y = Y, h = h, xXh = x0Xh, XXh = XjXih,
                   thetaXh = thetaXh)

    } else stop( "Version '", Version, "' gibt es nicht!")

#     # Minimierungsversuch 2 "besser" als 1?
#    if( amse_min_2 < amse_min) sig_opt_x0 <- sig_opt_x0_2

     # Minimizer und Minimum des Schaetzers von AMSE_x0(\sigma)
     # durch "grid search":
    sig_opt_x0 <- sig.grid[ which.min( An)]

    W <- Wn( sigma = sig_opt_x0, h = h, xXh = x0Xh, thetaXh = thetaXh)
    mnhatx0 <- W %*% Y

    c( mnhatx0 = mnhatx0, mnx0_hopt = mnx0_hopt, lox0 = lox0, lox0.h = lox0.h)
 }
    } )

   lox0.meanh <- mean( Sim[ "lox0.h",]);     Sim <- Sim[ -4,]
   mx0 <- mm( x0);     Means <- rowMeans( Sim);     Biases <- Means - mx0
   Variances <- apply( Sim, 1, var);     MSEs <- Biases * Biases + Variances

   Resultate[ paste( "x0", x0, sep = "="), paste( "n", n, sep = "="),
              sd = paste( "sd", sd.eps, sep = "="),
              c( "Bias", "Var", "MSE"), ] <- rbind( Biases, Variances, MSEs)

if( !(n > 25 || sd.eps == 1)) {
    ## Grafische Darstellung der Resultate (Boxplots)
    ##************************************************
   Bim1 <- as.list( as.data.frame( t( Sim)))
   Bim2 <- list( "'Idealistic' - Adjusted" = Sim[ "mnx0_hopt",] - Sim[ "mnhatx0",],
                 "'Classic' - Adjusted"    = Sim[ "lox0",]      - Sim[ "mnhatx0",])

   pch <- c( 1, 4, 5);    cols <- c( "red", "orange", "green");     cex <- 2
   pdf( paste( "Plots/AMSEopt_x", x0, "_n", formatC( n, width = 3, flag = 0),
               "_s", formatC( 10*sd.eps, width = 2, flag = 0), "_D", now,
               ".pdf", sep = ""),
        paper = "a4", width = 0, height = 0)
   par( pl <- list( oma = c( 0,0,8,0), mar = c( 3.5, 3.3, 1, 3)+.1,
                    mgp = c( 2, 0.5, 0), tcl = -0.3, lab = c( 5, 10, 0)))
    boxplot( Bim1, names = c( "Adjusted ...", "'Idealistic' ...", "'Classic' ..."),
             ylab = expression( hat(m)[n](x[0])~~~~"bzw."~~~~m[n](x[0])
                                               ~~~~"bzw."~~~~lokerns(x[0])),
             cex.axis = 1.2, cex.lab = 1.2,
             xlab = expression( "kernel regression estimate values at"~x[0]))

    abline( h = mx0, col = "blue")
    axis( 4, at = mx0, tick = FALSE, labels = expression( m(x[0])), line = 0,
          las = 2, cex.axis = 1.2, col.axis = "blue")

    points( Means, col = cols, pch = pch, cex = cex, lwd = 3)
    legend( "top", pch = pch, col = cols, pt.lwd = 3, bty = "n", cex = 1.4,
            legend = c( expression( bar(hat(m)[n](x[0]))),
                        expression( bar(m[n](x[0]))),
                        expression( bar(lokerns(x[0]))) ), y.intersp = 1.7)
    TitleAndTopAxis( at = 1:3 - 0.18, cex.axis = 1.2)


    boxplot( Bim2, cex.axis = 1.1, cex.lab = 1.1,
             ylab = expression( m[n](x[0]) - hat(m)[n](x[0])~~~~~~"bzw."~~~~~~
                                lokerns(x[0]) - hat(m)[n](x[0])),
             xlab = bquote( "Difference between kernel regression estimates at "~
                            x[0] == .(x0)~"with"~
                            sigma(paste(Y,"|",X==x)) == .(sd.eps)))
    abline( h = 0, col = "blue")
    points( c( Means[ "mnx0_hopt"] - Means[ "mnhatx0"],
               Means[ "lox0"] - Means[ "mnhatx0"]),
            col = c( "darkgreen", "brown"), pch = c( 8, 9), cex = cex, lwd = 2)
    legend( "topleft", pch = c( 8, 9, -1), col = c( "darkgreen", "brown", NA),
            pt.lwd = 2, bty = "n", cex = 1.4,
            legend = c( bquote( bar(m[n](x[0]) - hat(m)[n](x[0])) == 
                                .(round( Means[ "mnx0_hopt"] - Means[ "mnhatx0"],
                                         5))),
                        bquote( bar(lokerns(x[0]) - hat(m)[n](x[0]))== 
                                .(round( Means[ "lox0"] - Means[ "mnhatx0"], 5))),
                        expression() ), y.intersp = 1.7)
    TitleAndTopAxis( at = c( 0.833, 1.5, 2.165) - 0.12, cex.axis = 1.2)
   dev.off()
 }
   } # sd.eps-loop
  } # n-loop
 } # x0-loop

# rm( list = ls()); dev.off()


today <- format( Sys.time(), "%Y%m%d")
sink( paste( "Resultate_", today, ".txt", sep = ""))
 print( ftable( Resultate))
sink()

#round( ftable( Resultate[ c( "x0=4", "x0=8", "x0=2"),
#                          c( "n=25", "n=50", "n=100"),
#                          c( "sd=0.5", "sd=2"), , ],
#               row.vars = c( "x0", "sd", "Schaetzer") ), 4)

library( Hmisc)
punkte <- c( "minimum" = "x0=2", "point of inflection" = "x0=4",
             "'arbitrary' point" = "x0=5", "saddle point" = "x0=8")
Tab <- lapply( punkte, function( v)
                        ftable( Resultate[ v, c( "n=50", "n=25", "n=100"),
                                           c( "sd=0.5", "sd=2"), , ],
                                row.vars = c( "sd", "Schaetzer")) )
names( Tab) <- paste( "Performance in", names( punkte), 
                      paste( "$(", gsub( "x0", "x_0", punkte, fixed = TRUE),
                             ")$", sep = ""),
                      "based on", MCsamplesize, "Monte Carlo samples" )

CV <- attr( Tab[[1]], "col.vars")
CV$n <- paste( "$", CV$n, "$", sep = "")
CV$Methode2 <- c( "$\\hat m_n$", "$m_n$", "\\texttt{lokerns}")

RV <- attr( Tab[[1]], "row.vars")
RV$sd <- paste( "$", gsub( "sd", "\\sigma", RV$sd, fixed = TRUE), "$", sep = "")
RV$Schaetzer <- paste( "$\\widehat{\\text{", RV$Schaetzer, "}}$", sep = "")

w <- latex( object = Tab, file = paste( "Tab_", today, ".tex", sep = ""),
            title = "", where = "!htbp", dec = 3, first.hline.double = FALSE, 
            cgroup = CV$n, colheads = rep( CV$Methode2, length( CV$n)),
            extracolheads = rep( CV$Methode, length( CV$n)),
            rgroup = RV$sd, rowname = rep( RV$Schaetzer, length( RV$sd))
            )
