# =============================================================================
#'
#' @title Analytic Temperature Field with Harmonic Variation
#'
#' @description Provides a template for bbuilding nicely documented R functions
#'
#' @details Feel free to leave out things that are not needed.
#'
#' @param z Depth [m] below the surface, positive downwards. Must be positive 
#'        number, can be a vector or a scalar.  
#' @param t Time [s] after start of period for which output is desired. 
#'        Can be a vector or a scalar
#' @param Tmean Mean temperature [C] at the surface. Must be a scalar.
#' @param As Amplitude [C] of surface temperature oscillation. Must be larger 
#'        than zero and a scalar.
#' @param P Period [s] of surface temperature oscillation. Must be larger 
#'        than zero and a scalar.
#' @param kappa Thermal difussivity [m2/s] of material.
#' 
#' @return Returns an array of temperatures [C] with one column for each z
#'         specified and one row for each t specified.   
#' 
#' @export
#' @examples
#' #temperature for only one depth and time
#' temp <- TempAnalyticHarmonic(0, 10, 1, 3600000, 24*365*3600, 5e-07)
#' 
#' #temperature for three depths and one instant in time
#' temp <- TempAnalyticHarmonic(0, 10, c(1,2,3), 3600000, 24*365*3600, 5e-07)
#' 
#' #temperature for three depths and many time steps, plot 
#' time <- seq(from = 0, to = 24*365*3600, length.out = 200)
#' temp <- TempAnalyticHarmonic(0, 10, c(1,2,3), time, 24*365*3600, 5e-07)
#' plot( temp[,1], type="l")
#' lines(temp[,2], col=2)
#' lines(temp[,3], col=3)
#' 
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================
TempAnalyticHarmonic <- function(Tmean, As, z, t, P, kappa) { 	

  #check input
  if (As <= 0)            stop("Parameter 'As' must be larger than zero.")
  if (min(z)  <  0)       stop("Parameter 'z' must be larger or equal to zero.")
  if (P <= 0)             stop("Parameter 'omega' must be larger than zero.")  
  if (length(Tmean) != 1) stop("Parameter 'Tmean' must be a scalar.")  
  if (length(As) != 1)    stop("Parameter 'As' must be a scalar.") 
  if (length(P) != 1) stop("Parameter 'omega' must be a scalar.") 
  if (length(kappa) != 1) stop("Parameter 'kappa' must be a scalar.") 
   
  #make arrays of z and t to have final dimensions of nz columns and nt rows
  nz <- length(z)
  nt <- length(t)
  z  <- t(array(rep(z, nt), dim = c(nz, nt)))
  t  <-   array(rep(t, nz), dim = c(nt, nz))
        
  #substitution and conversion
  omega <- 2 * pi / P
  sub   <- -z * sqrt(omega/(2 * kappa)) 
   
  #calculation 
  res <- Tmean + As * exp(sub) * sin(omega * t + sub) 
 
  return(res) 
 }
