# =============================================================================
#'
#' @title Compute total water content, ice content, and liquid water content of
#'        soil following Dall'Amico et al. (2011).
#'
#' @description Compute total water content, ice content, and liquid water 
#'              content of soil based on van Genuchten parameters, temperature 
#'              and soil matric potential. 
#'
#' @details Details are found in: Dall’Amico, M., Endrizzi, S., Gruber, S., & 
#'          Rigon, R. (2011). A robust and energy-conserving model of freezing 
#'          variably-saturated soil. The Cryosphere, 5(2), 469-484. 
#'          doi:10.5194/tc-5-469-2011
#'
#'          Typical values of van Genuchten parameters are found in Table 2 of
#'          Gubler, S., Endrizzi, S., Gruber, S., & Purves, R. S. (2013). 
#'          Sensitivities and uncertainties of modeled ground temperatures in 
#'          mountain environments. Geoscientific Model Development, 6(4), 
#'          1319–1336. doi:10.5194/gmd-6-1319-2013
#'
#' @param alpha van Genuchten alpha [mm-1]
#' @param n van Genuchten n [-]
#' @param theta.sat Saturated water content [m3/m3]
#' @param theta.res Residual water content [m3/m3]
#' @param Tg Ground temperature [K], can be a vector.
#' @param Tm Melting temperature of water at atmospheric pressure [K], 
#'           usually 273.15 K.
#' @param psi0 Soil matric potential [m]. Matric potential of 0  corresponds 
#'             to saturated conditions, smaller values imply unsaturated 
#'             conditions. 
#' 
#' @return Returns a data frame with three columns: theta.v (total water content, 
#'         relative to soil volume, as a function of psi0), theta.i (ice content,
#'         relative to soil volume), and theta.w (liquid water content, relative
#'         to soil volume. The data frame has as manu rows as Tg has elements.   
#' 
#' @export
#' @examples DallAmico(0.001, 1.6, 0.5, 0.02, 272, 273.15, 0)
#'
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================
DallAmico <- function(alpha, n, theta.sat, theta.res, Tg, Tm, psi0) {

  #=== PREPARATION ============================================================
  g  <- 9.81   #m/s2
  Lf <- 333700 #J/kg
  m  <- 1 - 1 / n
  alpha <- alpha * 1000 #convert from [mm-1] to [m-1]
  Tg <- as.numeric(Tg)
  Tm <- as.numeric(Tm)
  
  #catch temperature in C
  Tm <- ifelse(Tm < 150, Tm + 273.15, Tm)
  Tg <- ifelse(Tg < 150, Tg + 273.15, Tg)
  
  #psi0 above 0
  psi0 <- ifelse(psi0 > 0, 0, psi0)  
  
  #catch theta_sat < theta_res
  stopifnot(sum(theta.sat < theta.res) == 0)
  
  #=== DEPRESSED MELTING TEMPERATURE IN UNSATURATED CONDITIONS ================
  #Equation 17 of Dall'Amico et al 2011
  #partially frozen
   T.star <- Tm + (g * Tm) / Lf * psi0
   	
  #=== WATER PRESSURE UNDER FREEZING CONDITIONS ===============================
  #Equation 20 or Dall'Amico et al 2011
  #all liquid conditions
  Tg <- ifelse(Tg >= T.star, T.star, Tg)

  #partially frozen
  psiT <- psi0 + Lf / (g * T.star) * (Tg - T.star)
	
  #=== TOTAL WATER CONTENT ====================================================
  #Equation 22 or Dall'Amico et al 2011
  theta.v <- theta.res + (theta.sat - theta.res) * (1 + (-alpha * psi0)^n)^(-m)

  #=== LIQUID WATER CONTENT UNDER FREEZING CONDITIONS =========================
  #Equation 23 or Dall'Amico et al 2011
  theta.w <- theta.res + (theta.sat - theta.res) * 
             (1 + (-alpha * psiT)^n)^(-m)

  #=== ICE CONTENT UNDER FREEZING CONDITIONS ==================================
  # Equation 24 or Dall'Amico et al 2011
  theta.i <- theta.v - theta.w

  return(data.frame(theta.v = theta.v,
                    theta.i = theta.i,
                    theta.w = theta.w))
}