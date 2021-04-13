# =============================================================================
#'
#' @title Unfozen water in soil
#'
#' @description Computes the proportion of unfozen water and ice in soil as 
#'              well as the temperature derivative of liquid water content. 
#'
#' @details This function applies invariant freezing parameters with depth and
#'          assumes saturated conditions without movement of water.
#' 
#'
#' @param mat Data frame with relevant ground properties.  
#' @param unfrozen.type Character string indicating the type of function to be
#'                      used; the standard is "INTERVAL". The full options are: 
#"
#'                      "DALLAMICO": Dall’Amico, M., Endrizzi, S., Gruber, S., 
#'                      & Rigon, R. (2011). A robust and energy-conserving 
#'                      model of freezing variably-saturated soil. The 
#'                      Cryosphere, 5(2), 469–484. doi:10.5194/tc-5-469-2011
#'                      Typical values of van Genuchten parameters are found 
#'                      in Table 2 of Gubler, S., Endrizzi, S., Gruber, S., 
#'                      & Purves, R. S. (2013). Sensitivities and uncertainties 
#'                      of modeled ground temperatures in mountain environments. 
#'                      Geoscientific Model Development, 6(4), 1319–1336. 
#'                      doi:10.5194/gmd-6-1319-2013
#'
#'                      "MOTTAGHY": Mottaghy, D., & Rath, V. (2006). Latent heat 
#'                      effects in subsurface heat transport modelling and their
#'                      impact on palaeotemperature reconstructions. Geophysical 
#'                      Journal International, 164(1), 236-245.
#'
#'                      "INTERVAL": Phase change takes place in an interval of
#'                      specified with (unfrozen.par) below 0C and at a 
#'                      constant rate.  	
#'
#' @param unfrozen.par Parameter set for the chosen unfrozen water function.
#'
#'                     "DALLAMICO": unfrozen.par[1]: van Genuchten alpha [mm-1], 
#'                     unfrozen.par[2]: van Genuchten n [-], unfrozen.par[3]: 
#'                     residual water content [m3/m3]. The saturated water 
#'                     content is given by the input (mat$wat) and can thus 
#'                     vary with depth. 
#'
#'                     "MOTTAGHY": unfrozen.par[1]: width of freezing 
#'                     interval [K], unfrozen.par[2]: omega.  
#'
#'                     "INTERVAL" unfrozen.par[1]: width of freezing 
#'                     interval [K].  
#' 
#' @return Returns the input data frame with three columns updated: liq [m3/m3]
#'         (liquid water content, relative to soil volume), ice [m3/m3] (ice
#'         content, relative to soil volume), and dice [m3/C] (change in ice
#'         content, relative to soil volume per degree Celsius). If these columns
#'         do not exist in the input data frame, they are created.   
#' 
#' @export
#' @examples
#' mat <- data.frame(Tj  = (-500:200)/100, 
#'                   wat = rep(0.5,701))
#'
#' #Mottahgy and Rath (2006) example
#' mat <- Unfrozen(mat, unfrozen.type = "MOTTAGHY", unfrozen.par = c(0, 0.5))
#' plot(mat$Tj, mat$dice, type="l", lty = 1, col = "black", ylim =c(0,2),
#'      xlab = "Temperature [C]", 
#'      ylab = "Liquid water content (solid), derivation (dashed)")
#' lines(mat$Tj, mat$liq, lty = 2, col = "black")
#'
#' #interval function example
#' mat <- Unfrozen(mat, unfrozen.type = "INTERVAL", unfrozen.par = 1)
#' lines(mat$Tj, mat$dice, lty = 1, col = "blue")
#' lines(mat$Tj, mat$liq,  lty = 2, col = "blue")
#'
#' #Dall'Amico et al. (2011) example, assuming saturated conditions.
#' mat <- Unfrozen(mat, unfrozen.type = "DALLAMICO", 
#'                 unfrozen.par = c(0.001, 1.4, 0.05))
#' lines(mat$Tj, mat$dice, lty = 1, col = "red")
#' lines(mat$Tj, mat$liq,  lty = 2, col = "red")
#'
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================
Unfrozen <- function(mat, unfrozen.type = "INTERVAL", unfrozen.par = 0.5) {
  
  unfrozen.ok <- FALSE
  
  #--- Interval function ------------------------------------------------------
  if (toupper(unfrozen.type) == "INTERVAL") {
  	#simple freezing interval
  	#derivative of fractional volumetric ice content
  	inrange  <- (mat$Tj < 0) * (mat$Tj > (-unfrozen.par))
  	mat$dice <- mat$wat / unfrozen.par * inrange  

    #fractional volumetric ice content
  	mat$ice <- mat$Tj/-unfrozen.par
  	mat$ice <- ifelse(mat$ice < 0, 0, mat$ice)
  	mat$ice <- ifelse(mat$ice > 1, 1, mat$ice)
  	mat$ice <- mat$wat * mat$ice
  	mat$liq <- mat$wat - mat$ice
  	unfrozen.ok <- TRUE
  }
  
  #--- Mottaghy and Rath (2006) -----------------------------------------------
  if (toupper(unfrozen.type) == "MOTTAGHY") {
  	#Mottaghy, D., & Rath, V. (2006). Latent heat effects in subsurface heat
  	#transport modelling and their impact on palaeotemperature reconstructions.
  	#Geophysical Journal International, 164(1), 236-245.
  	#Equations 3 and 8; unfrozen.par[1]: TL, unfrozen.par[2]: omega
  
    #fractional volumetric ice content, Eq. 3
    mat$liq <- exp(-((mat$Tj - unfrozen.par[1]) / unfrozen.par[2])^2)
    mat$liq <- ifelse(mat$Tj > unfrozen.par[1], 1, mat$liq)
  	mat$liq <- mat$wat * mat$liq
  	mat$ice <- mat$wat - mat$liq
  	
  	#derivative of fractional volumetric ice content, Eq 8
  	mat$dice <- exp(-((mat$Tj-unfrozen.par[1])/unfrozen.par[2])^2)  
  	mat$dice <- mat$dice * -2 * (mat$Tj-unfrozen.par[1]) / unfrozen.par[2]
  	mat$dice <- ifelse(mat$Tj > unfrozen.par[1], 0, mat$dice)

  	unfrozen.ok <- TRUE
  }

  #--- Dall'Amico et al. (2011) -----------------------------------------------
  if (toupper(unfrozen.type) == "DALLAMICO") {
  	#Dall’Amico, M., Endrizzi, S., Gruber, S., & Rigon, R. (2011). A robust and 
  	#energy-conserving model of freezing variably-saturated soil. The 
  	#Cryosphere, 5(2), 469–484. doi:10.5194/tc-5-469-2011
  	#We assume saturated conditions here.
  	alpha     <- rep(unfrozen.par[1], length(mat$Tj))
  	n         <- rep(unfrozen.par[2], length(mat$Tj))
  	theta.sat <- mat$wat
  	theta.res <- rep(unfrozen.par[3], length(mat$Tj))
  	
  	#catch zero water content problems
  	theta.res <- ifelse(theta.res > theta.sat, theta.sat, theta.res)
  	
  	#compute
  	res0 <- DallAmico(alpha, n, theta.sat, theta.res, mat$Tj+273,     273, mat$Tj*0) 
    res1 <- DallAmico(alpha, n, theta.sat, theta.res, mat$Tj+273.001, 273, mat$Tj*0) 

    #fractional volumetric ice/water content
  	mat$ice <- res0$theta.i
  	mat$liq <- res0$theta.w
  	mat$dice<- (res1$theta.w - res0$theta.w) * 1000 
  	unfrozen.ok <- TRUE
  }
  
  #feedback if function type was not recognized
  if (unfrozen.ok == FALSE) stop("unfrozen.type not recognized") 
  
  return(mat)
}
