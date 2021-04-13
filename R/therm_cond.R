# =============================================================================
#'
#' @title Termal conductivity of soil/rock parametrized from its constituents
#'        and of snow parameterizad from density or constituents.
#'
#' @description Computes an estimated thermal condictivity for snow based on snow
#'              density or on liquid water and ice content in the snow pack. 
#'              Differing parameterizations are available. The result is an 
#'              effective thermal conductivity as various processes (conduction,
#'              advection, diffusion, radiative transfer) are parameterised. In the
#'              ground, differing mixing model are applied to the proprtions
#'              and thermo-physical properties of ground constituents. 
#'
#' @details Positive are to the right while left shifts are expressed as a
#'          negative number. All shifts are circular. Elements shifted off one
#'          end wrap around and are shifted onto the other end. This function
#'          mimicks the behaviour of SHIFT in IDL. 
#'
#' @param mat Data frame with relevant ground properties. Needed: mat$ice (ice
#'            content, volumetric fraction) and mat$liq (liquid water
#'            content, volumetric fraction).  
#'
#' @param type.sno Character string indicating the type of parameterization to be
#'             used for snow; the standard is CALONNE. Available are:
#'
#'             CALONNE: Calonne, N., Flin, F., Morin, S., Lesaffre, B., du 
#'             Roscoat, S. R., & Geindreau, C. (2011). Numerical and experimental
#'             investigations of the effective thermal conductivity of snow. 
#'             Geophysical Research Letters, 38(23), doi:10.1029/2011GL049234  
#'             "YEN": Yen, Y.-C. (1981). Review of thermal properties of 
#'              snow, ice and sea ice (34 pages). Hanover, NH, USA.
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             STURM: Sturm, M., J. Holmgren, M. König, and K. Morris (1997),
#'             The thermal conductivity of seasonal snow, Journal of 
#'             Glaciology, 43(143), 26–41.
#'
#'             JORDAN: Jordan, R. E., Andreas, E. L., & Makshtas, A. P. 
#'             (1999). Heat budget of snow-covered sea ice at North Pole 4. 
#'             Journal of Geophysical Research, 104(C4), 7785. 
#'             doi:10.1029/1999JC900011
#'
#'             WILLIAMS: Figure 4.11 in "The Frozen Earth: fundamentals of 
#'             geocryology" by P. Williams and M. Smith, 1989.
#' 
#' @param type.gnd Character string indicating the type of parameterization to be
#'             used for ground; the standard is "COSENZA". Available are:
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             GEOMETRIC: Geometric mean, intermediate mixed conductivity model, 
#'             approximation of randomly oriented consituent elements. 
#'
#'             ARITHMETIC: Arithmetic mean, high mixed conductivity model, 
#'             approximation of consituent elements layered parallel to 
#'             temperature gradient.
#'
#'             HARMONIC: Harmonic mean, low mixed conductivity model, 
#'             approximation of consituent elements layered normal to 
#'             temperature gradient.
#'
#' @param layers.sno Index to identify snow layers within 'mat'.
#'
#' @return Returns the input data frame with the column kj updated in [W m-1 K-1]
#'         If this column does not exist in the input it is created.   
#' 
#' @export
#' @examples
#' rho.i <-  917  # [kg m-3]
#' rho.l <- 1000  # [kg m-3]
#' mat   <- data.frame(ice = (0:600)/1000, 
#'                     liq = rep(0.0,601))
#' dens  <- mat$ice*rho.i + mat$liq*rho.l
#' 
#' mat <- TermCond(mat, type.sno = "CALONNE", layers.sno = 1:601)
#' 
#' plot(dens, mat$kj, type = "l", col = "black", main = 
#'      "Comparison of parameterization for snow thermal conductivity", 
#'      xlab = "Density [kg m-3]", ylab = "Thermal conductivity [W m-1 K-1]")
#' 
#' mat <- TermCond(mat , type.sno = "YEN", layers.sno = 1:601)
#' lines(dens, mat$kj, col = "red")
#' 
#' mat <- TermCond(mat , type.sno = "COSENZA", layers.sno = 1:601)
#' lines(dens, mat$kj, col = "blue")
#' 
#' mat <- TermCond(mat , type.sno = "STURM", layers.sno = 1:601)
#' lines(dens, mat$kj, col = "green")
#' 
#' mat <- TermCond(mat , type.sno = "WILLIAMS", layers.sno = 1:601)
#' lines(dens, mat$kj, col = "orange")
#' 
#' mat <- ermCond(mat , type.sno = "JORDAN", layers.sno = 1:601)
#' lines(dens, mat$kj, col = "darkorchid2")
#' 
#' legend("topleft", inset=.05, seg.len = 2, lty = 1, lwd =1,
#'   	   c("CALONNE","YEN","COSENZA", "STURM", "WILLIAMS", "JORDAN"), 
#'   	   col = c("black", "red", "blue", "green", "orange", "darkorchid2"))
#'
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================
TermCond <- function(mat, layers.sno = c(0), 
                     type.sno = "CALONNE", type.gnd = "COSENZA") {


  #process snow layers ========================================================
  mat.sno <- mat[ layers.sno,]
  if (nrow(mat.sno) > 0) {
  	type.sno.ok <- FALSE

    if (toupper(type.sno) == "CALONNE") {
  	  mat.sno <- Calonne(mat.sno)
  	  type.sno.ok <- TRUE
    }
  	
    if (toupper(type.sno) == "YEN") {
  	  mat.sno <- Yen(mat.sno)
  	  type.sno.ok <- TRUE
    }  
    
    if (toupper(type.sno) == "COSENZA") {
    	  mat$k.soi <- 1 #dummy value
   	  mat.sno <- Cosenza(mat.sno)
  	  type.sno.ok <- TRUE
    }    
  
    if (toupper(type.sno) == "STURM") {
      mat.sno <- Sturm(mat.sno)
  	  type.sno.ok <- TRUE
    }     
    
    if (toupper(type.sno) == "JORDAN") {
      mat.sno <- Jordan(mat.sno)
  	  type.sno.ok <- TRUE
    }      

    if (toupper(type.sno) == "WILLIAMS") {
      mat.sno <- Williams(mat.sno)
  	  type.sno.ok <- TRUE
    }
    
    #feedback if function type was not recognized
    if (type.sno.ok == FALSE) stop("'type.sno' not recognized") 
  }
  
  
  #process ground layers ====================================================
  mat.gnd <- mat[-layers.sno,]
  if (layers.sno[1] <= 0) mat.gnd <- mat
  if (nrow(mat.gnd) > 0) {
  	type.gnd.ok <- FALSE
  	
  	if (toupper(type.gnd) == "COSENZA") {
   	  mat.gnd <- Cosenza(mat.gnd)
  	  type.gnd.ok <- TRUE
    } 
    
    if (toupper(type.gnd) == "GEOMETRIC") {
   	  mat.gnd <- Geometric(mat.gnd)
  	  type.gnd.ok <- TRUE
    } 
    
    if (toupper(type.gnd) == "ARITHMETIC") {
   	  mat.gnd <- Arithmetic(mat.gnd)
  	  type.gnd.ok <- TRUE
    } 
    
    if (toupper(type.gnd) == "HARMONIC") {
   	  mat.gnd <- Harmonic(mat.gnd)
  	  type.gnd.ok <- TRUE
    } 

    #feedback if function type was not recognized
    if (type.gnd.ok == FALSE) stop("'type.gnd' not recognized") 
  }
 
  return(rbind(mat.sno,mat.gnd))
}


#--- Calonne et al. (2011) --------------------------------------------------
Calonne <- function(mat) {
  #Calonne, N., Flin, F., Morin, S., Lesaffre, B., du Roscoat, S. R., & 
  #Geindreau, C. (2011). Numerical and experimental investigations of the 
  #effective thermal conductivity of snow. Geophysical Research Letters, 
  #38(23), n/a–n/a. doi:10.1029/2011GL049234
  rho.i <-  917  # [kg m-3]
  rho.l <- 1000  # [kg m-3]
  dens <- mat$ice*rho.i + mat$liq*rho.l
  mat$kj <- 2.5e-06 * dens^2 - 1.23e-04 * dens + 0.024 # Equation 12
  return(mat)
}
    
#--- Yen (1981) -------------------------------------------------------------
Yen <- function(mat) {
  #Yen, Y.-C. (1981). Review of thermal properties of Review of thermal 
  #properties of snow, ice and sea ice (p. 34). Hanover, NH, USA.
  k.ice <- 2.24  # [W m-1 K-1]
  rho.i <-  917  # [kg m-3]
  rho.l <- 1000  # [kg m-3]
  dens <- mat$ice*rho.i + mat$liq*rho.l
  mat$kj <- k.ice * (dens/rho.l)^1.885 # Equation 34
  return(mat)
}  
  
#--- Sturm at al. (1997) -----------------------------------------------------
Sturm <- function(mat) {
  #Sturm, M., J. Holmgren, M. König, and K. Morris (1997), The thermal 
  #conductivity of seasonal snow, J. Glaciol., 43(143), 26–41.
  mat$kj <- 10.0^( 2.650 / 1000.0 * dens - 1.652)
  
  	# if(density<0.156){
	# kt=0.023+0.234*density;
  # }else{
			# kt=0.138-1.01*density+3.233*pow(density,2.0);
	# }
  return(mat)
}     
    
#--- Jordan at al. (1999) -----------------------------------------------------
Jordan <- function(mat) {
  #Jordan, R. E., Andreas, E. L., & Makshtas, A. P. (1999). Heat budget of 
  #snow-covered sea ice at North Pole 4. Journal of Geophysical Research, 
  #104(C4), 7785. doi:10.1029/1999JC900011
  k.ice <- 2.24  # [W m-1 K-1]
  k.air <- 0.025 # [W m-1 K-1]
  rho.i <-  917  # [kg m-3]
  rho.l <- 1000  # [kg m-3]
  dens <- mat$ice*rho.i + mat$liq*rho.l
  mat$kj <- k.air + (7.75e-05 *dens + 1.105e-6 * dens^2) * (k.ice - k.air)
  return(mat)
}    
    
#--- Williams (1989) ----------------------------------------------------------
Williams <- function(mat) {
  #Figure 4.11 in "The Frozen Earth: fundamentals of geocryology" by P. 
  #Williams and M. Smith, 1989. Referenced to: 
  #Goodrich, L. E. (1982). The influence of snow cover on the ground thermal 
  #regime. Canadian Geotechnical Journal, 19(4), 421–432. doi:10.1139/t82-047
  rho.i <-  917  # [kg m-3]
  rho.l <- 1000  # [kg m-3]
  dens <- mat$ice*rho.i + mat$liq*rho.l
  mat$kj <- 2.9e-06 * dens^2
  return(mat)
}  

#--- Cosenza at al. (2003) ---------------------------------------------------
Cosenza <- function(mat) {
  #Cosenza, P., Guerin, R., & Tabbagh, A. (2003). Relationship between thermal 
  #conductivity and water content of soils using numerical modelling. European 
  #Journal of Soil Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x	
  k.liq <- 0.56  # [W m-1 K-1]
  k.ice <- 2.24  # [W m-1 K-1]
  k.air <- 0.025 # [W m-1 K-1]
  air <- 1 - mat$liq - mat$ice - mat$soi
  mat$kj <- (mat$liq * sqrt(k.liq) +  
             mat$ice * sqrt(k.ice) +
             mat$soi * sqrt(mat$k.soi) +
                 air * sqrt(k.air))^2     	
  return(mat)
}   

#--- Geometric mean ---------------------------------------------------
Geometric <- function(mat) {	
  k.liq <- 0.56  # [W m-1 K-1]
  k.ice <- 2.24  # [W m-1 K-1]
  k.air <- 0.025 # [W m-1 K-1]
  air <- 1 - mat$liq - mat$ice - mat$soi
  mat$kj <- k.liq^mat$liq *  
            k.ice^mat$ice * 
            k.air^air *
            mat$k.soi^mat$soi   	
  return(mat)
}  

#--- Arithmetic mean ---------------------------------------------------
Arithmetic <- function(mat) {	
  k.liq <- 0.56  # [W m-1 K-1]
  k.ice <- 2.24  # [W m-1 K-1]
  k.air <- 0.025 # [W m-1 K-1]
  air <- 1 - mat$liq - mat$ice - mat$soi
  mat$kj <- k.liq*mat$liq + 
            k.ice*mat$ice + 
            k.air*air +
            mat$k.soi*mat$soi   	
  return(mat)
}  

#--- Harmonic mean ---------------------------------------------------
Harmonic <- function(mat) {	
  k.liq <- 0.56  # [W m-1 K-1]
  k.ice <- 2.24  # [W m-1 K-1]
  k.air <- 0.025 # [W m-1 K-1]
  air <- 1 - mat$liq - mat$ice - mat$soi
  mat$kj <- 1/(mat$liq/k.liq + 
               mat$ice/k.ice + 
               air    /k.air +
               mat$soi/mat$k.soi)   	
  return(mat)
}  