# =============================================================================
#'
#' @title Heat capacity 
#'
#' @description Apparent heat capacity of soil or rock
#'
#' @details Computes the apparent heat capacity of soil or rock based on
#'          the heat capacity of the consituents and the apparent heat 
#'          capacity parameterized using the derivative of the liquit 
#'          water content (see funtion "Unfrozen"). 
#'
#' @param mat Data frame with relevant ground properties. Needed: mat$ice (ice
#'             content, volumetric fraction) and mat$liq (liquid water
#'             content, volumetric fraction).  
#' 
#' @return Returns the input data frame with the column cj updated in [J m-3 K-1]
#'         If this column does not exist in the input it is created.   
#' 
#' @export
#' @examples
#' mat <- HeatCapacity(mat)
#'
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================
HeatCapacity <- function(mat) {
  c.liq <- 4180000  # volumetric heat capacity of water [J m-3 K-1]
  c.ice <- 1925700  # volumetric heat capacity of ice [J m-3 K-1
  c.air <-    1200  # volumetric heat capacity of air [J m-3 K-1]
  Lf    <-   3.3e8  # volumetric latent heat of fusion [J m-3]
  mat$air <- 1 - mat$liq - mat$ice - mat$soi
  
  mat$cj <- c.liq * mat$liq + 
            c.ice * mat$ice + 
            c.air * mat$air +
            mat$c.soi * mat$soi 

  #apparent part
  mat$cj <- mat$cj + mat$dice * Lf * mat$wat
      	 
  return(mat)
}

