# =============================================================================
#'
#' @title Discretize soil vertically
#'
#' @description Generate thickness and depth of soil layers for numerical model  
#'
#' @details Makes a soil layer spacing based on the simple model 
#'          dz = dzmin * base^(n) where n is the layer number.
#'
#' @param dzmin Minimal z-spacing [m] ("how fine is the grid") 
#' @param zmax Depth of lowermost node [m] center ("how large is the domain")
#' @param base Resolution reduction ("how strong is the grid corsened with depth"); 
#'             base = 1 corresponds to equal spacing, base = 1.5 implies each 
#'             successive layer with depth begin 1.5 times thicker than the one 
#'             above. 
#' 
#' @return data frame with layer thickness (dz) and depth of layer centre (z)
#'         below the surface   
#' 
#' @export
#' @examples
#' soil <- SoilDiscretize(0.1, 20, 1.25)
#' 
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================

SoilDiscretize <- function(dzmin, zmax, base) {
	if (base < 1) stop("Parameter 'base' nneds to be >= 1")
	#layer thicknesses	
	dz <- dzmin * base^(0:2005)
	
	#depth of layer lower boundary
	z_bot <- cumsum(dz)
	
	#depth of layer upper boundary
	z_top <- z_bot - dz
	#depth of layer center (node)
	z <- (z_bot + z_top) / 2
	
	#data frame
	discr <- data.frame(dz=dz,z=z,z_bot=z_bot,z_top=z_top)
	
	#restrict to maximum depth
	discr<-subset(discr,z < zmax)
	nz<-length(discr$z)
	
	discr$dz[nz] <- zmax-discr$z_bot[nz-1] 	
	discr$z[nz]  <- discr$z_bot[nz-1] + discr$dz[nz]/2    
	
	#drop auxiliary columns
	discr<-discr[1:nz,1:2] 
	return(discr)
}
