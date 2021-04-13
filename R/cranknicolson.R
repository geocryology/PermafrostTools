# =============================================================================
#'
#' @title Cranknicolson
#'
#' @description Implicit solution of the 1D heat conduction equation. 
#'
#' @details 1-dimensional finite-difference (crank-nicolson) heat conduction
#'	        scheme that operates on variable depth-discretization. The 
#'          boundary conditions can be selected as either Dirichlet-type
#'          (prescribed temperature) or Neumann-type (prescribed energy flux).
#'          
#'
#' @param Tj  Temperature [K or C], array with dimensions (Z)
#' @param cj  Volumetric heat capacity [J m-3 K-1], array with dimensions (Z) 
#' @param kj  Thermal conductivity [W m-1 K-1], array with dimensions (Z)
#' @param gj  Volumetric heat source [W], array with dimensions (Z)
#'             This is ignored for boundary nodes with a Dirichlet condition.
#' @param zj  Z-levels increasing upwards [m], array with dimensions (Z)
#' @param dt  Time step [s]
#' @param bcl Lower boundary condition as either [W m-2] (Neumann) or [C or K] 
#'            (Dirichlet).   
#' @param bcu Upper boundary condition as either [W m-2] (Neumann) or [C or K] 
#'            (Dirichlet). 
#' @param bcutype Type of upper boundary condition used. Can be "NEUMANN" 
#'                (standard, fixed heat flux) or "DIRICHLET" (fixed temperature).  
#' @param bcltype Type of upper boundary condition used. Can be "NEUMANN"
#'                (standard, fixed heat flux) or "DIRICHLET" (fixed temperature).  
#' 
#' @return Returns An array (Z) of temperatures [C or K] resulting from
#'         forcing the system one time step further with the given bounday 
#'         conditions.   
#' 
#' @export
#' @examples
#' Tj <- c(  1,    1,    1,    1)
#' kj <- c(  2,    2,    2,    2)
#' cj <- c(2e6,  2e6,  2e6,  2e6)
#' gj <- c(  0,    0,    0,    0)
#' zj <- c(  0,   -1,   -2,   -3)
#' dt <- 30000
#' bcl <-  1
#' bcu <- -1
#' 
#' #set up plot
#' plot(Tj, zj, type="l", xlim=c(-1,2))
#' 
#' #loop over heat conduction and add new lines
#' for (n in 1:100) {
#'   Tj <- Cranknicolson(Tj, kj, cj, gj, zj, dt, bcl, bcu, 
#'                       bcutype = "DIRICHLET", bcltype = "NEUMANN")
#'   lines(Tj, zj)
#' }
#' 
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#' TODO: introduce heat transfer coefficient for boundary condition
# =============================================================================

Cranknicolson <- function(Tj, kj, cj, gj, zj, dt, bcl, bcu, bcutype = 
                          "NEUMANN", bcltype = "NEUMANN") {

  #PREPARATION ========================================================
  m <- length(Tj) # number of z-discretizations

  #COMPUTE MATRIX TERMS FROM SUBSTITUTIONS ============================
  #volumetric heat capacity divided by time step
  cdt <- cj / dt
       
  #z-difference of upper and lower node (artificial negative sign) 
  dz2      <- ColShiftM1(zj) - ColShiftP1(zj)                                 
  
  dz2[m] <- zj[m] - zj[m-1] #lowest  node
  dz2[1] <- zj[2] - zj[1]   #highest node

  #mean k over delta-z for upper space
  kdu    = (ColShiftP1(kj)  + kj) / (ColShiftP1(zj) - zj) / 2
  kdu[1] = kj[1] / (zj[1] - zj[2])

  #Upper diagonal term [N,Z], Note: C[*,m] fixed for BC
  c = ColShiftM1(kdu) / dz2

  #Lower diagonal term [N,Z], Note: a[ ,0]=0
  a = kdu / dz2

  #Center diagonal term [N,Z]
  b <- cdt - a - c

  #Right-hand side (on input) [N,Z]
  d <- (cdt + a + c) * Tj -  #node i,  center diagonal
       c * ColShiftM1(Tj) -  #node i+1, upper diagonal
       a * ColShiftP1(Tj) -  #node i-1, lower diagonal
       gj / dz2 * 2          #source term (normalized to W)

  #LOWER BOUNDARY CONDITION (NEUMANN OR DIRICHLET) ===================
  if (toupper(bcltype) == "DIRICHLET") {
  	 #lower boundary condition = Dirichlet
     d[m-1] <- d[m-1] - c[m-1] * bcl 
     c[m-1] <- 0 
  } else {
  	#lower boundary condition = Neumann        
    c[m] <- 0
    b[m] <-  cdt[m] - a[m]
    d[m] <- (cdt[m] + a[m]) * Tj[m] -        #node i, center diagonal
                          a[m]  * Tj[m-1] +  #node i+1, upper diagonal
            (2.0 * (-bcl) + gj[m]) / dz2[m]  #boundary condition
  }

  #SURFACE BOUNDARY CONDITION (NEUMANN OR DIRICHLET) =================
  if (toupper(bcutype) == "DIRICHLET") {
  	 #upper boundary condition = Dirichlet
     d[2] <- d[2] - a[2] * bcu 
     a[2] <- 0
  } else {
  	 #upper boundary condition = Neumann        
     a[1] <- 0 
     b[1]  =  cdt[1] - c[1]
     d[1]  = (cdt[1] + c[1]) * Tj[1] -   #node i,  center diagonal
                       c[1]  * Tj[2] +   #node i+1, upper diagonal
             (2.0 * (-bcu) + gj[1]) / dz2[1]  #boundary condition
  }

  #simultaneously solve tridiagonal systems for all n points/pixels ===
  Tj <- Trisol(a, b, c, d)

  #insert Dirichlet boundary condition if needed and return
  if (toupper(bcutype) == "DIRICHLET") Tj[1] <- bcu
  if (toupper(bcltype) == "DIRICHLET") Tj[m] <- bcl

  return(Tj)
}


# helper function: shift array by plus 1, like shift(zj, 0, 1)
ColShiftP1 <- function(arr) {
	nc <- length(arr)
	return(c(arr[nc],arr[1:(nc-1)]))
}

# helper function: shift array by minus 1, like shift(zj, 0, -1)
ColShiftM1 <- function(arr) {
	nc <- length(arr)
	return(c(arr[2:nc],arr[1]))
}
