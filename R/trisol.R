# =============================================================================
#'
#' @title Trisol
#'
#' @description Solves tri-diagonal matrix equations.
#'
#' @details Solves matrix with Thomas algorithm, Z discretizations with depth. 
#'
#' @param a Lower diagonal term [Z], Note: a[,0] is zero
#' @param b Center diagonal term [Z]
#' @param c Upper diagonal term [Z], Note: c[.Z-1] is zero
#' @param d Right-hand side (on input) / solution (on output) [Z]
#' 
#' @return Solution vector result. In case of heat conduction this is
#'	       the temperature for the new time step. 
#' 
#' @export
#' @examples
#' a <- c(0, 1, 1, 1)
#' b <- c(2, 2, 2, 2)
#' c <- c(1, 1, 1, 0)
#' d <- c(1, 1, 1, 1)
#' d <- Trisol(a, b, c, d)
#' 
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================

Trisol <- function(a, b, c, d) {     
  #get dimensions
  m <- length(a)

  #establish upper tridiagonal matrix
  for (i in 2:m) { #loop over z-discretizations
    r <- a[i] / b[i-1]
    b[i] <- b[i] - r * c[i-1]
    d[i] <- d[i] - r * d[i-1]
  }

  #back substitution
  #example m=4
  #IDL 0,1,2,3 i: 1,2,3 j: 2,1,0
  #R   1,2,3,4 i: 2,3,4 j: 3,2,1
  d[m] <- d[m] / b[m]
  for (i in 2:m) { #loop over z-discretizations
    j <- m-i+1
    d[j] <- (d[j] - c[j] * d[j+1]) / b[j]
  }  
  return(d)
} 
