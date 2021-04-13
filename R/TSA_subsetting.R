# =============================================================================
#'
#' @title Return list with subset dataframes
#'
#' @description Give one data.frame of an index-column containing unique values. It returns a list of 
#'              the subset-dataframes
#'
#' @details 
#'
#' @param obs_stat data frame with a index-column
#'      
#' 
#'
#' @return a list containing all the subset-dataframes
#' 
#' @keywords internal          
#' 
#' @export
#' @examples
#' snow_periodlist <- TSA_subsetting(obs_stat)
#' 
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_subsetting <- function(obs_stat){
  v <- obs_stat$index
  
  r <- rle(v > 0)
  r <- r$lengths[r$values]
  
  (pos <- v[v > 0])
  
  
  dates <- lapply(r, function(x) {
    out <- pos[1:x]
    pos <<- pos[-(1:x)]
    out
  })
  
  dates <- lapply(dates, function(x) f.select(x,obs_stat))
  
  
  return(dates)
}


# function to get the data out of the obs_stat data frame for each major temperature increase
f.select <- function(dates, obs_stat){
  incdates <- obs_stat[obs_stat$index %in% dates,]
  return(incdates)
}