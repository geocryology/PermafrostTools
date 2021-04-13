# =============================================================================
#'
#' @title Return an empty data frame containing loc_name, Snow.Start and Snow.End
#'
#' @description Creates an empty data frame containing loc_name, Snow.Start and Snow.End.
#'              The Loc.name is inserted.
#'
#' @details 
#'
#' @param obs_stat data frame with a index-column
#'      
#' 
#'
#' @return a dataframe with loc_name, Snow.Start=NA, Snow.End=NA
#' 
#' @keywords internal          
#' 
#' @export 
#' @examples
#' table <- TSA_row(location_name)
#' 
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_row <- function(location_name){

table <- as.data.frame(matrix(ncol = 3, nrow = 1))
names(table) <- c("loc_name","Snow.Start","Snow.End")

table$loc_name <-   location_name  #add location name
        
table$Snow.Start <- NA

table$Snow.End <- NA

return(table)

}