# =============================================================================
#'
#' @title Find start and end date of snow cover 
#'
#' @description Give one data.frame of a location or a list of data.frames of multiple 
#'              locations and the function returns the start and end date as well as  
#'              the MDr in a output table
#'
#' @details At locations with perennial snow cover, the snow cover is splitted into water years periods.
#'
#' @param obs_stat One dataframe of one location or a list of dataframes of multiple locations.
#'                 A dataframe needs at least the following columns:
#'                 Location name as loc_name;
#'                 Mean daily GST as agg_avg;
#'                 Minimum daily GST as agg_min;
#'                 Maximum daily GST as agg_max;
#'                 Daily sd as agg_sd;
#'                 Date as time (in POSIXct and the following format: YYYY-MM-DD)
#'      
#' 
#' @param v1    Max. daily standard deviation indicating snow for POSITIVE Ground Surface Temperature. 
#'              Default is 0.1.
#'              
#' @param v2    Max. daily standard deviation indicating snow for NEGATIVE Ground Surface Temperature. 
#'              Default is 0.3.
#'              
#' @param lengthsnow    Number of days that the snow cover should at least have to be selected as a valid snow cover.
#'                     Default is 5 days
#'
#' @param MDr.sd   Threshold for the mean daily standard deviation for a snow period. 
#'
#' @return Dataframe with fields: loc_name; Snow.Start; Snow.End; MDr
#'           
#' 
#' @export 
#' 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#' Snow_table <- TSA_snow_cover(BC_SO01_01, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3)
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_snow_cover <- function(obs_stat, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3){
  if(class(obs_stat)=="data.frame"){obs_stat <- list(obs_stat)}
  
  #create an error dataframe
  error.frame <- as.data.frame(matrix(ncol = 4, nrow = 1))
  names(error.frame) <- c("loc_name","Snow.Start","Snow.End","MDr")
  
  
  #run the snow_cover1 function over the whole dataset
  snow.cover.list <- lapply(obs_stat, FUN = function(X) tryCatch(f.snow.cover1(X, v1, v2, lengthsnow, MDr.sd), error=function(e) error.frame))
  
  #create a dataframe from the output list
  snow.cover <- do.call(rbind, snow.cover.list)
  
  #delete the Row.names column
  snow.cover$Row.names <- NULL

  return(snow.cover)
}

#function that calculates the snow cover for one location
f.snow.cover1 <- function(obs_stat, v1, v2, lengthsnow, MDr.sd){
  
  "&" <- function(...) UseMethod("&")
  "&.default" <- .Primitive("&")
  "&.character" <- function(...) paste(...,sep="")
  
  #
  m1 = 3     #maximum temperature where snow appears
  m2 = 0.5   #temperature where snow melt appears
  

  # create table
  table <- TSA_row(obs_stat$loc_name[1])
  
  table$MDr <- NA
  
  
  if(length(obs_stat$agg_avg)==0){return(table)}

  
  ### snow index
  obs_stat$snow.ind[obs_stat$agg_sd < ifelse(obs_stat$agg_max<0,yes=v2,no=v1)] <- 1
  obs_stat$snow.ind[obs_stat$agg_max >m1]<-0
  obs_stat$snow.ind[is.na(obs_stat$snow.ind)] <- 0
  
  ### temp index
  obs_stat$temp.ind[obs_stat$agg_max<m2]<-1
  obs_stat$temp.ind[is.na(obs_stat$temp.ind)] <- 0
  
  ## gaps with temperatures lower 0.5deg get filled up
  for(ii in 1:(length(obs_stat$snow.ind)-1)){
    if(obs_stat$snow.ind[ii] == 1 & obs_stat$snow.ind[ii + 1] == 0 & obs_stat$temp.ind[ii + 1] == 1) obs_stat$snow.ind[ii+1] <- 1
  }

  
  #gaps filling for max temp outliers
  for(a in 2:(length(obs_stat$snow.ind)-1)){
    
    if(obs_stat$agg_max[a] > 0 & obs_stat$agg_max[a-1] < 0 & obs_stat$agg_max[a+1] < 0 & obs_stat$agg_avg[a] < 0){obs_stat$snow.ind[a] <- 1}
    
  }
  
  
   
  # creating subsets of all snow periods
  obs_stat$index <- 0
  for(i in 1:length(obs_stat$index)){
    
    if(obs_stat$snow.ind[i]==1){obs_stat$index[i]<-obs_stat$time[i]}
    
  }
  
  snow.periods.list <- TSA_subsetting(obs_stat)
  
  #select snow periods by MDr
  snow.periods <- lapply(snow.periods.list, function(X) f.select.By.MDr(X, MDr.sd))
  snow.periods <- snow.periods[!sapply(snow.periods, is.null)]
  
  if(length(snow.periods)==0){table$MDr <- "MDr < 0";return(table)}
  
  
  # select by length
  snow.periods <- lapply(snow.periods, function(X) f.select.by.length.snow(X, lengthsnow))
  snow.periods <- snow.periods[!sapply(snow.periods, is.null)]
  
  if(length(snow.periods)==0){table$MDr <-"too short" ;return(table)}
  
  
  snow.periods <- lapply(snow.periods, function(X) split_perennial_snow(X))
  
  snow.periods <- unlist(lapply(snow.periods, function(x) if (class(x) == "data.frame") list(x) else x), recursive=FALSE)
  
  #extract start and end date and the MDr 
  table <- lapply(snow.periods, function(X) f.select.start.end.date.snow(X,table))
  
  table <- do.call(rbind.data.frame, table)
  
  #table$Snow.Start <- as.POSIXct(table$Snow.Start, origin="1970-01-01")
  #table$Snow.End <- as.POSIXct(table$Snow.End, origin="1970-01-01")
  
  
  return(table)
}

# function to exclude all snow periods with a too small MDr
f.select.By.MDr <- function(snow_period, MDr.sd){

  #calculate mean of the mean daily sd
  sd.day  <- mean(snow_period$agg_sd)
  
  MDr <- MDr.sd - sd.day
  
  snow_period$MDr <- round(MDr,digits=3)
  
  if(MDr>=0){return(snow_period)}
  else{snow_period <- NULL}
  
  
}


# function to select snow periods longer than lengthsnow 
f.select.by.length.snow <- function(snow_period, lengthsnow){
  
  if(length(snow_period$loc_name)>5){return(snow_period)}
  else{snow_period <- NULL}
  
}



# function to create a row with all important informations for each snow cover period
f.select.start.end.date.snow <- function(snow_period, table){
  
  table$Snow.Start <- head(snow_period$time,n=1)
  
  table$Snow.End <- tail(snow_period$time,n=1)
  
  table$MDr <- head(snow_period$MDr, n=1)
  
  
  return(table)
}


split_perennial_snow <- function(obs_stat){
  
  year.start <-  as.numeric(format(as.Date(head(obs_stat$time,n=1)), "%Y"))
  year.end <- as.numeric(format(as.Date(tail(obs_stat$time,n=1)), "%Y"))
  
  month.start <-  as.numeric(format(as.Date(head(obs_stat$time,n=1)), "%m"))
  month.end <- as.numeric(format(as.Date(tail(obs_stat$time,n=1)), "%m"))
  
  
  yeardiff <- year.end-year.start
  
  if(yeardiff == 0){return(obs_stat)}
  if(yeardiff == 1 & month.end < month.start){return(obs_stat)}else{
    
    dates.posix = as.POSIXlt(obs_stat$time)
    # Year offset
    offset = ifelse(dates.posix$mon >= 10 - 1, 1, 0)
    # Water year
    obs_stat$wateryear = dates.posix$year + 1900 + offset
    
    
    new.list <- as.list(
      by(obs_stat, obs_stat$wateryear, function(X){
        
        X$wateryear <- NULL
        X
        
      }))
    
    return(new.list)  
    
  }
  
  
}
