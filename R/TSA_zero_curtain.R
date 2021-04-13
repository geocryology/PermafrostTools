# =============================================================================
#'
#' @title Detect zero curtain periods during a snow cover period
#'
#' @description Give one data.frame of a location or a list of data.frames of multiple 
#'              locations and the function returns a table with the start and end dates of zero curtain
#'              periods during snow cover periods.
#'
#' @details It is also possible to adjust the parameters for the TSA_snow_cover function in this function, 
#'          as the warming periods are calculated for each snow cover period. 
#'          Therefore, the TSA_snow_cover is run as well in this function.
#'
#' @param obs_stat one dataframe of one location of a list of dataframes of multiple locations
#'                 A dataframe needs at least the following columns:
#'                 Location name as loc_name;
#'                 Mean daily GST as agg_avg;
#'                 Maximum daily GST as agg_max;
#'                 Daily sd as agg_sd (not necessarily needed);
#'                 Date as time (in POSIXct and the following format: YYYY-MM-DD)
#'      
#' @param temp    The temperature boundary within the mean agg_avg of a zero curtain period has to be. 
#'                Default: 0.2 degrees Celsius     
#' 
#' @param slopesd   Fraction of the standard deviation of the slope used as threshold for the slope. Default: 0.1
#'             
#' 
#' @param tempsd     Fraction of the mean daily standard deviation used as threshold for the daily sd. Default: 0.01
#'               
#' @param lengthzc     Length that a zero curtain period should at least have. Default: 2 days.
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
#' @param snow.period   In case the snow period is already available for one location, it can be inserted here.
#'                      The data frame of the now period needs to have the following columns: 
#'                      loc_name, Snow.Start, Snow.End, MDr
#'
#' @return Data fame with the following fields: 
#'         loc_name: location name;
#'         Snow.Start: Start date of the analyzed snow period; 
#'         Snow.End: End date of the analyzed snow period;
#'         Zero.Curtain.Start: Start date of a zero curtain period;
#'         Zero.Curtain.End: End date of a zero curtain period;
#'         Median.Temperature: Median temperature of a zero curtain period [C];
#'         sd.Temperature: Standard deviation of the temperature of a zero curtain period;
#'         
#'         
#' 
#' 
#' 
#' @export 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#' ZeroCurtainTable <- TSA_zero_curtain(BC_SO01_01, temp= 0.2, slopesd= 0.1, tempsd= 0.01, lengthzc=2, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, snow.period=NULL)
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_zero_curtain <- function(obs_stat, temp= 0.2, slopesd= 0.1, tempsd= 0.01, lengthzc= 2, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, snow.period=NULL){
  if(class(obs_stat)=="data.frame"){obs_stat <- list(obs_stat)}
  
  #create an error dataframe
  error.frame <- as.data.frame(matrix(ncol = 7, nrow = 1))
  names(error.frame) <- c("loc_name","Snow.Start","Snow.End","Zero.Curtain.Start","Zero.Curtain.End","Median.Temperature","sd.Temperature")
  
  
  #run the f.zero.curtain1 function over the whole dataset
  zero.curtain.list <- lapply(obs_stat, FUN = function(X) tryCatch(f.zero.curtain1(X, temp, slopesd, tempsd, lengthzc, v1, v2, lengthsnow, MDr.sd, snow.period), error=function(e) error.frame))
  
  #create a dataframe from the output list
  zero.curtain <- do.call(rbind, zero.curtain.list)
  
  #delete the Row.names column
  zero.curtain$Row.names <- NULL

  
  return(zero.curtain)
  
}


#function to calculate all the zero curtain periods for a single location
f.zero.curtain1 <- function(obs_stat, temp, slopesd, tempsd, lengthzc, v1, v2, lengthsnow, MDr.sd, snow.period){
  
  if(is.null(snow.period)){snow.period <- TSA_snow_cover(obs_stat, v1, v2, lengthsnow, MDr.sd)}
  
  # create table
  table <- TSA_row(obs_stat$loc_name[1])
  
  table$Zero.Curtain.Start <- NA
  
  table$Zero.Curtain.End <- NA
  
  table$Median.Temperature <- NA
  
  table$sd.Temperature <- NA
  
  if(is.na(snow.period$Snow.Start[1])){return(table)}
  
  
  # calculating slope of the maximum temperature
  obs_stat$slope <- append(diff(obs_stat$agg_max), 0, after = 0)
  
  # subset the dataframe using the snow periods
  obs_stat$index <- 0
  
  for(i in 1:length(snow.period$loc_name)){
    
    obs_stat$index[snow.period$Snow.Start[i] <= obs_stat$time & obs_stat$time <= snow.period$Snow.End[i]] <- 1
    
  }
  
  for(ii in 1:length(obs_stat$index)){
    
    if(obs_stat$index[ii]==1){obs_stat$index[ii]<-obs_stat$time[ii]}
    
  }
  
  snow <- TSA_subsetting(obs_stat)
  
  
  #apply calculation function on the snow periods
  table1 <- lapply(snow, function(X) tryCatch(f.zero.curtain.calc(X,temp, slopesd, tempsd, lengthzc, table),error=function(e) table))
  
  #create data frame
  zero.curtain.table <- do.call(rbind, table1)
  
  #zero.curtain <- zero.curtain.table
  zero.curtain <- zero.curtain.table[complete.cases(zero.curtain.table), ]
  
  #make sure all times are in correct format
  zero.curtain$Zero.Curtain.Start <- as.POSIXct(zero.curtain$Zero.Curtain.Start, origin="1970-01-01")
  zero.curtain$Zero.Curtain.End <- as.POSIXct(zero.curtain$Zero.Curtain.End, origin="1970-01-01")
  zero.curtain$Snow.Start <- as.POSIXct(zero.curtain$Snow.Start, origin="1970-01-01")
  zero.curtain$Snow.End <- as.POSIXct(zero.curtain$Snow.End, origin="1970-01-01")
  
  return(zero.curtain)
}
 

#function to calculate the zero curtains for one snow period
f.zero.curtain.calc <- function(obs_stat, temp, slopesd, tempsd, lengthzc, table){
  
  table$Snow.Start <- head(obs_stat$time,n=1)
  
  table$Snow.End <- tail(obs_stat$time,n=1)

  
  # define threshold for the slope using the standard deviation
  standarddev <- slopesd*sd(obs_stat$slope)
  
  obs_stat$slope1 <-0
  obs_stat$slope1[obs_stat$slope<standarddev & obs_stat$slope > -standarddev]<-1
  
  
  
  # temperature standard deviation index
  obs_stat$sd <- 0
  obs_stat$sd[obs_stat$agg_sd<tempsd] <- 1
  
  
  # combine the standard deviation index and slope index
  obs_stat$sd <- obs_stat$sd * obs_stat$slope1
  
  # subset all the zero curtain periods according to the zero curtain index
  obs_stat$index <- obs_stat$agg_sd * obs_stat$sd
  
  zero.curtains <- TSA_subsetting(obs_stat)
  
  if(length(zero.curtains)==0){return(table)}
  
  # select by maximum temperature that is reached 
  zero.curtains <- lapply(zero.curtains, function(X) f.select.by.maxtemp1(X,temp))
  zero.curtains <- zero.curtains[!sapply(zero.curtains, is.null)]
  
  if(length(zero.curtains)==0){return(table)}
  
  # select by length
  zero.curtains <- lapply(zero.curtains, function(X) f.select.by.length(X, lengthzc))
  zero.curtains <- zero.curtains[!sapply(zero.curtains, is.null)]
  
  if(length(zero.curtains)==0){return(table)}
  
  # create a row with all the important information for each zero curtain period
  table <- lapply(zero.curtains, function(X) f.select.start.end.date(X,table))
  
  table <- do.call(rbind.data.frame, table)
  
  return(table)
  
}

  






# function to select zero curtains by lenth
f.select.by.length <- function(obs_stat, lengthzc){
  
  if(length(obs_stat$loc_name)>lengthzc){return(obs_stat)}
  else{obs_stat <- NULL}
  
}


# check that the warming period reaches a certain temperature limit
f.select.by.maxtemp1 <- function(obs_stat,temp){
  
  if(mean(obs_stat$slope)<0){
    
    if(mean(obs_stat$agg_avg)>-temp & mean(obs_stat$agg_avg)<temp) {return(obs_stat)}
    else{obs_stat <- NULL}
  }
  
  else{
    if(tail(obs_stat$agg_avg,n=1)>-0.3 & tail(obs_stat$agg_avg,n=1)<0.1) {return(obs_stat)}
    else{obs_stat <- NULL}
    
  }
  
  
  
  
}






# function to create a row with all important informations for each zero curtain period
f.select.start.end.date <- function(obs_stat, table){
  
  start.date <- head(obs_stat$time,n=1)-1
  
  start.date <- as.POSIXct(substring(start.date,1,10))
  
  table$Zero.Curtain.Start <- start.date
  
  table$Zero.Curtain.End <- tail(obs_stat$time,n=1)
  
  table$Median.Temperature <- round(median(obs_stat$agg_avg),digits=3)
  
  table$sd.Temperature <- round(sd(obs_stat$agg_avg),digits=4)
  
  return(table)
}

