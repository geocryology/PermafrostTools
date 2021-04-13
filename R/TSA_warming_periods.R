# =============================================================================
#'
#' @title Detect major warming periods during a snow cover period
#'
#' @description Give one data.frame of a location or a list of data.frames of multiple 
#'              locations and the function returns a table with the start dates of the major warming periods  
#'              during a winter as well as quality values.
#'
#' @details It is also possible to adjust the parameters for the TSA_snow_cover function in this function, 
#'          as the warming periods are calculated for each snow cover period. 
#'          Therefore, the TSA_snow_cover is run as well in this function.
#'
#' @param obs_stat One dataframe of one location or a list of dataframes of multiple locations
#'                 A dataframe needs at least the following columns:
#'                 Location name as loc_name;
#'                 Mean daily GST as agg_avg;
#'                 Minimum daily GST as agg_min;
#'                 Maximum daily GST as agg_max;
#'                 Daily sd as agg_sd; 
#'                 Date as time (in POSIXct and the following format: YYYY-MM-DD)
#'      
#' @param tperc    Percentage of the 10-percent-quantile of the winter temperatures that the 
#'                     warming period should at least reach. Default is 40 percent.
#' 
#' @param sdoutlier   The number implies how many multiples of the standard deviation of the slope are taken to define a threshold.  
#'             The threshold is used to detect slope outliers.
#' 
#' @param sdsteep    The number implies what fraction of the standard deviation of the slope 
#'               is taken to select the steepest part of a warming period.
#' 
#' @param v1    Max. daily standard deviation indicating snow for POSITIVE Ground Surface Temperature. 
#'              Default is 0.1.
#'              
#' @param v2    Max. daily standard deviation indicating snow for NEGATIVE Ground Surface Temperature. 
#'              Default is 0.3.
#'              
#' @param lengthsnow    Number of days that the snow cover should at least have to be selected as a valid snow cover.
#'                     Default is 5 days.
#'
#' @param MDr.sd   Threshold for the mean daily standard deviation for a snow period. 
#'
#' @param snow.period   In case the snow period is already available for one location, it can be inserted here.
#'                      The data frame of the now period needs to have the following columns: 
#'                      loc_name, Snow.Start, Snow.End, MDr
#'
#'
#' @return Datafame with the following fields: 
#'         loc_name: location name;
#'         Snow.Start: Start date of the analyzed snow period; 
#'         Snow.End: End date of the analyzed snow period;
#'         Warming.Period.Start: start date of the warming period;
#'         RD: indicates if the warming period is assumed to be the Ripening date; 
#'         Temp.diff: Temperature of the steepest part of the warming period [C];
#'         Mean.Slope: Mean slope of the steepest part of the warming period; 
#'         Temp.To.Zero: How much below/above zero is the end of the steepest part of the warming period [C];
#'         
#' 
#' 
#' 
#' @export 
#' 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#' WarmingPeriodsTable <- TSA_warming_periods(BC_SO01_01,tperc = 40, sdoutlier = 3, sdsteep = 0.75, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, snow.period=NULL)
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_warming_periods <- function(obs_stat, tperc = 40, sdoutlier = 3, sdsteep = 0.75, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, snow.period=NULL){
    if(class(obs_stat)=="data.frame"){obs_stat <- list(obs_stat)}
    
    #create an error dataframe
    error.frame <- as.data.frame(matrix(ncol = 8, nrow = 1))
    names(error.frame) <- c("loc_name", "Snow.Start","Snow.End","Warming.Period.Start","RD","Temp.Diff","Mean.Slope","Temp.To.Zero")
  
    #run the f.warming.periods function over the whole dataset
    warming.period.list <- lapply(obs_stat, FUN = function(X) tryCatch(f.warming.periods(X, tperc, sdoutlier, sdsteep, v1, v2, lengthsnow, MDr.sd, snow.period), error=function(e) error.frame))
    
    #create a dataframe from the output list
    warming.period <- do.call(rbind, warming.period.list)
    
    #delete the Row.names column
    warming.period$Row.names <- NULL

    
    return(warming.period)
}



# Function that calculates all the parameters
f.warming.periods <- function(obs_stat, tperc, sdoutlier, sdsteep, v1, v2, lengthsnow, MDr.sd, snow.period){
  
  if(is.null(snow.period)){snow.period <- TSA_snow_cover(obs_stat, v1, v2, lengthsnow, MDr.sd)}

  # create table
  table <- TSA_row(obs_stat$loc_name[1])
  
  table$Warming.Period.Start <- NA
  
  table$RD <- NA
  
  table$Temp.Diff <- NA
  
  table$Mean.Slope <- NA
  
  table$Temp.To.Zero <- NA
  
  table$Mean.Slope.after.warming <- NA
  
  
  # if no snow, then no warming periods selected
  if(is.na(snow.period$Snow.Start[1])){return(table)}
  
  # 4. calculating slope of the maximum temperature
  obs_stat$slope <- append(diff(obs_stat$agg_max), 0, after = 0)
  
  # 5. subset the dataframe using the snow periods
  obs_stat$index <- 0
  
  for(i in 1:length(snow.period$loc_name)){
    
    obs_stat$index[snow.period$Snow.Start[i] <= obs_stat$time & obs_stat$time <= snow.period$Snow.End[i]] <- 1
    
  }
  
  for(ii in 1:length(obs_stat$index)){
    
    if(obs_stat$index[ii]==1){obs_stat$index[ii]<-obs_stat$time[ii]}
    
  }
  
  snow <- TSA_subsetting(obs_stat)
  
  
  #apply calculation function on the snow periods
  table1 <- lapply(snow, function(X) tryCatch(f.wamring.period.calc(X, dataset= obs_stat, tperc, sdoutlier, sdsteep, table),error=function(e) table))
  
  #create data frame
  warming.period.table <- do.call(rbind, table1)
  
  #warming.periods <- warming.period.table
  warming.periods <- warming.period.table[complete.cases(warming.period.table), ]
  
  #make sure that all times ar in POSIXct
  warming.periods$Warming.Period.Start <- as.POSIXct(warming.periods$Warming.Period.Start, origin="1970-01-01")
  warming.periods$Snow.Start <- as.POSIXct(warming.periods$Snow.Start, origin="1970-01-01")
  warming.periods$Snow.End <- as.POSIXct(warming.periods$Snow.End, origin="1970-01-01")
  
  
  return(warming.periods)
}
  
  
  
  




#function to calculate the warming periods for one snow period
f.wamring.period.calc <- function(obs_stat, dataset, tperc, sdoutlier, sdsteep, table){
  
  table$Snow.Start <- head(obs_stat$time,n=1)
  
  table$Snow.End <- tail(obs_stat$time,n=1)
  
  # calculating the max temperature that should be at least reached by the warming period
  quantiles <- quantile(obs_stat$agg_avg, probs = seq(0, 1, by= 0.1))
  
  max_temp <- as.numeric(quantiles[2]/100 * tperc)
  
  # calculating the standard deviation of the slope to find outliers
  standarddev <- sdoutlier*sd(obs_stat$slope)
  
  obs_stat$slope1 <-0
  obs_stat$slope1[obs_stat$slope>standarddev]<-1
  
  # 8. set all negative slopes to 0
  obs_stat$slope[obs_stat$slope<0]<-0
  
  # 9. extract all major warming periods (warming periods that have a at least one value higher than standarddev)
  obs_stat$index <- 0
  for(i in 1:length(obs_stat$index)){
    
    if(obs_stat$slope[i]!=0){obs_stat$index[i]<-obs_stat$time[i]}
    
  }
  
  
  major.increases <- TSA_subsetting(obs_stat)
  
  
  major.increases <- lapply(major.increases, function(x) {
    if(sum(x$slope1)!=0){return(x)}
    else{x <- NULL}
    
  })
  
  major.increases <- major.increases[!sapply(major.increases, is.null)]
  
  # extract the steepest parts of the warming periods
  steepest.part <- lapply(major.increases, function(X) f.detect.steep(X, sdsteep))
  
  # select all warming periods that reach the max temp
  steepest.part <- lapply(steepest.part, function(X) f.select.by.maxtemp(X,max_temp))
  
  steepest.part <- steepest.part[!sapply(steepest.part, is.null)]
  
  if(length(steepest.part)==0){return(table)}
  
  # exclude all warming periods that start in positive temperatures
  steepest.part <- lapply(steepest.part, function(X) f.select.by.start.temp(X))
  
  steepest.part <- steepest.part[!sapply(steepest.part, is.null)]
  
  if(length(steepest.part)==0){return(table)}
  
  
  # calculate the mean slope for the 10 days after the steepest part of the warming period
  steepest.part <- lapply(steepest.part, function(X) f.calc.period.after.warming(X, dataset))
  
  
  # Create data table with all the information
  table <- lapply(steepest.part, function(X) f.start.and.stats(X,table))
  
  table <- do.call(rbind, table)
  
  
  # select all the RD dates
  table$RD <- F
  
  table <- do.call("rbind", as.list(
    by(table, table$Snow.Start, function(X){
      true.row <- which.max(X$Mean.Slope.after.warming)
      X$RD[true.row] <- T
      X
    })))
  
  row.names(table) <- c(1:length(table$loc_name))
  

  return(table)
  
}






# function to detect the steepest parts of a temperature increase period
f.detect.steep <- function(obs_stat, sdsteep){
  
  if(length(obs_stat$loc_name)==1){return(obs_stat)}
  
  std <- sdsteep*(sd(obs_stat$slope))
  
  obs_stat <- obs_stat[!obs_stat$slope<std,]
  
  return(obs_stat)
}





# function to check that the warming period reaches a certain temperature limit
f.select.by.maxtemp <- function(obs_stat,max_temp){
  
  if(tail(obs_stat$agg_max,n=1)>max_temp){return(obs_stat)}
  else{obs_stat <- NULL}
}


# function to exclude all warming periods that start above 0
f.select.by.start.temp <- function(obs_stat){
  
  if(head(obs_stat$agg_avg,n=1)<0){return(obs_stat)}
  else{obs_stat <- NULL}
  
}


# bring all important parameters toghether for one warming period
f.start.and.stats <- function(obs_stat, table){
  
  start.date <- head(obs_stat$time,n=1)-1
  
  start.date <- as.POSIXct(substring(start.date,1,10),origin="1970-01-01")
  
  table$Warming.Period.Start <- start.date #add the warming period start date
  
  table$RD <- NA
  
  table$Temp.Diff <- round(sum(obs_stat$slope),digits=3) #add the temperature difference within warming period
  
  table$Mean.Slope <- round(mean(obs_stat$slope),digits=3) #add mean slope over warming period
  
  table$Temp.To.Zero <- round(tail(obs_stat$agg_max, n=1),digits=3) #add the temperature difference to zero
  
  table$Mean.Slope.after.warming <- round(tail(obs_stat$Mean.Slope.after.warming, n=1), digits=3)
  
  return(table)
}


f.calc.period.after.warming <- function(X, dataset, table){
  
  
  warming.end <- tail(X$time, n=1)
  
  row.in.dataset <- which(dataset$time==warming.end)
  
  after.warming <- row.in.dataset + 10
  
  after.warming.date <- dataset$time[after.warming]
  
  X$Mean.Slope.after.warming <- mean(dataset$slope[dataset$time >= warming.end & dataset$time <= after.warming.date])
  
  return(X)
  
}


