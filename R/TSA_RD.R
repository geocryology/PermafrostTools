# =============================================================================
#'
#' @title Detect RD date using the theory and function by Schmid et al. (2012)
#'
#' @description Give one data.frame of a location or a list of data.frames of multiple 
#'              locations and the function returns the Ripening Date for the snow periods.
#'
#' @details It is also possible to adjust the parameters for the TSA_snow_cover function in this function, 
#'          as the warming periods are calculated for each snow cover period. 
#'          Therefore, the TSA_snow_cover is run as well in this function.
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
#' @param v 	  Threshold zero curtain (if dailymax & dailymin is within +/- v, this day is a zero curtain day)
#'              Default: 0.25
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
#' @return Dataframe with fields: loc_name; Snow.Start; Snow.End; MDr
#'           
#' 
#' @export 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#' RD_table <- TSA_RD(BC_SO01_01, v= 0.25, m1 = 3, m2 = 0.5, v1 = 0.1, v2 = 0.3, lenghtsnow = 5, MDr.sd = 0.3, snow.period=NULL)
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_RD <- function(obs_stat, v= 0.25, v1 = 0.1, v2 = 0.3, lenghtsnow = 5, MDr.sd = 0.3, snow.period=NULL){
  if(class(obs_stat)=="data.frame"){obs_stat <- list(obs_stat)}
  
  #create an error dataframe
  error.frame <- as.data.frame(matrix(ncol = 7, nrow = 1))
  names(error.frame) <- c("loc_name","Snow.Start","Snow.End","RD.Date")
  
  #run the f.RD.marcol1 function over the whole dataset
  RD.table.list <- lapply(obs_stat, FUN = function(X) tryCatch(f.RD.marcol1(X, v, v1, v2, lenghtsnow, MDr.sd, snow.period), error=function(e) error.frame))
  
  #create a dataframe from the output list
  RD.table <- do.call(rbind, RD.table.list)
  
  #delete the Row.names column
  RD.table$Row.names <- NULL
  
  return(RD.table)
}




#function to calculate the RD for one location
f.RD.marcol1 <- function(obs_stat, v, v1, v2, lenghtsnow, MDr.sd, snow.period){
  
  "&" <- function(...) UseMethod("&")
  "&.default" <- .Primitive("&")
  "&.character" <- function(...) paste(...,sep="")
  
  if(is.null(snow.period)){snow.period <- TSA_snow_cover(obs_stat, v1, v2, lenghtsnow, MDr.sd)}
  
  #create data table
  table <- TSA_row(obs_stat$loc_name[1])
  
  table$RD.Date <- NA
  
  
  if(is.na(snow.period$Snow.Start[1])){return(table)}

  
  
  ### freezing index (FI)	
  FI<-sum(obs_stat$agg_avg[which(obs_stat$agg_avg<0)])
  
  ## snow periods
  obs_stat$snow.ind <- 0
  
  for(i in 1:length(snow.period$loc_name)){
    
    obs_stat$snow.ind[snow.period$Snow.Start[i] <= obs_stat$time & obs_stat$time <= snow.period$Snow.End[i]] <- 1
    
  }
  
  obs_stat$index <- obs_stat$snow.ind*obs_stat$agg_sd  
  
  ### zero curtain index
  obs_stat$zc.ind <- 0
  obs_stat$zc.ind[obs_stat$agg_min <= v&obs_stat$agg_min>-v&obs_stat$agg_max <= v&obs_stat$agg_max>-v] <- 1
  
  
  ## freezing index
  obs_stat$freez.ind <- 0
  obs_stat$freez.ind[obs_stat$agg_avg<(-v)]<-1
  
  
  # make subset of the snow periods
  for(ii in 1:length(obs_stat$index)){
    
    if(obs_stat$index[ii]==1){obs_stat$index[ii]<-obs_stat$time[ii]}
    
  }
  
  snow <- TSA_subsetting(obs_stat)
  
  #calculate RD for each snow period
  table1 <- lapply(snow, function(X) tryCatch(f.RDcalculate(X,table, FI),error=function(e) table))
  
  RD.table <- do.call(rbind, table1)
  
  RD.table <- RD.table[complete.cases(RD.table), ]
  
  #make sure all times are in right format
  RD.table$RD.Date <- as.POSIXct(RD.table$RD.Date, origin="1970-01-01")
  RD.table$Snow.Start <- as.POSIXct(RD.table$Snow.Start, origin="1970-01-01")
  RD.table$Snow.End <- as.POSIXct(RD.table$Snow.End, origin="1970-01-01")
  
  
  return(RD.table)
}


#function to calculate the RD for one snow period
f.RDcalculate <- function(obs_stat,table, FI){
  
  table$Snow.Start <- head(obs_stat$time,n=1)
  
  table$Snow.End <- tail(obs_stat$time,n=1)
  
  
  t.FI<-	-50
  
  if(sum(obs_stat$freez.ind)!=0){
    j <- 0
    k <- 0
    for(ii in 1:(length(obs_stat$freez.ind)-1)){
      if(obs_stat$freez.ind[ii] == 1 & obs_stat$freez.ind[ii + 1] != 0){j <- j + 1}
      if(obs_stat$freez.ind[ii] == 1 & obs_stat$freez.ind[ii + 1] == 0 & k > j){j <- 0}
      if(obs_stat$freez.ind[ii] == 1 & obs_stat$freez.ind[ii + 1] == 0 & k <= j){ k <- j; 
      freez.end <- ii; j <- 0}
    }   
  
  
    obs_stat$zc.ind<-obs_stat$zc.ind*obs_stat$snow.ind
    obs_stat$zc.ind[1:freez.end]<-0
  
  }else{
    obs_stat$zc.ind<-obs_stat$zc.ind*obs_stat$snow.ind
  }

  if(sum(obs_stat$zc.ind)==0){return(table)}

  # subsetting the zero curtain
  #onset of melting period
  j <- 0
  k <- 0
  for(a in 1:(length(obs_stat$zc.ind)-1)){
    if(obs_stat$zc.ind[a] == 1){j <- j + 1}
    if(obs_stat$zc.ind[a+1] == 0 & k <= j){ k <- j; ind.end <- a; j <- 0} 
  }      
  ind.start <-  ind.end -k + 1

  RD  <- 	obs_stat$time[ind.start]



  # Fill "table" data table with all the information
  #add snow end date

  table$RD.Date <- as.POSIXct(substring(RD,1,10)) #add the warming period start date


  if(sum(obs_stat$zc.ind)==0){table$RD.Date <- NA}   #no date if no zero curtain
  if(FI>t.FI){table$RD.Date <- NA}                   #date only if FI is big enough
  if(table$RD.Date>table$Snow.End){table$RD.Date <- NA}                     #date only if basal ripening is before melt out date

  return(table)
}
