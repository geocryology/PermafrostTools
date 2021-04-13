# =============================================================================
#'
#' @title Return data frame or a list of data frames with the daily aggregated data for a location
#'        or several locations.
#'
#' @description Give a data.frame, a list or a character of location-names and the function returns a data frame or a list 
#'              of data frames with daily aggregated data for each location.
#'
#' @details The default settings are set as to import the whole dataset of a location. By adjusting the time_b and 
#'          time_e parameters also only single years or winter periods can be imported.
#'
#' @param con Connection to Permafrost data base
#'      
#' @param inventory   A data frame with a column "name", a list or a character containing all the location names from which
#'                    the data should be imported    
#' 
#' @param time_b     Begin date for the import. Format: "1950-01-01 00:00:00+00"
#'             
#' @param time_e     End date for the import. Format: "2050-01-01 00:00:00+00"
#'
#'
#' @return List of data fames with the following fields: 
#'         Location name as loc_name;
#'         Observation height/depth as height [m];
#'         Mean daily GST as agg_avg [C];
#'         Number of daily observations as agg_cnt;
#'         Minimum daily GST as agg_min [C];
#'         Maximum daily GST as agg_max [C];
#'         Mean daily standard deviation as agg_sd;
#'         Date as time (in POSIXct and the following format: YYYY-MM-DD);
#'         
#'         
#' 
#' 
#' 
#' @export 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_data_import <- function(con, inventory,time_b="1950-01-01 00:00:00+00",time_e="2050-01-01 00:00:00+00"){
  
  options(stringsAsFactors = FALSE)
  
  #retrieve location names and store it as vector
  if(class(inventory)=="character"){loc_names <- inventory}
  else if(class(inventory)=="data.frame"){loc_names <- inventory$name}
  else if(class(inventory)=="list"){loc_names <- inventory}
  else if(class(inventory)=="vector"){loc_names <- inventory}
  
  
  #lapply the f.aggregate function over the loc_name list
  loc_data.list <- lapply(loc_names, FUN = function(X) f.aggregate(X,time_b,time_e,con))
  
  return(loc_data.list)
  
}


### Function to load and aggregate the GST data from the database  
f.aggregate <- function(location_name,time_b,time_e,con){ 
  
  # make averaging period [s]
  period = 3600 * 24
  verbose = FALSE
  unit_of_measurement = "C"
  
  
  # make query for aggregation of the datalogger from database (avg, count, min, max, sd)
  q <- paste0("SELECT locations.name AS loc_name, ",
              "observations.height_min_metres AS height, ",
              "AVG(observations.numeric_value) as agg_avg, ",
              "COUNT(observations.numeric_value) as agg_cnt, ",
              "MIN(observations.numeric_value) as agg_min, ",
              "MAX(observations.numeric_value) as agg_max, ",
              "STDDEV(observations.numeric_value) as agg_sd, ",
              "TO_TIMESTAMP(FLOOR(EXTRACT('epoch' FROM observations.corrected_utc_time AT TIME ZONE 'UTC') / ",
              period, ") * ", period,
              ") AT TIME ZONE 'UTC' AS time ",	            
              "FROM observations INNER JOIN ",
              "locations ON observations.location = locations.coordinates ",
              "WHERE observations.corrected_utc_time BETWEEN ",
              "'", time_b, "' AND '", time_e, "' AND ",
              "locations.name = ANY('{", paste(location_name, collapse=", ") ,"}'::text[]) ",
              "AND observations.unit_of_measure = '", unit_of_measurement, "' ",
              "GROUP BY observations.height_min_metres, locations.name, time ",
              "ORDER BY loc_name ASC, height DESC;")
  
  #query
  if (verbose == TRUE) {
    print("=== Query sent:")
    print(q)
  }
  obs_stat <- as.data.frame(dbGetQuery(con, q))
  
  
  obs_stat$time <-  as.POSIXct(substring(obs_stat$time,1,10))    
  
  #filter days with less than 10 values
  obs_stat <- subset(obs_stat, obs_stat$agg_cnt > 10)
  
  
  return(obs_stat)
  
}
