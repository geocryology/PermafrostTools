# =============================================================================
#'
#' @title Plot of the time series of a location
#'
#' @description Give one data.frame of a location or a list of data.frames of multiple 
#'              locations and the function returns a plot containing the minimum, maximum, mean daily GST
#'              as well as the snow cover, RD, warming periods and zero curtains. The output is a
#'              pdf-file, saved in a own defined output directory
#'
#' @details 
#'
#' @param obs_stat one dataframe of one location or a list of dataframes of multiple locations.
#'                 A dataframe needs at least the following columns:
#'                 Location name as loc_name;
#'                 Mean daily GST as agg_avg;
#'                 Minimum daily GST as agg_min;
#'                 Maximum daily GST as agg_max;
#'                 Daily sd as agg_sd;
#'                 Date as time (in POSIXct and the following format: YYYY-MM-DD)
#' 
#' @param out.path path to directory where the plot should be saved. e.g: "~/Desktop/"
#' 
#' @param v1    max. daily standard deviation indicating snow for POSITIVE Ground Surface Temperature. 
#'              Default is 0.1.
#'              
#' @param v2    max. daily standard deviation indicating snow for NEGATIVE Ground Surface Temperature. 
#'              Default is 0.3.
#'              
#' @param lengthsnow    number of days that the snow cover should at least have to be selected as a valid snow cover.
#'                     Default is 5 days
#'
#' @param MDr.sd   threshold for the mean daily standard deviation for a snow period. 
#'
#' @param v 	  threshold zero curtain (if dailymax & dailymin is within +/- v, this day is a zero curtain day)
#'              Default: 0.25
#'
#' @param tperc    percentage of the 10-percent-quantile of the winter temperatures that the 
#'                     warming period should at least reach. Default is 40 percent.
#' 
#' @param sdoutlier   The number implies how many multiples of the standard deviation of the slope are taken to define the threshold.  
#'             The threshold is used to detect slope outliers.
#' 
#' @param sdsteep    The number implies what fraction of the standard deviation of the slope 
#'               is taken to select the steepest part of a warming period.
#'               
#' @param temp    the temperature boundary within the mean agg_avg of a zero curtain period has to be. 
#'                Default: 0.2 degrees Celsius     
#' 
#' @param slopesd     Fraction of the standard deviation of the slope used as threshold for the slope. Default: 0.1
#'             
#' 
#' @param tempsd     Fraction of the mean daily standard deviation used as threshold for the daily sd. Default: 0.01
#'               
#' @param lengthzc     Length that a zero curtain period should at least have. Default: 2 days.
#'
#'
#' @return Pdf-file saved in a own defined directory. Plot contains:
#'         min-, max-, mean- daily temperatures as well as snow cover, RD, warming periods, zero curtains
#'           
#' 
#' @export 
#' @examples
#' con <- dbpf_con()
#' BC_SO01_01 <- TSA_data_import(con, "BC16-SO01_01")
#' TSA_snow_cover(BC_SO01_01, "~/Desktop/", v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, v= 0.25, tperc = 40, sdoutlier = 3, sdsteep = 0.75, temp= 0.2, slopesd= 0.1, tempsd= 0.01, lengthzc= 2)
#'
#' @author Thomas Knecht <t.knecht@@hotmail.com>
# =============================================================================

TSA_plot <- function(obs_stat,out.path, v1 = 0.1, v2 = 0.3, lengthsnow = 5, MDr.sd = 0.3, v= 0.25, tperc = 40, sdoutlier = 3, sdsteep = 0.75, temp= 0.2, slopesd= 0.1, tempsd= 0.01, lengthzc= 2){
  if(class(obs_stat)=="data.frame"){obs_stat <- list(obs_stat)}
  
  invisible(lapply(obs_stat, FUN = function(X) tryCatch(f.single.plot(X, out.path, v1, v2, lengthsnow, MDr.sd, v, tperc, sdoutlier, sdsteep, temp, slopesd, tempsd, lengthzc), error=function(e) NA)))
  
}


#plot function for one location
f.single.plot <- function(obs_stat, out.path, v1, v2, lengthsnow, MDr.sd, v, tperc, sdoutlier, sdsteep, temp, slopesd, tempsd, lengthzc){
 
  graphics.off()
  
  #define line boundaries
  Zeroheight <- (mean(obs_stat$agg_avg[obs_stat$agg_avg>0]))/2
  
  snowheight <- (max(obs_stat$agg_avg))/2
  
  yheight1 <- min(obs_stat$agg_avg)
  yheight2 <- max(obs_stat$agg_avg)
  
  #calculate all the dates
  snow.period <- TSA_snow_cover(obs_stat, v1, v2, lengthsnow, MDr.sd)
  
  RD.marcol <- TSA_RD(obs_stat, v, snow.period=snow.period)
  
  warming.periods <- TSA_warming_periods(obs_stat, tperc, sdoutlier, sdsteep, snow.period=snow.period)
  
  zero_curtains <- TSA_zero_curtain(obs_stat, temp, slopesd, tempsd, lengthzc, snow.period=snow.period)
  
  title <- as.character(obs_stat$loc_name[1])
  
  #plot the dates if they exist
    if(is.na(warming.periods$Warming.Period.Start[1])){
      warmingplot <- geom_segment()
    }else{
      warmingplot <- geom_segment(aes(x=warming.periods$Warming.Period.Start,xend=warming.periods$Warming.Period.Start,y=yheight1,yend=yheight2, color="Warming_Start"),size=0.5)
    }

    if(is.na(zero_curtains$Zero.Curtain.Start[1])){
      zeroplot <- geom_segment()
    }else{
      zeroplot <- geom_segment(aes(x=zero_curtains$Zero.Curtain.Start,xend=zero_curtains$Zero.Curtain.End,y=Zeroheight,yend=Zeroheight, color="Zero_Curtains"),size=5)
    }        
  
    if(is.na(snow.period$Snow.Start[1])){
      snowplot <- geom_segment()
    }else{
      snowplot <- geom_segment(aes(x=snow.period$Snow.Start,xend=snow.period$Snow.End,y=snowheight,yend=snowheight, color="Snow_Cover"),size=5) 
    } 
    
  
    if(is.na(RD.marcol$RD.Date[1])){
      rdplot <- geom_segment()
    }else{
      rdplot <- geom_segment(aes(x=RD.marcol$RD.Date,xend=RD.marcol$RD.Date,y=yheight1,yend=yheight2, color="RD_MarcOl"))
    }
  
  
  myplot <- ggplot() + geom_line(data=obs_stat, aes(time, agg_avg,color="Mean_Temp")) +
            geom_line(data=obs_stat, aes(time, agg_max,color="Max_Temp")) +
            geom_line(data=obs_stat, aes(time, agg_min,color="Min_Temp")) +
            snowplot +        
            zeroplot +
            warmingplot +
            rdplot +
            ylab(expression("Temperature"~(degree*C)*"")) +
            xlab("Date") +
            ggtitle(title) +
            scale_colour_manual("", values = c(Mean_Temp="black", Max_Temp="blue", Min_Temp="grey", Snow_Cover="green", Zero_Curtains="red", Warming_Start="black", RD_MarcOl="brown"))+ 
            guides(colour = guide_legend(override.aes = list(size=0.5))) +
            theme_bw()
  
  #export as pdf file  
  mypath <- file.path(out.path,paste(obs_stat$loc_name[1],".pdf",sep="") )
  pdf(file=mypath,width=7,height=5)
  print(myplot)
  dev.off()  

}






