###===============================###===============================###
### Guillaume Evin
### 16/01/2023, Grenoble
###  INRAE-UR ETNA
### guillaume.evin@inrae.fr
###
###  Provide a function to disaggregate simulated daily precipitation using oberved
### hourly precipitation.
###
### Wilks, D.S. (1998) "Multisite generalization of a daily stochastic
### precipitation generation model", J Hydrol, 210: 178-191
###
### Breinl, Korbinian, and Giuliano Di Baldassarre. 2019. Space-Time Disaggregation
### of Precipitation and Temperature across Different Climates and Spatial Scales
### Journal of Hydrology: Regional Studies
### https://doi.org/10.1016/j.ejrh.2018.12.002.

###===============================###===============================###
#' disagTempD2H
#'
#' This function applies the method of fragments to disaggregate daily simulations
#' of temperature, as done in several papers.
#' @useDynLib MethodOfFragments, .registration = TRUE
#' @param YObsXXH matrix of observed temperatures at a XXH time step: (nTobs*npdt) x nStation
#' @param YObs24H matrix of observed temperatures at a 24H time step: nTobs x nStation
#' @param YSim24H matrix of simulated temperaturess at a 24H time step: nTsim x nStation
#' @param timeObs24H vector of time corresponding to YObs24H (used to obtain the corr. seasons)
#' @param timeSim24H vector of time corresponding to YSim24H (used to obtain the corr. seasons)
#' @param npdt Number of time steps per day for YObsXXH (e.g. 24 for hourly observations)
#'
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp evalCpp
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#' @references
#' \insertAllCited{}
#' @export
#' @author Guillaume Evin
disagTempD2H = function(YObsXXH, # matrix of observed intensities at XXH: (nTobs*24) x nStation
                        YObs24H, # matrix of observed intensities at 24h: nTobs x nStation
                        YSim24H, # matrix of simulated intensities at a daily scale: nTsim x nStation
                        timeObs24H, # vector of time corresponding to YObs24H
                        timeSim24H, # vector of time corresponding to YSim24H
                        npdt=24){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(YObs24H))
  nTsim = as.integer(nrow(YSim24H))
  nTobsXXh = as.integer(nrow(YObsXXH))

  #### check if npdt matches the dimensions
  if(npdt!=(nTobsXXh/nTobs)){
    stop("npdt must correspond to the ratio between the number of lines in YObsXXH and YObs24H")
  }

  ###### number of stations
  nStat = as.integer(ncol(YObsXXH))

  ##### season: first agreement between seasons
  seasObs24H = month2season(as.numeric(format(timeObs24H,'%m')))
  seasSim24H = month2season(as.numeric(format(timeSim24H,'%m')))


  ##### call Fortran function
  disag.out = disagTempMOF(Npdt=as.integer(24),
                           mObs=seasObs24H,
                           mSim=seasSim24H,
                           YobsXX=YObsXXH,
                           Yobs24=YObs24H,
                           Ysim24=YSim24H)
  return(disag.out)
}

#==============================================================================
#' month2season
#'
#' transform vector of months to seasons
#'
#' @param vecMonth a vector of months given as integers 1:12
#'
#' @author Guillaume Evin
month2season = function(vecMonth){
  iSeason = c(1,1,2,2,2,3,3,3,4,4,4,1)
  return(iSeason[vecMonth])
}
