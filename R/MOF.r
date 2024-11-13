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

###===============================###===============================###
#' disagPrecD2H
#'
#' This function applies the method of fragments to disaggregate daily simulations
#' of precipitation, as done in several papers, e.g.
#' \insertCite{breinlSpacetimeDisaggregationPrecipitation2019}{MethodOfFragments}.
#' @useDynLib MethodOfFragments, .registration = TRUE
#' @param YObsXXH matrix of observed temperatures at a XXH time step: (nTobs*npdt) x nStation
#' @param YObs24H matrix of observed temperatures at a 24H time step: nTobs x nStation
#' @param YSim24H matrix of simulated temperaturess at a 24H time step: nTsim x nStation
#' @param timeObs24H vector of time corresponding to YObs24H (used to obtain the corr. seasons)
#' @param timeSim24H vector of time corresponding to YSim24H (used to obtain the corr. seasons)
#' @param npdt Number of time steps per day for YObsXXH (e.g. 24 for hourly observations)
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#' @param nlagscore integer indicating the number of preceeding days used to calcule the scores, can be 0, 1 or 2
#'
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp evalCpp
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#' @references
#' \insertAllCited{}
#' @export
#' @author Guillaume Evin
disagPrecD2H = function(YObsXXH, # matrix of observed intensities at XXH: (nTobs*24) x nStation
                        YObs24H, # matrix of observed intensities at 24h: nTobs x nStation
                        YSim24H, # matrix of simulated intensities at a daily scale: nTsim x nStation
                        timeObs24H, # vector of time corresponding to YObs24H
                        timeSim24H, # vector of time corresponding to YSim24H
                        npdt=24,
                        prob.class = c (0.5, 0.75, 0.9, 0.99), # vector of probabilities indicating class of "similar" mean intensities
                        nlagscore=1 # integer indicating the number of preceeding days used to calcule the scores
){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(YObs24H))
  nTsim = as.integer(nrow(YSim24H))
  nTobsXXh = as.integer(nrow(YObsXXH))

  #### check if npdt matches the dimensions
  if(npdt!=(nTobsXXh/nTobs)){
    stop("npdt must correspond to the ratio between the number of lines in YObsXXH and YObs24H")
  }

  #### check if npdt matches the dimensions
  if(!nlagscore%in%c(0,1,2)){
    stop("nlagscore must be equal to 0, 1 or 2")
  }

  ###### number of stations
  nStat = as.integer(ncol(YObsXXH))

  ##### season: first agreement between seasons
  seasObs24H = month2season(as.numeric(format(timeObs24H,'%m')))
  seasSim24H = month2season(as.numeric(format(timeSim24H,'%m')))


  ##### class: classification of precipitation events according
  # to the mean precipitation over all stations. 4 classes
  classObs = vector(length = nTobs)
  classSim =  vector(length = nTsim)
  # for each season
  for(i.s in 1:4){
    # mean obs
    iObs.s = seasObs24H==i.s
    Yobs.s = YObs24H[iObs.s,]
    if(nStat==1){
      mean.s = Yobs.s
    }else{
      mean.s = apply(Yobs.s,1,mean,na.rm=T)
    }
    # 4 breaks by default: small, moderate, high, extremes precipitation
    # they are obtained from the observations
    q.mean.s = stats::quantile(mean.s, probs=prob.class, na.rm=T)
    if(any(q.mean.s==0)){
      q.mean.s = q.mean.s[q.mean.s!=0]
    }
    breaks = c(0,q.mean.s,max(Yobs.s,na.rm=T))

    # observed class
    class.s = cut(mean.s, breaks=breaks, labels = FALSE, include.lowest = T)
    class.s[is.na(class.s)] = -9999
    classObs[iObs.s] = class.s

    # simulated class
    iSim.s = seasSim24H==i.s
    Ysim.s = YSim24H[iSim.s,]
    if(nStat==1){
      mean.s = Ysim.s
    }else{
      mean.s = apply(Ysim.s,1,mean,na.rm=T)
    }
    class.s = cut(mean.s, breaks=breaks, labels = FALSE, include.lowest = T)
    classSim[iSim.s] = class.s
  }

  ###### replace NA values by -9999 (can be processed in Fortran)
  YObsXXH[is.na(YObsXXH)] = -9999
  YObs24H[is.na(YObs24H)] = -9999

  ##### call C++ function
  disag.out = disagPrecMOF(Npdt=as.integer(24),
                           mObs=seasObs24H,
                           mSim=seasSim24H,
                           cObs=as.integer(classObs),
                           cSim=as.integer(classSim),
                           YobsXX=YObsXXH,
                           Yobs24=YObs24H,
                           Ysim24=YSim24H,
                           nLagScore=as.integer(nlagscore))

  return(list(Ysim=disag.out$Ysim,codeDisag=disag.out$codeDisag))
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
