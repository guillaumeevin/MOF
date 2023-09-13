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

###===============================###===============================###
#' disag.3D.to.1D
#'
#' This function applies the method of fragments to disaggregate daily simulations
#' of precipitation, as done in several papers, e.g.
#' \insertCite{breinlSpacetimeDisaggregationPrecipitation2019}{MethodOfFragments}.
#' @useDynLib MethodOfFragments, .registration = TRUE
#' @param YObs1H matrix of observed intensities at 1h: (nTobs*24) x nStation
#' @param YObs24H matrix of observed intensities at 24h: nTobs x nStation
#' @param YSim24H matrix of simulated intensities per 3-day period: nTsim x nStation
#' @param timeObs24H vector of time corresponding to YObs24H (used to obtain the corr. seasons)
#' @param timeSim24H vector of time corresponding to YSim24H (used to obtain the corr. seasons)
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#' @importFrom Rdpack reprompt
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#' @references
#' \insertAllCited{}
#' @export
#' @author Guillaume Evin
disag.3D.to.1D = function(YObs1H, # matrix of observed intensities at 1h: (nTobs*24) x nStation
  YObs24H, # matrix of observed intensities at 24h: nTobs x nStation
  YSim24H, # matrix of simulated intensities per 3-day period: nTsim x nStation
  timeObs24H, # vector of time corresponding to YObs24H
  timeSim24H, # vector of time corresponding to YSim24H
  prob.class = c (0.5, 0.75, 0.9, 0.99) # vector of probabilities indicating class of "similar" mean intensities
){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(YObs24H))
  nTsim = as.integer(nrow(YSim24H))

  ###### number of stations
  nStat = as.integer(ncol(YObs1H))

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
    q.mean.s = stats::quantile(mean.s, probs=prob.class, na.rm=T)
    if(any(q.mean.s==0)){
      q.mean.s = q.mean.s[q.mean.s!=0]
    }
    # observed class
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
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
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
    classSim[iSim.s] = class.s
  }

  ###### replace NA values by -9999 (can be processed in Fortran)
  YObs1H[is.na(YObs1H)] = -9999
  YObs24H[is.na(YObs24H)] = -9999

  ##### call Fortran function
  disag3Day.out = .Fortran("disagPrec_F",  PACKAGE="MethodOfFragments",
                           Yobs24=YObs24H, Yobs1=YObs1H, mObs=as.integer(seasObs24H), cObs=as.integer(classObs),
                           Ysim24=YSim24H, mSim=as.integer(seasSim24H), cSim=as.integer(classSim),
                           nTobs=nTobs, nStat=nStat, nTsim=nTsim, nLagScore=as.integer(1),
                           Ysim=matrix(0,nrow=nTsim*24,ncol=nStat),codeDisag=matrix(0,nrow=nTsim,ncol=nStat))
  return(list(Ysim=disag3Day.out$Ysim,codeDisag=disag3Day.out$codeDisag))
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
