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
### Breinl, Korbinian, and Giuliano Di Baldassarre. 2019. “Space-Time Disaggregation
### of Precipitation and Temperature across Different Climates and Spatial Scales.”
### Journal of Hydrology: Regional Studies 21 (February): 126–46.
### https://doi.org/10.1016/j.ejrh.2018.12.002.

###===============================###===============================###

###===============================###===============================###
#' disag.3D.to.1D
#' @useDynLib MethodOfFragments, .registration = TRUE
#' @noRd
#' @param Yobs1H matrix of observed intensities at 1h: (nTobs*24) x nStation
#' @param YObs24H matrix of observed intensities at 24h: nTobs x nStation
#' @param mObs24H vector of season corresponding to Yobs24H
#' @param YSim24H matrix of simulated intensities per 3-day period: nTsim x nStation
#' @param mSim24H vector of season corresponding to the period simulated
#' @param prob.class vector of probabilities indicating class of "similar" mean intensities
#'
#' @return \item{list}{Ysim matrix of disagregated daily precipitation, codeDisag matrix of disagregation codes}
#'
#' @author Guillaume Evin
disag.3D.to.1D = function(Yobs1H, # matrix of observed intensities at 1h: (nTobs*24) x nStation
  YObs24H, # matrix of observed intensities at 24h: nTobs x nStation
  YSim24H, # matrix of simulated intensities per 3-day period: nTsim x nStation
  timeObs24H, # vector of time corresponding to YObs24H
  timeSim24H, # vector of time corresponding to YSim24H
  prob.class =c (0.5, 0.75, 0.9, 0.99) # vector of probabilities indicating class of "similar" mean intensities
){
  ###### number of 3-day periods simulated
  nTobs = as.integer(nrow(Yobs24H))
  nTsim = as.integer(nrow(YSim24H))

  ###### number of stations
  nStat = as.integer(ncol(Yobs1H))

  ##### season: first agreement between seaons
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
    Yobs.s = Yobs24H[iObs.s,]
    mean.s = apply(Yobs.s,1,mean,na.rm=T)
    # 4 breaks by default: small, moderate, high, extremes precipitation
    q.mean.s = quantile(mean.s, probs=prob.class)
    if(any(q.mean.s==0)){
      q.mean.s = q.mean.s[q.mean.s!=0]
    }
    # observed class
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
    classObs[iObs.s] = class.s
    # simulated class
    iSim.s = seasSim24H==i.s
    Ysim.s = YSim24H[iSim.s,]
    mean.s = apply(Ysim.s,1,mean,na.rm=T)
    class.s = cut(mean.s, breaks=c(0,q.mean.s,max(mean.s)), labels = FALSE, include.lowest = T)
    classSim[iSim.s] = class.s
  }

  ###### replace NA values by -9999 (can be processed in Fortran)
  Yobs[is.na(Yobs1H)] = -9999
  Yobs24H[is.na(Yobs24H)] = -9999

  ##### call Fortran function
  disag3Day.out = .Fortran("disagPrec_F",  PACKAGE="MethodOfFragments",
                           Yobs=Yobs1H, Y3obs=Yobs24H, mObs=as.integer(seasObs24H), cObs=as.integer(classObs),
                           Y3sim=YSim24H, mSim=as.integer(seasSim24H), cSim=as.integer(classSim),
                           nTobs=nTobs, nStat=nStat, nTsim=nTsim, nLagScore=as.integer(1),
                           Ysim=matrix(0,nrow=nTsim*3,ncol=nStat),codeDisag=matrix(0,nrow=nTsim,ncol=nStat))
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
