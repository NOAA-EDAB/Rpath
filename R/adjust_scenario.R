#'Adjust Rsim.scenario parameters
#'
#'Modifies the various parameters of the \code{rsim.scenario} object. Parameters 
#'that can be adjusted using this function are:
#'
#'@family Adjust functions
#'
#'@inheritParams adjust.fishing
#'@param group The model group that the parameter change will affect.  Note that 
#'       a value of \emph{'all'} will affect all groups associate with the groupto
#'       variable.
#'@param groupto The corresponding group who's parameter is affecting the group
#'       variable.
#'@return Returns an \code{Rsim.scenario} object with the new parameter.
#'@export 
adjust.scenario <- function(Rsim.scenario, parameter, group, groupto = NA, value){
  #Lookup group numbers
  if(group == 'all'){
    groupnum <- 0:Rsim.scenario$params$NUM_GROUPS
  } else {
    groupnum <- Rsim.scenario$params$spnum[which(Rsim.scenario$params$spname 
                                                 == group)]
  }
  if(!is.na(groupto)){
    groupnumto <- Rsim.scenario$params$spnum[which(Rsim.scenario$params$spname 
                                                   == groupto)]
  }
  
  #Lookup parameter number
  param.num <- which(names(Rsim.scenario$params) == parameter)
  
  #Modify parameter
  if(parameter %in% c('B_BaseRef', 'MzeroMort', 'UnassimRespFrac', 'ActiveRespFrac',
                      'FtimeAdj', 'FtimeQBOpt', 'PBopt', 
                      'NoIntegrate', 'HandleSelf', 'ScrambleSelf')){
    Rsim.scenario$params[[param.num]][groupnum + 1] <- value
  }
  
  if(parameter %in% c('QQ', 'DD', 'VV', 'HandleSwitch', 'PredPredWeight', 
                      'PreyPreyWeight')){
    linknum <- which(Rsim.scenario$params$PreyFrom %in% groupnum &
                       Rsim.scenario$params$PreyTo == groupnumto)  
    Rsim.scenario$params[[param.num]][linknum] <- value
  }
  return(Rsim.scenario)
}