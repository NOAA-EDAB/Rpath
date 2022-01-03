#'Fishing Mortality Table
#'
#'Creates a table of fishing mortalities by species group and gear for an 
#'\code{Rsim.scenario} object.
#'
#'@family Rpath functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a data table of F values for each species/gear combination.
#'
#'@export
#' 
frate.table <- function(Rsim.scenario){
  #Need to define variables to eliminate check() note about no visible binding
  Group <- Gear <- Q <- NULL
  
  fish <- data.table(Group = Rsim.scenario$params$FishFrom,
                     Gear  = Rsim.scenario$params$FishThrough,
                     Q     = Rsim.scenario$params$FishQ)
  group <- unique(fish[Group > 0, Group])
  gear  <- unique(fish[Gear  > 0, Gear])
  
  group.name <- Rsim.scenario$params$spname[unique(fish[Group > 0, Group]) + 1]
  gear.name  <- Rsim.scenario$params$spname[unique(fish[Gear  > 0, Gear ]) + 1]
  
  fish.out <- c()
  for(i in 1:length(group)){
    fish.group <- fish[Group == group[i], ]
    fish.all.gear <- data.table(Group = group.name[i])
    for(j in 1:length(gear)){
      f.gear <- data.table(Group = group.name[i], V1 = fish.group[Gear == gear[j], sum(Q)])
      setnames(f.gear, 'V1', gear.name[j])
      fish.all.gear <- merge(fish.all.gear, f.gear, by = 'Group')
    }
    fish.out <- rbindlist(list(fish.out, fish.all.gear))
  }
  Total <- fish.out[, rowSums(.SD), .SDcols = gear.name]
  fish.out <- cbind(fish.out, Total)
  return(fish.out)
}

#'Adjust Fishing Mortality
#'
#'Modifies the fishing mortality value for a species by a particular gear.  
#'Parameters that can be adjusted using this function are: \emph{ForcedEffort},
#'\emph{ForcedFRate}, or \emph{ForcedCatch}. 
#'
#'@family Adjust functions
#'
#'@inheritParams rsim.run
#'@param parameter Parameters to be modified (See Description)
#'@param group Name of the group whose parameter is being changed.
#'@param sim.year Year of the simulation that should be modified.  Can be a range of years.
#'@param sim.month Month of the year that should be modified.  If set to 0, all months of
#'             the year are modified.
#'@param value New value for the parameter.
#'@return Returns an \code{Rsim.scenario} object with the new fishing parameter 
#'    values.
#'@export 
adjust.fishing <- function(Rsim.scenario, parameter, group = NA, sim.year = 1, 
                           sim.month = 0, value){
  #Check that parameter and group exist
  if(!parameter %in% c('ForcedEffort', 'ForcedFRate', 'ForcedCatch')){stop("Fishing parameter not found")}
  
  if(!all(group %in% Rsim.scenario$params$spname)){
     stop("Groups not found:",group[!(group %in% Rsim.scenario$params$spname)])
  }
  #if(!group %in% Rsim.scenario$params$spname){stop("Group not found")}
  
  #Create index in case the number of values is equal to length(sim.year) * length(sim.month)
  ivalue <- 0 
  
  #Loop over years if more than 1 sim.year provided
  for(iyear in seq_along(sim.year)){
    for(imonth in seq_along(sim.month)){
      ivalue <- ivalue + 1
      #look-up what rows correspond to the year
      #Regex used to account for year.month row names in Effort Matrix
      year.row <- which(gsub("\\..*", "", rownames(Rsim.scenario$fishing[[parameter]]))
                        == sim.year[iyear])
      
      #identify what rows correspond to the sim.months - 0 indicates the whole year
      if(sim.month[1] != 0){
        year.row <- year.row[1:12 %in% sim.month[imonth]]
        }
      
      #Apply the value to the correct row
      if(length(value) == 1){
        #If only 1 value is supplied for multiple years need to only point to that value
        Rsim.scenario$fishing[[parameter]][year.row, group] <- value
        }else if(sim.month[1] == 0){
          Rsim.scenario$fishing[[parameter]][year.row, group] <- value[iyear]
          }else if(length(value) == length(sim.month)){
            Rsim.scenario$fishing[[parameter]][year.row, group] <- value[imonth]
            }else {
              Rsim.scenario$fishing[[parameter]][year.row, group] <- value[ivalue]
            }
      }
    }

  return(Rsim.scenario)
}
  
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

#'Adjust Forcing Parameters
#'
#'Modifies the various forcing parameters of the rsim scenario object.
#'
#'@family Rpath functions
#'
#'@param Rsim.scenario Object generated by rsim.scenario.
#'@param parameter The Rsim.scenario forcing parameter to be modified. 
#'@param group The model group that the parameter change will affect.
#'@param sim.year Year of the simulation that should be modified.  Can be a range of years.
#'@param sim.month Month of the year that should be modified.  If set to 0, all months of
#'             the year are modified.
#'@param bymonth Boolean value that denotes whether to use sim.year/sim.month combo
#'               or just sim.month as a sequential vector starting at 1.
#'@param value The new value for the parameter.
#'@return Returns an Rsim.scenario object with the new parameter.
#'@export 
adjust.forcing <- function(Rsim.scenario, parameter, group, sim.year = 1, sim.month = 0, 
                           bymonth = F, value){
  #Check that parameter and group exist
  if(!parameter %in% c('ForcedPrey', 'ForcedMort', 'ForcedRecs', 'ForcedSearch', 'ForcedActresp',  
                       'ForcedMigrate', 'ForcedBio')){stop("Forcing parameter not found")}
  if(!all(group %in% Rsim.scenario$params$spname)){
    stop("Groups not found:",group[!(group %in% Rsim.scenario$params$spname)])
  }
  #if(!group %in% Rsim.scenario$params$spname){stop("Group not found")}
  
  #Create vector of values if only one supplied for multiple years
  if(length(value) == 1 & length(sim.year) > 1) value <- rep(value, length(sim.year))
  
  if(bymonth == F){
    #Forcing matrices are by month not year so need to convert year/month combo
    #Loop over years if more than 1 sim.year provided
    for(iyear in seq_along(sim.year)){
      #look-up what rows correspond to the year
      #Regex used to account for year.month row names in Forcing Matrices
      year.row <- which(gsub("\\..*", "", rownames(Rsim.scenario$forcing[[parameter]]))
                        == sim.year[iyear])
      
      #If more than 1 month is provided than need to identify which rows
      #correspond with the year month combo...also if month does not equal 0.
      if(length(sim.month) > 1){
        year.row <- year.row[1:12 %in% sim.month]
      }else{
        if(sim.month > 0) year.row <- year.row[1:12 == sim.month]
      }
      
      Rsim.scenario$forcing[[parameter]][year.row, group] <- value[iyear]
    }
  }else{
    Rsim.scenario$forcing[[parameter]][sim.month, group] <- value
  }
  
  
  return(Rsim.scenario)
}

# KYA - Some "adjustment" functions for Rsim.  Used "set" instead of "adjust" as
# that's more of a standard naming convention.  (In future make this a method
# of objects?)

#'Set Rsim.scenario parameters
#'
#'Modifies the various parameters of the \code{rsim.scenario} object. Parameters 
#'that can be adjusted using this function are: 
#'
#'@family Set functions
#'
#'@inheritParams rsim.run
#'@inheritParams rsim.fishing
#'@param start_state Rsim starting values object generated by \code{rsim.state}
#'@param forcing Rsim forcing matrix object generated by \code{rsim.forcing} 
#'@param fishing Rsim fishing matrix object generated by \code{rsim.fishing}
#'@param stanzas Rsim stanza parameters object generated by \code{rsim.stanzas}
#'
#'@return Returns an \code{Rsim.scenario} object with the new parameter.
#'
#'@export
#'  
set.rsim.scene<-function(Rsim.scenario,params=NULL,start_state=NULL,forcing=NULL,fishing=NULL,stanzas=NULL){
  rsim <- list()
  class(rsim) <- 'Rsim.scenario'
  attr(rsim, 'eco.name') <- attr(Rsim.scenario, 'eco.name')
  # can add type checks later
    if (!is.null(params))     {rsim$params      <- params      } else {rsim$params      <- Rsim.scenario$params      }
    if (!is.null(start_state)){rsim$start_state <- start_state } else {rsim$start_state <- Rsim.scenario$start_state }
    if (!is.null(forcing))    {rsim$forcing     <- forcing     } else {rsim$forcing     <- Rsim.scenario$forcing     }
    if (!is.null(fishing))    {rsim$fishing     <- fishing     } else {rsim$fishing     <- Rsim.scenario$fishing     }
    if (!is.null(stanzas))    {rsim$stanzas     <- stanzas     } else {rsim$stanzas     <- Rsim.scenario$stanzas     }
  return(rsim)
}

#'Retrieve Rsim scenario parameters
#'
#'Helper function that will retrieve parameters that were used in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.params} object.
#'
#'@export
#'    
get.rsim.params<-function(Rsim.scenario){
  return(Rsim.scenario$params)
}

#'Retrieve Rsim starting state values
#'
#'Helper function that will retrieve the starting state values that were used 
#'in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.state} object.
#'
#'@export
#'   
get.rsim.start_state<-function(Rsim.scenario){
  return(Rsim.scenario$start_state)
}

#'Retrieve Rsim forcing parameters
#'
#'Helper function that will retrieve the forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.forcing} object.
#'
#'@export
#'   
get.rsim.forcing<-function(Rsim.scenario){
  return(Rsim.scenario$forcing)
}

#'Retrieve Rsim fishing forcing parameters
#'
#'Helper function that will retrieve the fishing forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.fishing} object.
#'
#'@export
#'  
get.rsim.fishing<-function(Rsim.scenario){
  return(Rsim.scenario$fishing)
}

#'Retrieve Rsim stanza parameters
#'
#'Helper function that will retrieve the stanza parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.stanzas} object.
#'
#'@export
#' 
get.rsim.stanzas<-function(Rsim.scenario){
  return(Rsim.scenario$stanzas)
}







