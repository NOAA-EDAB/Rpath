#'Fishing Mortality Table
#'
#'Creates a table of fishing mortalities by species group and gear for an 
#'\code{rsim.scenario()} object.
#'
#'@family Rpath functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a data table of F values for each species/gear combination.
#'\item{Group}{Group/node}
#'\item{Fishery/gearType}{Gear type for the fishery} ...
#'\item{Total}{Sum over all gear types}
#'
#'
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' rates <- frate.table(Rsim.scenario)
#' # display the head of the data frame
#' head(rates)
#'
#'@export



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
#'@param parameter Parameters to be modified (\code{ForcedEffort},\code{ForcedCatch},\code{ForcedFRate})
#'@param group Name of the group whose parameter is being changed. Valid values are found in the `Group` field of the object
#' created from running \code{rpath()}
#'@param sim.year Year of the simulation that should be modified.  Can be a range of years.
#'@param sim.month Month of the year that should be modified.  If set to 0, all months of
#'             the year are modified.
#'@param value New value for the parameter.
#'
#'@return Returns an \code{Rsim.scenario()} object with the new fishing parameter 
#'    values.
#'    
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' # Change value of forcedFRate for Squids in years 3 through 5 to the value of 2 (for all months)
#' Rsim.scenario.adjusted.fishing <- adjust.fishing(Rsim.scenario,parameter="ForcedFRate",group="Squids",sim.year=3:5,value = 2)
#'    
#'    
#'    
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
#'Modifies the various parameters of the \code{rsim.scenario()} object. Parameters that can be adjusted using this function are: 
#'\emph{B_BaseRef}, \emph{MzeroMort},\emph{UnassimRespFrac}, \emph{ActiveRespFrac}, \emph{FtimeAdj},
#'\emph{FtimeQBOpt}, \emph{PBopt}, \emph{NoIntegrate},\emph{HandleSelf}, \emph{ScrambleSelf}, \emph{QQ},
#' \emph{DD}, \emph{VV}, \emph{HandleSwitch}, \emph{PredPredWeight}, \emph{PreyPreyWeight}
#'
#'@family Adjust functions
#'
#'@inheritParams adjust.fishing
#'
#'@param parameter Parameters to be modified (Choose from: \code{B_BaseRef, MzeroMort, 
#' UnassimRespFrac, ActiveRespFrac, FtimeAdj, FtimeQBOpt, PBopt, NoIntegrate,
#' HandleSelf, ScrambleSelf, QQ, DD, VV, HandleSwitch, PredPredWeight, PreyPreyWeight})
#'@param group The model group that the parameter change will affect.  Note that 
#'       a value of \emph{'all'} will affect all groups associated with the `groupto`
#'       variable. Valid values are found in the `Group` field of the object created
#'      from running \code{rpath()}
#'@param groupto The corresponding group who's parameter is affecting the group
#'       variable. Required for parameters \code{QQ}, \code{DD}, \code{VV}, \code{HandleSwitch},\code{PredPredWeight},
#'       \code{PreyPreyWeight}
#'@return Returns an \code{rsim.scenario()} object with the new parameter.
#'
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' # Adjust the PBopt parameter for cod. Set to value = 2
#' Rsim.scenario.adjusted <- adjust.scenario(Rsim.scenario, parameter="PBopt",group = "cod", groupto = "all", value = 2)   
#'
#'
#'@export 

adjust.scenario <- function(Rsim.scenario, parameter, group, groupto = NA, value){
  #Lookup group numbers
  if(group[1] == 'all'){
    groupnum <- 0:Rsim.scenario$params$NUM_GROUPS
  } else {
    groupnum <- Rsim.scenario$params$spnum[which(Rsim.scenario$params$spname 
                                                 %in% group)]
  }
  if(!is.na(groupto)){
    groupnumto <- Rsim.scenario$params$spnum[which(Rsim.scenario$params$spname 
                                                 %in% groupto)]
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
#'@family Adjust functions
#'
#'@inheritParams adjust.fishing
#'
#'@param bymonth Boolean value that denotes whether to use sim.year/sim.month combo
#'               or just sim.month as a sequential vector starting at 1.
#'               
#'@return Returns an Rsim.scenario object with the new parameter.
#'
#'
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' # Adjust the ForcedPrey parameter for cod in year 1 for all months. Change the value to 10
#' Rsim.scenario.adjusted <- adjust.forcing(Rsim.scenario, parameter="ForcedPrey",group = "cod", sim.year = 1, sim.month=0,value=10)   
#' head(Rsim.scenario.adjusted$forcing$ForcedPrey)
#'
#'
#'
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
  
  if(bymonth){
    Rsim.scenario$forcing[[parameter]][sim.month, group] <- value
  }else {
    
  #Create index in case the number of values is equal to length(sim.year) * length(sim.month)
  ivalue <- 0 
  
  #Loop over years if more than 1 sim.year provided
  for(iyear in seq_along(sim.year)){
    for(imonth in seq_along(sim.month)){
      ivalue <- ivalue + 1
      #look-up what rows correspond to the year
      #Regex used to account for year.month row names in Effort Matrix
      year.row <- which(gsub("\\..*", "", rownames(Rsim.scenario$forcing[[parameter]]))
                        == sim.year[iyear])
      
      #identify what rows correspond to the sim.months - 0 indicates the whole year
      if(sim.month[1] != 0){
        year.row <- year.row[1:12 %in% sim.month[imonth]]
      }
      
      #Apply the value to the correct row
      if(length(value) == 1){
        #If only 1 value is supplied for multiple years need to only point to that value
        Rsim.scenario$forcing[[parameter]][year.row, group] <- value
      }else if(sim.month[1] == 0){
        Rsim.scenario$forcing[[parameter]][year.row, group] <- value[iyear]
      }else if(length(value) == length(sim.month)){
        Rsim.scenario$forcing[[parameter]][year.row, group] <- value[imonth]
      }else {
        Rsim.scenario$forcing[[parameter]][year.row, group] <- value[ivalue]
      }
    }
  }
  }
  
  return(Rsim.scenario)
}

# KYA - Some "adjustment" functions for Rsim.  Used "set" instead of "adjust" as
# that's more of a standard naming convention.  (In future make this a method
# of objects?)

#'Set Rsim.scenario parameters
#'
#'Modifies the various parameters of the \code{rsim.scenario()} object. Parameters 
#'that can be adjusted using this function are: \code{params},\code{start_state},
#'\code{forcing},\code{fishing},\code{stanzas}
#'
#'@family Adjust functions
#'
#'@inheritParams rsim.run
#'@inheritParams rsim.fishing
#'@param start_state Rsim starting values object generated by \code{rsim.state()}
#'@param forcing Rsim forcing matrix object generated by \code{rsim.forcing()} 
#'@param fishing Rsim fishing matrix object generated by \code{rsim.fishing()}
#'@param stanzas Rsim stanza parameters object generated by \code{rsim.stanzas()}
#'
#'@return Returns an \code{Rsim.scenario} object with the new parameter.
#'
#'@examples 
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' stanzas <- rsim.stanzas(AB.params)
#' Rsim.scenario.new <- set.rsim.scene(Rsim.scenario,stanzas=stanzas)
#'
#'
#'
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

#'Retrieve parameters from an Rsim scenario
#'
#'Helper function that will retrieve the parameters that were used 
#'in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a `params` object.
#'
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' params <- get.rsim.params(Rsim.scenario)
#' names(params)
#' 

#'
#'
#'@export
#'    
get.rsim.params<-function(Rsim.scenario){
  return(Rsim.scenario$params)
}

#'Retrieve starting state values from an Rsim scenario
#'
#'Helper function that will retrieve the starting state values that were used 
#'in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an `start_state` object.
#'
#'@examples
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' params <- get.rsim.start_state(Rsim.scenario)
#' names(params) 
#'
#'
#'
#'@export
#'   
get.rsim.start_state<-function(Rsim.scenario){
  return(Rsim.scenario$start_state)
}

#'Retrieve forcing parameters from an Rsim scenario
#'
#'Helper function that will retrieve the forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a `forcing` object.
#'
#'@examples 
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' params <- get.rsim.forcing(Rsim.scenario)
#' names(params)
#'
#'
#'@export
#'   
get.rsim.forcing<-function(Rsim.scenario){
  return(Rsim.scenario$forcing)
}

#'Retrieve fishing forcing parameters from an Rsim scenario
#'
#'Helper function that will retrieve the fishing forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a `fishing` object.
#'
#'@examples 
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' params <- get.rsim.fishing(Rsim.scenario)
#' names(params)
#'
#'
#'
#'@export
#'  
get.rsim.fishing<-function(Rsim.scenario){
  return(Rsim.scenario$fishing)
}

#'Retrieve stanza parameters from an Rsim scenario
#'
#'Helper function that will retrieve the stanza parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns a `stanzas` object.
#'
#'@examples 
#' # Read in Rpath parameter file and balance model
#' Rpath <- rpath(AB.params)
#' # Create a 50 yr Rpath scenario
#' Rsim.scenario <- rsim.scenario(Rpath, AB.params, years = 1:50)
#' params <- get.rsim.stanzas(Rsim.scenario)
#' names(params)
#'
#'
#'
#'@export
#' 
get.rsim.stanzas<-function(Rsim.scenario){
  return(Rsim.scenario$stanzas)
}







