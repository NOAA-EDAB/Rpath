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