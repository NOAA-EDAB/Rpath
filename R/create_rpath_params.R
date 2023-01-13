## Rpath functions for creating and checking parameter files  

#'Creates a shell for the Rpath.params object
#'
#'Creates a shell of the Rpath.params list object which contains the model, diet,
#'multistanza, and pedigree parameters.
#'
#'@family Rpath functions
#'
#'@param group Vector of group names.
#'@param type Numeric vector of group type. Living = 0, Producer = 1, Detritus = 2,
#'  Fleet = 3.
#'@param stgroup Vector of multistanza group names.  Include NA for non-stanza groups.
#'
#'@return Outputs a list object of Rpath parameters which are populated with values 
#'  of NA or logical default values.  Values can then be filled in using
#'  R.  Use check.rpath.params() to ensure parameter files are filled out
#'  correctly (NOTE: This does not ensure data is correct just that it is 
#'  in the right places).
#'@import data.table
#'@export
create.rpath.params <- function(group, type, stgroup = NA){
  #Need to define variables to eliminate check() note about no visible binding
  Group <- DetInput <- V1 <- StGroupNum <- NULL
  
  Rpath.params <- list()
  
  pred.group  <- group[which(type < 2)]
  prey.group  <- group[which(type < 3)]
  det.group   <- group[which(type == 2)]
  fleet.group <- group[which(type == 3)]
  
  #Model parameters
  model <- data.table(Group       = group, 
                      Type        = type, 
                      Biomass     = as.numeric(NA),
                      PB          = as.numeric(NA),
                      QB          = as.numeric(NA),
                      EE          = as.numeric(NA),
                      ProdCons    = as.numeric(NA),
                      BioAcc      = as.numeric(NA),
                      Unassim     = as.numeric(NA),
                      DetInput    = as.numeric(NA))

  #Add detritial groups
  for(i in 1:length(det.group)){
    model[Group %in% det.group, DetInput := 0]
    model[, V1 := as.numeric(NA)]
    setnames(model, "V1", det.group[i])
  }
    
  #Add fleets twice - Landings and Discards
  for(i in 1:length(fleet.group)){
    model[, V1 := c(rep(0, length(group) - length(fleet.group)), 
                    rep(NA, length(fleet.group)))]
    setnames(model, "V1", fleet.group[i])
  }
  for(i in 1:length(fleet.group)){
    model[, V1 := c(rep(0, length(group) - length(fleet.group)), 
                    rep(NA, length(fleet.group)))]
    setnames(model, "V1", paste(fleet.group[i], '.disc', sep = ''))
  }
  Rpath.params$model <- model
  
  #Diet matrix
  diet <- data.table(Group = c(prey.group, 'Import'))
  for(i in 1:length(pred.group)){
    diet[, V1 := as.numeric(NA)]
    setnames(diet, "V1", pred.group[i])
  }
  Rpath.params$diet <- diet
  
  #Multistanza parameters
  if(length(stgroup) > 1){
    #Group Parameters
    StanzaGroups  <- unique(stgroup[!is.na(stgroup)])
    nstanzas      <- as.vector(table(stgroup)[StanzaGroups])
    NStanzaGroups <- length(StanzaGroups)
    Rpath.params$stanzas$NStanzaGroups <- NStanzaGroups
    
    stgroups <- data.table(StGroupNum  = 1:NStanzaGroups,
                           StanzaGroup = StanzaGroups,
                           nstanzas    = nstanzas,
                           VBGF_Ksp    = NA,
                           VBGF_d      = 0.66667,
                           Wmat        = NA,
                           BAB         = 0,
                           RecPower    = 1)
    
    #Individual Stanza Parameters
    ind.stanza.group <- model[!is.na(stgroup), Group]
    ieco <- which(!is.na(stgroup))
    stindiv <- data.table(StGroupNum = rep(stgroups[, StGroupNum], 
                                           stgroups[, nstanzas]),
                          StanzaNum  = as.integer(0),
                          GroupNum   = ieco,
                          Group      = ind.stanza.group,
                          First      = NA,
                          Last       = NA,
                          Z          = NA,
                          Leading    = NA)
  
    } else {
      
    Rpath.params$stanzas$NStanzaGroups <- 0
    
    stgroups <- data.table(StGroupNum  = NA,
                           StanzaGroup = NA,
                           nstanzas    = NA,
                           VBGF_Ksp    = NA,
                           VBGF_d      = NA,
                           Wmat        = NA,
                           RecPower    = NA)
    
    stindiv <- data.table(StGroupNum = NA,
                          StanzaNum  = NA,
                          GroupNum   = NA,
                          Group      = NA,
                          First      = NA,
                          Last       = NA,
                          Z          = NA,
                          Leading    = NA)
  }
  
  Rpath.params$stanzas$stgroups <- stgroups
  Rpath.params$stanzas$stindiv  <- stindiv
    
  #Pedigree
  pedigree <- data.table(Group   = group,
                         Biomass = 1,
                         PB      = 1,
                         QB      = 1,
                         Diet    = 1)
  #Add fleet pedigree
  for(i in 1:length(fleet.group)){
    pedigree[, V1 := 1]
    setnames(pedigree, "V1", fleet.group[i])
  }
  Rpath.params$pedigree <- pedigree
  class(Rpath.params) <- 'Rpath.params'
  return(Rpath.params)
}
