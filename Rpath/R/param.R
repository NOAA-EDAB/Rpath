## Rpath functions for creating and checking parameter files  

#'Creates a shell for the Rpath.param object
#'
#'Creates a shell of the Rpath.param list object which contains the model, diet,
#'multistanza, and pedigree parameters.
#'
#'@family Rpath functions
#'
#'@param group Vector of group names.
#'@param type Numeric vector of group type. Living = 0, Producer = 1, Detritus = 2,
#'  Fleet = 3.
#'@param stgroup Vector of multistanza group names.  Include NA for non-stanza groups.
#'@param nstanzas Numeric vector of the number of stanzas per multistanza groups.
#'
#'@return Outputs a list object of Rpath parameters which are populated with values 
#'  of NA or logical default values.  Values can then be filled in using
#'  R.  Use check.rpath.param() to ensure parameter files are filled out
#'  correctly (NOTE: This does not ensure data is correct just that it is 
#'  in the right places).
#'@import data.table
#'@export
create.rpath.param <- function(group, type, stgroup = NA){
  Rpath.param <- list()
  
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
  Rpath.param$model <- model
  
  #Diet matrix
  diet <- data.table(Group = prey.group)
  for(i in 1:length(pred.group)){
    diet[, V1 := as.numeric(NA)]
    setnames(diet, "V1", pred.group[i])
  }
  Rpath.param$diet <- diet
  
  #Multistanza parameters
  if(length(stgroup) > 1){
    #Group Parameters
    StanzaGroups  <- unique(stgroup[!is.na(stgroup)])
    nstanzas      <- as.vector(table(stgroups)[StanzaGroups])
    NStanzaGroups <- length(StanzaGroups)
    Rpath.param$stanzas$NStanzaGroups <- NStanzaGroups
    
    stgroups <- data.table(StGroupNum  = 1:NStanzaGroups,
                           StanzaGroup = StanzaGroups,
                           nstanzas    = nstanzas,
                           VBGF_Ksp    = NA,
                           VBGF_d      = 0.66667,
                           Wmat        = NA,
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
      
    Rpath.param$stanzas$NStanzaGroups <- 0
    
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
  
  Rpath.param$stanzas$stgroups <- stgroups
  Rpath.param$stanzas$stindiv  <- stindiv
    
  #Pedigree
  pedigree <- data.table(Group = group,
                         B     = 1,
                         PB    = 1,
                         QB    = 1,
                         Diet  = 1)
  #Add fleet pedigree
  for(i in 1:length(fleet.group)){
    pedigree[, V1 := 1]
    setnames(pedigree, "V1", fleet.group[i])
  }
  Rpath.param$pedigree <- pedigree
  
  return(Rpath.param)
}


#'Check Rpath parameter files
#'
#'Logical check that the parameter files are filled out correctly, i.e. data is entered where it is
#'expected.
#'
#'@family Rpath functions
#'
#'@param filename Name of the parameter file.  Can be the path to a .csv or an R
#'object. 
#'@param parameter The type of parameter file you are checking.  Choices include 
#"model", "diet", "juvenile", and "pedigree".
#'
#'@return Checks Rpath parameter files for consistency.  An error message will be produced if one of
#'  the logical checks fails.  Checks include: 
#'  (NOTE: This does not ensure data is correct just that it is in the right places).
#'@import data.table
#'@export
check.rpath.param <- function(Rpath.param){

  #Check to make sure all types are represented
  if(length(Rpath.param$model[Type == 0, ]) == 0) stop('Model must contain at least 1 consumer')
  if(length(Rpath.param$model[Type == 1, ]) == 0) stop('Model must contain a producer group')
  if(length(Rpath.param$model[Type == 2, ]) == 0) stop('Model must contain at least 1 detrital group')
  if(length(Rpath.param$model[Type == 3, ]) == 0) stop('Model must contain at least 1 fleet')
  
  #Check that there is the proper number of columns
  n.groups <- nrow(Rpath.param$model)
  n.dead   <- length(Rpath.param$model[Type == 2, Group])
  n.fleet  <- length(Rpath.param$model[Type == 3, Group])
  if(ncol(Rpath.param$model) != 10 + n.dead + 2 * n.fleet){
    stop('Model does not have the correct number of column.  There should be 10 
         columns plus one for each detrital group plus two for each fleet group 
         (landings and discards).  Please double check your columns')
  }
    
  #Check that either biomass or EE is entered and not both
  if(length(Rpath.param$model[is.na(Biomass) & is.na(EE) & Type < 2, Group]) > 0){
    stop(paste(Rpath.param$model[is.na(Biomass) & is.na(EE) & Type < 2, Group], 
               'are missing both Biomass and EE...must enter one \n', sep = ' '))
  }
  if(length(Rpath.param$model[!is.na(Biomass) & !is.na(EE) & Type < 2, Group]) > 0){
    stop(paste(Rpath.param$model[!is.na(Biomass) & !is.na(EE) & Type < 2, Group],
               'have both Biomass and EE...only one should be entered \n', sep = ' '))
    #consider making this a warning() since you can have both in EwE
  }

  #Check that Biomass / PB / QB / EE / ProdCons is not entered for types 3
  if(length(Rpath.param$model[Type == 3 & !is.na(Biomass), Group]) > 0){
    stop(paste(Rpath.param$model[Type == 3 & !is.na(Biomass), Group], 
               'is a fleet and should not have a biomass...set to NA \n', 
               sep = ' '))
  }
  if(length(Rpath.param$model[Type == 3 & !is.na(PB), Group]) > 0){
    stop(paste(Rpath.param$model[Type == 3 & !is.na(PB), Group], 
               'is a fleet and should not have a PB...set to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type > 1 & !is.na(QB), Group]) > 0){
    stop(paste(Rpath.param$model[Type > 1 & !is.na(QB), Group], 
               'are not living and should not have a QB...set to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type > 1 & !is.na(EE), Group]) > 0){
    stop(paste(Rpath.param$model[Type > 1 & !is.na(EE), Group], 
               'are not living and should not have a EE...set to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type > 1 & !is.na(ProdCons), Group]) > 0){
    stop(paste(Rpath.param$model[Type > 1 & !is.na(ProdCons), Group], 
               'are not living and should not have a ProdCons...set to NA \n', 
               sep = ' '))
  }
  
  #Check that types 0 and 1 have a PB
  if(length(Rpath.param$model[Type < 2 & is.na(PB), Group]) > 0){
    stop(paste(Rpath.param$model[Type < 2 & is.na(PB), Group],
               'are missing a PB...set to >= 0 \n', sep = ' '))
  }
  
  #Check that consumers have a QB or ProdCons but not both
  if(length(Rpath.param$model[is.na(QB) & is.na(ProdCons) & Type < 1, Group]) > 0){
    stop(paste(Rpath.param$model[is.na(QB) & is.na(ProdCons) & Type < 1, Group], 
               'are missing both QB and ProdCons...must enter one \n', sep = ' '))
  }
  if(length(Rpath.param$model[!is.na(QB) & !is.na(ProdCons) & Type < 1, Group]) > 0){
    stop(paste(Rpath.param$model[!is.na(QB) & !is.na(ProdCons) & Type < 1, Group],
               'have both QB and ProdCons...only one should be entered \n', 
               sep = ' '))
  }
  
  #Check that BioAcc / Unassim is NA for fleets and numeric for types < 3
  if(length(Rpath.param$model[Type == 3 & !is.na(BioAcc), Group]) > 0){
    stop(paste(Rpath.param$model[Type == 3 & !is.na(BioAcc), Group],
               'are fleets and should not have a BioAcc...set to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type == 3 & !is.na(Unassim), Group]) > 0){
    stop(paste(Rpath.param$model[Type == 3 & !is.na(Unassim), Group],
               'are fleets and should not have an Unassim...set to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type != 3 & is.na(BioAcc), Group]) > 0){
    stop(paste(Rpath.param$model[Type != 3 & is.na(BioAcc), Group],
               'must have a number for BioAcc...set to >= 0 \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type != 3 & is.na(Unassim), Group]) > 0){
    stop(paste(Rpath.param$model[Type != 3 & is.na(Unassim), Group],
               'must have a number for Unassim...set to >= 0 \n', sep = ' '))
  }
  
  #Check that only Type 2 has DetInput set
  if(length(Rpath.param$model[Type != 2 & !is.na(DetInput), Group]) > 0){
    stop(paste(Rpath.param$model[Type != 2 & !is.na(DetInput), Group],
               'are not detritus...set DetInput to NA \n', sep = ' '))
  }
  if(length(Rpath.param$model[Type == 2 & is.na(DetInput), Group]) >0){
    stop(paste(Rpath.param$model[Type == 2 & is.na(DetInput), Group],
               'are detritus...set DetInput to 0 \n', sep = ' '))
  }
  
  #Check detritus fate is numeric and sum to 1
  det.matrix <- Rpath.param$model[, 11:(10 + n.dead), with = F]
  test.rows  <- rowSums(det.matrix)
  if(length(setdiff(which(Rpath.param$model[, Type] == 2), which(test.rows != 1))) > 0){
    stop(paste(Rpath.param$model[, Group][setdiff(which(Rpath.param$model[, Type] == 2), 
                                          which(test.rows != 1))],
               'detrital fate does not sum to 1 \n', sep = ' '))
  }
  if(length(which(is.na(det.matrix))) > 0){
    na.group <- which(is.na(det.matrix))
    for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
      na.group[i] - n.groups
    na.group <- unique(na.group)
    stop(paste(Rpath.param$model[na.group, Group], 
               'one or more detrital fates are NA...set to >= 0 \n', sep = ' '))
  }
  
  #Check that landings and discards are numbers for type < 3
  fleet.matrix <- Rpath.param$model[1:(n.groups - n.fleet), (11 + n.dead):ncol(Rpath.param$model), 
                            with = F]
  if(length(which(is.na(fleet.matrix))) > 0){
    na.group <- which(is.na(fleet.matrix))
    for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
      na.group[i] - n.groups
    na.group <- unique(na.group)
    stop(paste(Rpath.param$model[na.group, Group], 
               'one or more catches are NA...set to >= 0 \n', sep = ' '))
  }
  
  #Check that fleets aren't catching other fleets
  fleet.matrix.2 <- Rpath.param$model[(n.groups - n.fleet + 1):n.groups, (11 + n.dead):
                                ncol(Rpath.param$model), with = F]
  if(length(which(!is.na(fleet.matrix.2))) > 0){
    not.na.group <- which(!is.na(fleet.matrix.2))
    for(i in 1:length(not.na.group)){
      while(not.na.group[i] > n.fleet){
        not.na.group[i] <- not.na.group[i] - n.fleet
        not.na.group <- unique(not.na.group)
      }
    }
    stop(paste(Rpath.param$model[not.na.group + (n.groups - n.fleet), Group], 
               'are catching another fleet...set to NA \n', sep = ' '))
  }

  #Diet
  #Check that columns sum to 1
  col.names <- names(Rpath.param$diet)[2:ncol(Rpath.param$diet)]
  col.sums <- Rpath.param$diet[, lapply(.SD, sum, na.rm = T), .SDcols = col.names]
  #Check types (>0 & <=1 are primary producers)
  types <- Rpath.param$model[Type < 2, Type]
  dctype <- col.sums + types
  if(length(which(dctype != 1)) > 0){
    for(i in 1:length(which(dctype != 1))){
    stop(paste(col.names[which(dctype != 1)][i], 'sum,', 
                  col.sums[, which(dctype !=1)[i], with = F], 
                  'is not 1...check DC or proportion of primary production'))
    }
  }

cat('Rpath parameter file is functional')
}

#'Read Rpath parameters from .csv files
#'
#'Creates an Rpath.param object from a series of .csv files.
#'
#'@family Rpath functions
#'
#'@param modfile file location of the flat file containing the model parameters.
#'@param dietfile file location of the flat file containing the diet parameters.
#'@param stanzagroupfile file location of the flat file containing the group parameters
#'  for multistanza groups.  If not specified a blank stanza list will be created.  
#'@param stanzafile file location of the flat file containing the individual stanza 
#'  parameters for multistanza groups.  If not specified a blank stanza list will 
#'  be created.
#'@param pedfile file location of the flat file containg the pedgigree parameters.
#'@return Outputs an Rpath.param object that can be used for Rpath and subsequently
#'  Rsim.  (NOTE: This does function does not ensure data is correct or in the 
#'  correct locations...run check.rpath.param to ensure the appropriate columns are
#'  present).
#'@export
read.rpath.param <- function(modfile, dietfile, pedfile = NA,
                              stanzagroupfile = NA, stanzafile = NA){
  Rpath.param <- list()
  Rpath.param$model <- as.data.table(read.csv(modfile,  header = T))
  Rpath.param$diet  <- as.data.table(read.csv(dietfile, header = T))
  
  if(!is.na(stanzagroupfile)){
    stanzagroup <- as.data.table(read.csv(stanzagroupfile, header = T))
    Rpath.param$stanzas$NStanzaGroups <- nrow(stanzagroup)
    Rpath.param$stanzas$stgroups      <- stanzagroup
    Rpath.param$stanzas$stindiv       <- as.data.table(read.csv(stanzafile, 
                                                                 header = T))
  } else {
    Rpath.param$stanzas$NStanzaGroups <- 0
    Rpath.param$stanzas$stgroups      <- data.table(StGroupNum  = NA,
                                                     StanzaGroup = NA,
                                                     nstanzas    = NA,
                                                     VBGF_Ksp    = NA,
                                                     VBGF_d      = NA,
                                                     Wmat        = NA,
                                                     RecPower    = NA)
    Rpath.param$stanzas$stindiv       <- data.table(StGroupNum   = NA,
                                                     StanzaNum   = NA,
                                                     GroupNum    = NA,
                                                     Group       = NA,
                                                     First       = NA,
                                                     Last        = NA,
                                                     Z           = NA,
                                                     Leading     = NA)
  }
  if(!is.na(pedfile)){
    Rpath.param$pedigree <- as.data.table(read.csv(pedfile, header = T))
  } else {
    Rpath.param$pedigree <- data.table(Group = Rpath.param$model$Group,
                                        B     = 1,
                                        PB    = 1,
                                        QB    = 1,
                                        Diet  = 1)
    fleets <- Rpath.param$model[Type == 3, Group]
    for(i in 1:length(fleets)){
      Rpath.param$pedigree[, V1 := 1]
      setnames(Rpath.param$pedigree, 'V1', fleets[i])
    }
  }
  return(Rpath.param)
}

#'Write Rpath parameters to .csv files
#'
#'Creates a series of .csv files from an Rpath.param object.
#'
#'@family Rpath functions
#'
#'@param Rpath.param R object containing the Rpath parameters.  Most likely this
#'  was created using create.rpath.param or read.rpath.param.
#'@param eco.name ecosystem name that will be included in all the file names.
#'@param path location for the output files.  
#'@return Outputs a series of .csv files named by the provided eco.name and the 
#'  parameters they represent.  For example the model parameters will be named 
#'  "eco.name_model.csv".
#'@export
write.rpath.param <- function(Rpath.param, eco.name, path = ''){
    
  write.csv(Rpath.param$model,
            file = file.path(path, paste(eco.name, '_model.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.param$diet,
            file = file.path(path, paste(eco.name, '_diet.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.param$pedigree, 
            file = file.path(path, paste(eco.name, '_pedigree.csv', sep = '')), 
            row.names = F)
  
  #Multistanza parameters are in several different files
  write.csv(Rpath.param$stanzas$stgroups, 
            file = file.path(path, paste(eco.name, '_stanza_groups.csv', sep = '')),
            row.names = F)
  
  write.csv(Rpath.param$stanzas$stindiv,  
            file = file.path(path, paste(eco.name, '_stanzas.csv', sep = '')),
            row.names = F)
  
  if(Rpath.param$stanzas$NStanzaGroups > 0){
    for(isp in 1:Rpath.param$stanzas$NStanzaGroups){
      write.csv(Rpath.param$stanzas$StGroup[[isp]], 
                file = file.path(path, paste(eco.name, '_', 
                             Rpath.param$stanzas$stgroups$StanzaGroup[isp], 
                             '.csv', sep = '')), row.names = F)
    }
  }
}
