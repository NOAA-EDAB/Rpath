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
#'@param nstanzas Numeric vector of the number of stanzas per multistanza groups.
#'
#'@return Outputs a list object of Rpath parameters which are populated with values 
#'  of NA or logical default values.  Values can then be filled in using
#'  R.  Use check.rpath.params() to ensure parameter files are filled out
#'  correctly (NOTE: This does not ensure data is correct just that it is 
#'  in the right places).
#'@import data.table
#'@export
create.rpath.params <- function(group, type, stgroup = NA){
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
  Rpath.params$pedigree <- pedigree
  
  return(Rpath.params)
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
check.rpath.params <- function(Rpath.params){
  w <- 0 #warning counter
  #Check to make sure all types are represented
  if(nrow(Rpath.params$model[Type == 0, ]) == 0){
    warning('Model must contain at least 1 consumer')
    w <- w + 1
  }
  if(nrow(Rpath.params$model[Type == 1, ]) == 0){
    warning('Model must contain a producer group')
    w <- w + 1
  }
  if(nrow(Rpath.params$model[Type == 2, ]) == 0){
    warning('Model must contain at least 1 detrital group')
    w <- w + 1
  }
  if(nrow(Rpath.params$model[Type == 3, ]) == 0){
    warning('Model must contain at least 1 fleet')
  }
  
  #Check that there is the proper number of columns
  n.groups <- nrow(Rpath.params$model)
  n.living <- length(Rpath.params$model[Type <= 1, Group])
  n.dead   <- length(Rpath.params$model[Type == 2, Group])
  n.fleet  <- length(Rpath.params$model[Type == 3, Group])
  if(ncol(Rpath.params$model) != 10 + n.dead + 2 * n.fleet){
    warning('Model does not have the correct number of column.  There should be 10 
         columns plus one for each detrital group plus two for each fleet group 
         (landings and discards).  Please double check your columns')
    w <- w + 1
  }
    
  #Check that either biomass or EE is entered and not both
  if(length(Rpath.params$model[is.na(Biomass) & is.na(EE) & Type < 2, Group]) > 0){
    warning(paste(Rpath.params$model[is.na(Biomass) & is.na(EE) & Type < 2, Group], 
               'are missing both Biomass and EE...must enter one \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[!is.na(Biomass) & !is.na(EE) & Type < 2, Group]) > 0){
    warning(paste(Rpath.params$model[!is.na(Biomass) & !is.na(EE) & Type < 2, Group],
               'have both Biomass and EE...Note that Rpath does not calculate BA 
               please enter a value for BA if appropriate \n', sep = ' '))
    w <- w + 1
  }

  #Check that Biomass / PB / QB / EE / ProdCons is not entered for types 3
  if(length(Rpath.params$model[Type == 3 & !is.na(Biomass), Group]) > 0){
    warning(paste(Rpath.params$model[Type == 3 & !is.na(Biomass), Group], 
               'is a fleet and should not have a biomass...set to NA \n', 
               sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type == 3 & !is.na(PB), Group]) > 0){
    warning(paste(Rpath.params$model[Type == 3 & !is.na(PB), Group], 
               'is a fleet and should not have a PB...set to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type > 1 & !is.na(QB), Group]) > 0){
    warning(paste(Rpath.params$model[Type > 1 & !is.na(QB), Group], 
               'are not living and should not have a QB...set to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type > 1 & !is.na(EE), Group]) > 0){
    warning(paste(Rpath.params$model[Type > 1 & !is.na(EE), Group], 
               'are not living and should not have a EE...set to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type > 1 & !is.na(ProdCons), Group]) > 0){
    warning(paste(Rpath.params$model[Type > 1 & !is.na(ProdCons), Group], 
               'are not living and should not have a ProdCons...set to NA \n', 
               sep = ' '))
    w <- w + 1
  }
  
  #Check that types 0 and 1 have a PB unless QB and ProdCons are entered
  if(length(Rpath.params$model[Type < 2 & is.na(PB), Group]) > 0){
    no.pb <- Rpath.params$model[Type < 2 & is.na(PB), Group]
    if(length(Rpath.params$model[Group %in% no.pb & (is.na(QB) | is.na(ProdCons)), Group]) > 0){
      warning(paste(Rpath.params$model[Group %in% no.pb & (is.na(QB) | is.na(ProdCons)), Group],
               'are missing a PB without a QB and PQ...set to >= 0 \n', sep = ' '))
    w <- w + 1
    }
    }
  
  #Check that consumers have a QB or ProdCons but not both unless missing PB
  if(length(Rpath.params$model[is.na(QB) & is.na(ProdCons) & Type < 1, Group]) > 0){
    warning(paste(Rpath.params$model[is.na(QB) & is.na(ProdCons) & Type < 1, Group], 
               'are missing both QB and ProdCons...must enter one \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[!is.na(QB) & !is.na(ProdCons) & Type < 1, Group]) > 0){
    both <- Rpath.params$model[!is.na(QB) & !is.na(ProdCons) & Type < 1, Group]
    if(length(Rpath.params$model[Group %in% both & !is.na(PB), Group]) > 0){
      warning(paste(Rpath.params$model[Group %in% both & !is.na(PB), Group],
                    'have PB, QB, and ProdCons...only two should be entered \n', 
                    sep = ' '))
      w <- w + 1  
    }
  }
  
  #Check that BioAcc / Unassim is NA for fleets and numeric for types < 3
  if(length(Rpath.params$model[Type == 3 & !is.na(BioAcc), Group]) > 0){
    warning(paste(Rpath.params$model[Type == 3 & !is.na(BioAcc), Group],
               'are fleets and should not have a BioAcc...set to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type == 3 & !is.na(Unassim), Group]) > 0){
    warning(paste(Rpath.params$model[Type == 3 & !is.na(Unassim), Group],
               'are fleets and should not have an Unassim...set to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type != 3 & is.na(BioAcc), Group]) > 0){
    warning(paste(Rpath.params$model[Type != 3 & is.na(BioAcc), Group],
               'must have a number for BioAcc...set to >= 0 \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type != 3 & is.na(Unassim), Group]) > 0){
    warning(paste(Rpath.params$model[Type != 3 & is.na(Unassim), Group],
               'must have a number for Unassim...set to >= 0 \n', sep = ' '))
    w <- w + 1
  }
  
  #Check that only Type 2 has DetInput set
  if(length(Rpath.params$model[Type != 2 & !is.na(DetInput), Group]) > 0){
    warning(paste(Rpath.params$model[Type != 2 & !is.na(DetInput), Group],
               'are not detritus...set DetInput to NA \n', sep = ' '))
    w <- w + 1
  }
  if(length(Rpath.params$model[Type == 2 & is.na(DetInput), Group]) >0){
    warning(paste(Rpath.params$model[Type == 2 & is.na(DetInput), Group],
               'are detritus...set DetInput to 0 \n', sep = ' '))
    w <- w + 1
  }
  
  #Check detritus fate is numeric and sum to 1
  det.matrix <- Rpath.params$model[, 11:(10 + n.dead), with = F]
  test.rows  <- rowSums(det.matrix)
  if(length(setdiff(which(Rpath.params$model[, Type] == 2), which(test.rows != 1))) > 0){
    warning(paste(Rpath.params$model[, Group][setdiff(which(Rpath.params$model[, Type] == 2), 
                                          which(test.rows != 1))],
               'detrital fate does not sum to 1 \n', sep = ' '))
    w <- w + 1
  }
  if(length(which(is.na(det.matrix))) > 0){
    na.group <- which(is.na(det.matrix))
    for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
      na.group[i] - n.groups
    na.group <- unique(na.group)
    warning(paste(Rpath.params$model[na.group, Group], 
               'one or more detrital fates are NA...set to >= 0 \n', sep = ' '))
    w <- w + 1
  }
  
  #Check that landings and discards are numbers for type < 3
  fleet.matrix <- Rpath.params$model[1:(n.groups - n.fleet), (11 + n.dead):ncol(Rpath.params$model), 
                            with = F]
  if(length(which(is.na(fleet.matrix))) > 0){
    na.group <- which(is.na(fleet.matrix))
    for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
      na.group[i] - n.groups
    na.group <- unique(na.group)
    warning(paste(Rpath.params$model[na.group, Group], 
               'one or more catches are NA...set to >= 0 \n', sep = ' '))
    w <- w + 1
  }
  
  #Check that fleets aren't catching other fleets
  fleet.matrix.2 <- Rpath.params$model[(n.groups - n.fleet + 1):n.groups, (11 + n.dead):
                                ncol(Rpath.params$model), with = F]
  if(length(which(!is.na(fleet.matrix.2))) > 0){
    not.na.group <- which(!is.na(fleet.matrix.2))
    for(i in 1:length(not.na.group)){
      while(not.na.group[i] > n.fleet){
        not.na.group[i] <- not.na.group[i] - n.fleet
        not.na.group <- unique(not.na.group)
      }
    }
    warning(paste(Rpath.params$model[not.na.group + (n.groups - n.fleet), Group], 
               'are catching another fleet...set to NA \n', sep = ' '))
    w <- w + 1
  }

  #Diet
  #Check that columns sum to 1
  col.names <- names(Rpath.params$diet)[2:ncol(Rpath.params$diet)]
  col.sums <- Rpath.params$diet[, lapply(.SD, sum, na.rm = T), .SDcols = col.names]
  
  #Check types (>0 & <=1 are primary producers)
  types <- Rpath.params$model[Type < 2, Type]
  dctype <- round(col.sums + types, 3)
  if(length(which(dctype != 1)) > 0){
    for(i in 1:length(which(dctype != 1))){
    warning(paste(col.names[which(dctype != 1)][i], 'sum,', 
                  col.sums[, which(dctype !=1)[i], with = F], 
                  'is not 1...check DC or proportion of primary production'))
      w <- w + 1
    }
  }

  #Check number of columns
  dietcol <- ncol(Rpath.params$diet)
  if(dietcol != (n.living + 1)){
    warning(paste(dietcol, ' is the incorrect number of columns in diet matrix.',
                  'There should be', n.living + 1))
    w <- w + 1
  }
  
  #Check that final row of diet is "Import"
  if(!Rpath.params$diet[nrow(Rpath.params$diet), 1] == 'Import' &
     !Rpath.params$diet[nrow(Rpath.params$diet), 1] == 'import'){
    
    warning('Diet matrix is missing the import row.  Please add "Import" as the 
            final row.  All entries can be 0 or NA.')
    w <- w + 1
  }

if(w == 0){
  cat('Rpath parameter file is functional')
} else {
    cat('Rpath parameter file needs attention!')
  }  
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
read.rpath.params <- function(modfile, dietfile, pedfile = NA,
                              stanzagroupfile = NA, stanzafile = NA){
  Rpath.params <- list()
  Rpath.params$model <- as.data.table(read.csv(modfile,  header = T))
  Rpath.params$diet  <- as.data.table(read.csv(dietfile, header = T))
  
  if(!is.na(stanzagroupfile)){
    stanzagroup <- as.data.table(read.csv(stanzagroupfile, header = T))
    Rpath.params$stanzas$NStanzaGroups <- nrow(stanzagroup)
    Rpath.params$stanzas$stgroups      <- stanzagroup
    Rpath.params$stanzas$stindiv       <- as.data.table(read.csv(stanzafile, 
                                                                 header = T))
  } else {
    Rpath.params$stanzas$NStanzaGroups <- 0
    Rpath.params$stanzas$stgroups      <- data.table(StGroupNum  = NA,
                                                     StanzaGroup = NA,
                                                     nstanzas    = NA,
                                                     VBGF_Ksp    = NA,
                                                     VBGF_d      = NA,
                                                     Wmat        = NA,
                                                     RecPower    = NA,
                                                     Wmat001     = NA,
                                                     Wmat50      = NA,
                                                     Amat001     = NA,
                                                     Amat50      = NA)
    Rpath.params$stanzas$stindiv       <- data.table(StGroupNum   = NA,
                                                     StanzaNum   = NA,
                                                     GroupNum    = NA,
                                                     Group       = NA,
                                                     First       = NA,
                                                     Last        = NA,
                                                     Z           = NA,
                                                     Leading     = NA)
  }
  if(!is.na(pedfile)){
    Rpath.params$pedigree <- as.data.table(read.csv(pedfile, header = T))
  } else {
    Rpath.params$pedigree <- data.table(Group = Rpath.params$model$Group,
                                        B     = 1,
                                        PB    = 1,
                                        QB    = 1,
                                        Diet  = 1)
    fleets <- as.character(Rpath.params$model[Type == 3, Group])
    for(i in 1:length(fleets)){
      Rpath.params$pedigree[, V1 := 1]
      setnames(Rpath.params$pedigree, 'V1', fleets[i])
    }
  }
  return(Rpath.params)
}

#'Write Rpath parameters to .csv files
#'
#'Creates a series of .csv files from an Rpath.params object.
#'
#'@family Rpath functions
#'
#'@param Rpath.params R object containing the Rpath parameters.  Most likely this
#'  was created using create.rpath.params or read.rpath.params.
#'@param eco.name ecosystem name that will be included in all the file names.
#'@param path location for the output files.  
#'@return Outputs a series of .csv files named by the provided eco.name and the 
#'  parameters they represent.  For example the model parameters will be named 
#'  "eco.name_model.csv".
#'@export
write.rpath.params <- function(Rpath.params, eco.name, path = ''){
    
  write.csv(Rpath.params$model,
            file = file.path(path, paste(eco.name, '_model.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.params$diet,
            file = file.path(path, paste(eco.name, '_diet.csv', sep = '')), 
            row.names = F)
  
  write.csv(Rpath.params$pedigree, 
            file = file.path(path, paste(eco.name, '_pedigree.csv', sep = '')), 
            row.names = F)
  
  #Multistanza parameters are in several different files
  write.csv(Rpath.params$stanzas$stgroups, 
            file = file.path(path, paste(eco.name, '_stanza_groups.csv', sep = '')),
            row.names = F)
  
  write.csv(Rpath.params$stanzas$stindiv,  
            file = file.path(path, paste(eco.name, '_stanzas.csv', sep = '')),
            row.names = F)
  
  if(Rpath.params$stanzas$NStanzaGroups > 0){
    for(isp in 1:Rpath.params$stanzas$NStanzaGroups){
      write.csv(Rpath.params$stanzas$StGroup[[isp]], 
                file = file.path(path, paste(eco.name, '_', 
                             Rpath.params$stanzas$stgroups$StanzaGroup[isp], 
                             '.csv', sep = '')), row.names = F)
    }
  }
}
