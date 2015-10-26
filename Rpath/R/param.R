## Rpath functions for creating and checking parameter files  

#'Create shells for the Rpath parameter files
#'
#'Creates a shell of the parameter files that can then be filled out in R, Excel, or
#'another spreadsheet program.
#'
#'@family Rpath functions
#'
#'@param group Vector of group names.  If parameter equals "juvenile", this should be a vector of 
#'  stanza groups only (Juvenile and Adults in one group).
#'@param type Numeric vector of group type. Living = 0, Producer = 1, Detritus = 2,
#'  Fleet = 3. Default NA is used for the juvenile and pedigree parameter files.
#'@param filename Name of the output file saved as a .csv. If NA the file will not be written.
#'@param parameter The type of parameter file you are creating.  Choices include "model",
#'  "diet", "juvenile", and "pedigree".
#'
#'@return Outputs a shell of the parameter file indicated by the parameter variable.  The shell
#'  is populated with values of NA or logical default values.  Values can then be filled in using
#'  R or any spreadsheet program.  Use check.rpath.param() to ensure parameter files are filled out
#'  correctly (NOTE: This does not ensure data is correct just that it is in the right places).
#'@import data.table
#'@export
create.rpath.param <- function(parameter = 'model', group = NA, type = NA, 
                               filename = NA){
  pred.group  <- group[which(type < 2)]
  prey.group  <- group[which(type < 3)]
  det.group   <- group[which(type == 2)]
  fleet.group <- group[which(type == 3)]
  
  #Base model
  if(parameter == 'model'){
    out <- data.table(Group       = group, 
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
      out[Group %in% det.group, DetInput := 0]
      out[, V1 := as.numeric(NA)]
      setnames(out, "V1", det.group[i])
    }
    
    #Add fleets twice - Landings and Discards
    for(i in 1:length(fleet.group)){
      out[, V1 := c(rep(0, length(group) - length(fleet.group)), rep(NA, length(fleet.group)))]
      setnames(out, "V1", fleet.group[i])
    }
    for(i in 1:length(fleet.group)){
      out[, V1 := c(rep(0, length(group) - length(fleet.group)), rep(NA, length(fleet.group)))]
      setnames(out, "V1", paste(fleet.group[i], '.disc', sep = ''))
    }
  }
  
  #Diet matrix
  if(parameter == 'diet'){
    out <- data.table(Group = prey.group)
    for(i in 1:length(pred.group)){
      out[, V1 := as.numeric(NA)]
      setnames(out, "V1", pred.group[i])
    }
  }
  
  #Juvenile file
  if(parameter == 'juvenile'){
    out <- data.table(StanzaName = group,
                      JuvNum     = as.numeric(NA),
                      AduNum     = as.numeric(NA),
                      RecAge     = as.numeric(NA),
                      RecMonth   = as.numeric(NA),
                      VonBK      = as.numeric(NA),
                      AduZ_BAB   = as.numeric(NA),
                      JuvZ_BAB   = as.numeric(NA),
                      VonBD      = 0.6667,
                      Wmat50     = as.numeric(NA),
                      Wmat001    = as.numeric(NA),
                      Amat50     = as.numeric(NA),
                      Amat001    = as.numeric(NA),
                      RecPower   = 1)
  }
  
  #Pedigree
  if(parameter == 'pedigree'){
    out <- data.table(Group = group,
                      B     = 1,
                      PB    = 1,
                      QB    = 1,
                      Diet  = 1)
    #Add fleet pedigree
    for(i in 1:length(fleet.group)){
      out[, V1 := 1]
      setnames(out, "V1", fleet.group[i])
    }
  }
if(!is.na(filename)){
  write.csv(out, file = filename, row.names = F)
}else{
  return(out)
}
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
#'@param type The type of parameter file you are checking.  Choices include "model",
#'  "diet", "juvenile", and "pedigree".
#'
#'@return Checks Rpath parameter files for consistency.  An error message will be produced if one of
#'  the logical checks fails.  Checks include: 
#'  (NOTE: This does not ensure data is correct just that it is in the right places).
#'@import data.table
#'@export
check.rpath.param <- function(filename, type = 'model'){
  if(is.character(filename)){
    parameter <- as.data.table(read.csv(filename))
  } else {
    parameter <- as.data.table(filename)
  }
  
  if(type == 'model'){
    #Check to make sure all types are represented
    if(length(parameter[Type == 0, ]) == 0) stop('Model must contain at least 1 consumer')
    if(length(parameter[Type == 1, ]) == 0) stop('Model must contain a producer group')
    if(length(parameter[Type == 2, ]) == 0) stop('Model must contain at least 1 detrital group')
    if(length(parameter[Type == 3, ]) == 0) stop('Model must contain at least 1 fleet')
    
    #Check that there is the proper number of columns
    n.groups <- nrow(parameter)
    n.dead   <- length(parameter[Type == 2, Group])
    n.fleet  <- length(parameter[Type == 3, Group])
    if(ncol(parameter) != 10 + n.dead + 2 * n.fleet){
      stop('Model does not have the correct number of column.  There should be 10 
           columns plus one for each detrital group plus two for each fleet group 
           (landings and discards).  Please double check your columns')
    }
      
    #Check that either biomass or EE is entered and not both
    if(length(parameter[is.na(Biomass) & is.na(EE) & Type < 2, Group]) > 0){
      stop(paste(parameter[is.na(Biomass) & is.na(EE) & Type < 2, Group], 
                 'are missing both Biomass and EE...must enter one \n', sep = ' '))
    }
    if(length(parameter[!is.na(Biomass) & !is.na(EE) & Type < 2, Group]) > 0){
      stop(paste(parameter[!is.na(Biomass) & !is.na(EE) & Type < 2, Group],
                 'have both Biomass and EE...only one should be entered \n', sep = ' '))
      #consider making this a warning() since you can have both in EwE
    }
  
    #Check that Biomass / PB / QB / EE / ProdCons is not entered for types 2 and 3
    if(length(parameter[Type > 1 & !is.na(Biomass), Group]) > 0){
      stop(paste(parameter[Type > 1 & !is.na(Biomass), Group], 
                 'are not living and should not have a biomass...set to NA \n', 
                 sep = ' '))
    }
    if(length(parameter[Type > 1 & !is.na(PB), Group]) > 0){
      stop(paste(parameter[Type > 1 & !is.na(PB), Group], 
                 'are not living and should not have a PB...set to NA \n', sep = ' '))
    }
    if(length(parameter[Type > 1 & !is.na(QB), Group]) > 0){
      stop(paste(parameter[Type > 1 & !is.na(QB), Group], 
                 'are not living and should not have a QB...set to NA \n', sep = ' '))
    }
    if(length(parameter[Type > 1 & !is.na(EE), Group]) > 0){
      stop(paste(parameter[Type > 1 & !is.na(EE), Group], 
                 'are not living and should not have a EE...set to NA \n', sep = ' '))
    }
    if(length(parameter[Type > 1 & !is.na(ProdCons), Group]) > 0){
      stop(paste(parameter[Type > 1 & !is.na(ProdCons), Group], 
                 'are not living and should not have a ProdCons...set to NA \n', 
                 sep = ' '))
    }
    
    #Check that types 0 and 1 have a PB
    if(length(parameter[Type < 2 & is.na(PB), Group]) > 0){
      stop(paste(parameter[Type < 2 & is.na(PB), Group],
                 'are missing a PB...set to >= 0 \n', sep = ' '))
    }
    
    #Check that consumers have a QB or ProdCons but not both
    if(length(parameter[is.na(QB) & is.na(ProdCons) & Type == 0, Group]) > 0){
      stop(paste(parameter[is.na(QB) & is.na(ProdCons) & Type == 0, Group], 
                 'are missing both QB and ProdCons...must enter one \n', sep = ' '))
    }
    if(length(parameter[!is.na(QB) & !is.na(ProdCons) & Type == 0, Group]) > 0){
      stop(paste(parameter[!is.na(QB) & !is.na(ProdCons) & Type == 0, Group],
                 'have both QB and ProdCons...only one should be entered \n', 
                 sep = ' '))
    }
    
    #Check that BioAcc / Unassim is NA for fleets and numeric for types < 3
    if(length(parameter[Type == 3 & !is.na(BioAcc), Group]) > 0){
      stop(paste(parameter[Type == 3 & !is.na(BioAcc), Group],
                 'are fleets and should not have a BioAcc...set to NA \n', sep = ' '))
    }
    if(length(parameter[Type == 3 & !is.na(Unassim), Group]) > 0){
      stop(paste(parameter[Type == 3 & !is.na(Unassim), Group],
                 'are fleets and should not have an Unassim...set to NA \n', sep = ' '))
    }
    if(length(parameter[Type != 3 & is.na(BioAcc), Group]) > 0){
      stop(paste(parameter[Type != 3 & is.na(BioAcc), Group],
                 'must have a number for BioAcc...set to >= 0 \n', sep = ' '))
    }
    if(length(parameter[Type != 3 & is.na(Unassim), Group]) > 0){
      stop(paste(parameter[Type != 3 & is.na(Unassim), Group],
                 'must have a number for Unassim...set to >= 0 \n', sep = ' '))
    }
    
    #Check that only Type 2 has DetInput set
    if(length(parameter[Type != 2 & !is.na(DetInput), Group]) > 0){
      stop(paste(parameter[Type != 2 & !is.na(DetInput), Group],
                 'are not detritus...set DetInput to NA \n', sep = ' '))
    }
    if(length(parameter[Type == 2 & is.na(DetInput), Group]) >0){
      stop(paste(parameter[Type == 2 & is.na(DetInput), Group],
                 'are detritus...set DetInput to 0 \n', sep = ' '))
    }
    
    #Check detritus fate is numeric and sum to 1
    det.matrix <- parameter[, 11:(10 + n.dead), with = F]
    test.rows  <- rowSums(det.matrix)
    if(length(setdiff(which(parameter[, Type] == 2), which(test.rows != 1))) > 0){
      stop(paste(parameter[, Group][setdiff(which(parameter[, Type] == 2), 
                                            which(test.rows != 1))],
                 'detrital fate does not sum to 1 \n', sep = ' '))
    }
    if(length(which(is.na(det.matrix))) > 0){
      na.group <- which(is.na(det.matrix))
      for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
        na.group[i] - n.groups
      na.group <- unique(na.group)
      stop(paste(parameter[na.group, Group], 
                 'one or more detrital fates are NA...set to >= 0 \n', sep = ' '))
    }
    
    #Check that landings and discards are numbers for type < 3
    fleet.matrix <- parameter[1:(n.groups - n.fleet), (11 + n.dead):ncol(parameter), 
                              with = F]
    if(length(which(is.na(fleet.matrix))) > 0){
      na.group <- which(is.na(fleet.matrix))
      for(i in 1:length(na.group)) while(na.group[i] > n.groups) na.group[i] <- 
        na.group[i] - n.groups
      na.group <- unique(na.group)
      stop(paste(parameter[na.group, Group], 
                 'one or more detrital fates are NA...set to >= 0 \n', sep = ' '))
    }
    
    #Check that fleets aren't catching other fleets
    fleet.matrix.2 <- parameter[(n.groups - n.fleet + 1):n.groups, (11 + n.dead):
                                  ncol(parameter), with = F]
    if(length(which(!is.na(fleet.matrix.2))) > 0){
      not.na.group <- which(!is.na(fleet.matrix.2))
      for(i in 1:length(not.na.group)){
        while(not.na.group[i] > n.fleet){
          not.na.group[i] <- not.na.group[i] - n.fleet
          not.na.group <- unique(not.na.group)
        }
      }
      stop(paste(parameter[not.na.group + (n.groups - n.fleet), Group], 
                 'are catching another fleet...set to NA \n', sep = ' '))
    }
  
  
  }
  
cat(paste(type, 'parameter file is functional'))
}


