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
create.rpath.param <- function(filename = NA, group, type = NA, parameter = 'model'){
  pred.group  <- group[which(type < 2)]
  prey.group  <- group[which(type < 3)]
  det.group   <- group[which(type == 2)]
  fleet.group <- group[which(type == 3)]
  
  #Base model
  if(parameter == 'model'){
    out <- data.table(Group    = group, 
                      Type     = type, 
                      Biomass  = NA,
                      PB       = NA,
                      QB       = NA,
                      EE       = NA,
                      ProdCons = NA,
                      BioAcc   = NA,
                      Unassim  = NA,
                      DetInput = NA)
    
    #Add detritial groups
    for(i in 1:length(det.group)){
      out[, V1 := NA]
      setnames(out, "V1", det.group[i])
    }
    
    #Add fleets twice - Landings and Discards
    for(i in 1:length(fleet.group)){
      out[, V1 := NA]
      setnames(out, "V1", fleet.group[i])
    }
    for(i in 1:length(fleet.group)){
      out[, V1 := NA]
      setnames(out, "V1", fleet.group[i])
    }
  }
  
  #Diet matrix
  if(parameter == 'diet'){
    out <- data.table(Group = prey.group)
    for(i in 1:length(pred.group)){
      out[, V1 := NA]
      setnames(out, "V1", pred.group[i])
    }
  }
  
  #Juvenile file
  if(parameter == 'juvenile'){
    out <- data.table(StanzaName = group,
                      JuvNum     = NA,
                      AduNum     = NA,
                      RecAge     = NA,
                      RecMonth   = NA,
                      VonBK      = NA,
                      AduZ_BAB   = NA,
                      JuvZ_BAB   = NA,
                      VonBD      = .6667,
                      Wmat50     = NA,
                      Wmat001    = NA,
                      Amat50     = NA,
                      Amat001    = NA,
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
#'@param filename Name of the parameter file.  Must be a .csv. 
#'@param parameter The type of parameter file you are checking.  Choices include "model",
#'  "diet", "juvenile", and "pedigree".
#'
#'@return Checks Rpath parameter files for consistency.  An error message will be produced if one of
#'  the logical checks fails.  Checks include: 
#'  (NOTE: This does not ensure data is correct just that it is in the right places).
#'@import data.table
#'@export
check.rpath.param <- function(filename, parameter = 'model'){
  parameter <- as.data.table(read.csv(filename))
  
}