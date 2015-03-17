## Rpath functions for creating and checking parameter files  

#'Create parameter file shells
#'
#'Creates a shell of the parameter files that can then be filled out in R, Excel, or
#'another spreadsheet program.
#'
#'@family Rpath functions
#'
#'@param group Vector of group names.
#'@param type Numeric vector of group type. Living = 0, Producer = 1, Detritus = 2,
#'  Fleet = 3.
#'@param file File name for output parameters. If NA the files will not be written
#'  to comma delimited files.
#'@param juvenile Logical value whether the juvenile file is being created.
#'
#'@return Outputs the model and diet parameter shells.  If a file name is provided the 
#'  extentions "_mod" and "_diet" will be appended to the file name.  If juvenile equals 
#'  True, the "_juv" file will be created instead.
#'@import data.table
#'@export
create.rpath.param <- function(file = NA, group = NA, type = NA, juvenile = F){
  if(juvenile == F){
    #Base model
    model <- data.table(Group    = group, 
                        Type     = type, 
                        Biomass  = 0,
                        PB       = 0,
                        QB       = 0,
                        EE       = 0,
                        ProdCons = 0,
                        BioAcc   = 0,
                        Unassim  = 0,
                        DetInput = 0)
    
    #Add detritial groups
    if(length(which(type == 2)) > 0){
      for(i in 1:length(which(type == 2))){
        model[, V1 := 0]
        setnames(model, "V1", model[Type == 2, Group][i])
      }
    }
    
    #Add fleets twice - Landings and Discards
    if(length(which(type == 3)) > 0){
      for(i in 1:length(which(type == 3))){
        model[, V1 := 0]
        setnames(model, "V1", model[Type == 3, Group][i])
      }
      for(i in 1:length(which(type == 3))){
        model[, V1 := 0]
        setnames(model, "V1", model[Type == 3, Group][i])
      }
    }
    
    #Diet matrix
    diet <- data.table(Group = group)
    if(!is.na(group)){
      for(i in 1:length(which(type < 2))){
        diet[, V1 := 0]
        setnames(diet, "V1", group[i])
      }
    }
  }
}