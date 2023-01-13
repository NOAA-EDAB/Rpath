#'Check Rpath parameter files
#'
#'Logical check that the parameter files are filled out correctly, i.e. data is entered where it is
#'expected.
#'
#'@family Rpath functions
#'
#'@inheritParams rpath
#'
#'@return Checks Rpath parameter files for consistency.  An error message will be produced if one of
#'  the logical checks fails.  Checks include: 
#'  (NOTE: This does not ensure data is correct just that it is in the right places).
#'@import data.table
#'@export
#'
check.rpath.params <- function(Rpath.params){
  #Need to define variables to eliminate check() note about no visible binding
  Type <- Group <- Biomass <- EE <- PB <- QB <- ProdCons <- BioAcc <- Unassim <- DetInput <- NULL
  
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
    cat('Rpath parameter file is functional. \n')
  } else {
    cat('Rpath parameter file needs attention! \n')
  }  
}