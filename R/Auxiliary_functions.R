#'Quantify mixed trophic impacts
#'
#'Uses a static \code{Rpath} model and creates a matrix of mixed trophic impacts.
#'
#'@family Rpath functions
#'
#'@param Rpath Rpath R object containing a static \code{Rpath} model.
#'@param Rpath.params R object containing the \code{Rpath} model parameters.  This is generated
#'  either by the \code{\link{create.rpath.params}()} or \code{\link{read.rpath.params}()} functions.
#'@param increase Logical value indicating whether a marginal increase is applied.
#'
#'@return Returns a matrix of mixed trophic impacts.
#'
#'@export 
MTI <- function(Rpath, Rpath.params, increase = T){
  #Need to define variables to eliminate check() note about no visible binding
  Group <- Type <- V1 <- NULL
  x <- copy(Rpath.params)
  y <- copy(Rpath)
  
  #MTI is a matrix where MTIij = DCij - FCji
  #Get number of fleets and detrital groups
  nfleets       <- length(which(x$model$Type == 3))
  ndetritus     <- length(which(x$model$Type == 2))
  ngroups       <- length(x$model$Group)
  fleetnames    <- x$model$Group[which(x$model$Type == 3)]
  detritusnames <- x$model$Group[which(x$model$Type == 2)]
  livingnames   <- x$model$Group[which(x$model$Type < 2)]
  allnames      <- x$model$Group
  
  #Set up DCij including detritus and fleet
  DC <- as.data.table(x$diet)
  
  #Remove import and re-standardize remaining DC
  DC <- DC[Group != 'Import']
  DC.types  <- x$model[Type < 2, Type]
  DC.colsum <- DC[, colSums(.SD, na.rm = T), .SDcols = livingnames]
  DC.check <- round(DC.colsum + DC.types, 3)
  for(ipred in 1:length(livingnames)){
    if(DC.check[ipred] != 1){
      DC[, ipred + 1 := .SD/DC.colsum[ipred], .SDcols = livingnames[ipred]]
    }
  }
  
  #Add fleet "prey" rows
  fleetrows <- as.data.table(matrix(rep(0, nfleets * ncol(DC)), nfleets, 
                                    ncol(DC)))
  fleetrows[, V1 := fleetnames]
  DC <- rbindlist(list(DC, fleetrows), use.names = F)
  
  #Add detritus and fleet "pred" columns
  detcols <- as.data.table(matrix(rep(0, ndetritus * ngroups), ngroups, ndetritus))
  setnames(detcols, paste0('V', seq_along(detritusnames)), detritusnames)
  DC <- cbind(DC, detcols)
  
  #Calculate proportion of catch for fishing "DC"
  totcatch <- y$Landings + y$Discards
  totcatch.sum <- colSums(totcatch)
  for(ifleet in 1:ncol(totcatch)){
    fleet.prop <- as.data.table(totcatch[, ifleet] / totcatch.sum[ifleet])
    setnames(fleet.prop, 'V1', fleetnames[ifleet])
    if(ifleet == 1){
      fleetcols <- fleet.prop
    }else{
      fleetcols <- cbind(fleetcols, fleet.prop)
    }
  }
  DC <- cbind(DC, fleetcols)
  
  #Ensure column order is correct
  DCi.order <- c('Group', DC$Group)
  setcolorder(DC, DCi.order)
  
  #Fix NAs
  DC[is.na(DC)] <- 0
  
  #Remove Group names
  DC[, Group := NULL]
  
  #FCji - Proportion of predation on j from predator i
  #Calculate consumption
  bio  <- y$BB
  BQB  <- bio * y$QB
  Tij <- DC * BQB[col(as.matrix(DC))]
  
  #Add fishery removals
  Tij <- Tij[, which(!names(Tij) %in% fleetnames), with = F]
  Tij <- cbind(Tij, totcatch)
  
  #Calculate net production
  Tim <- rowSums(Tij)
  
  #Calculate fraction of i's production consumed by pred j
  FCij <- c()
  for(iprey in 1:nrow(Tij)){
    fij <- Tij[iprey, ] / Tim[iprey]
    FCij <- rbind(FCij, fij)
  }
  FCji <- t(FCij)
  FCji[is.na(FCji)] <- 0
 
  #Merge pred and prey
  net.impact <- as.matrix(DC - FCji)
  
  #Create Identity Matrix
  identity.matrix <- diag(ncol(net.impact))
  
  #Calculate all Mixed Trophic Impacts
  MTI <- MASS::ginv(identity.matrix - net.impact) - identity.matrix
  
  
  return(MTI)
  
  
}
