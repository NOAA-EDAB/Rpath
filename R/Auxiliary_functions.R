#'Quantify mixed trophic impacts
#'
#'Uses a static \code{Rpath} model and creates a matrix of mixed trophic impacts.
#'
#'@family Rpath functions
#'
#'@param Rpath Rpath R object containing a static \code{Rpath} model.
#'@param Rpath.params R object containing the \code{Rpath} model parameters.  This is generated
#'  either by the \code{\link{create.rpath.params}()} or \code{\link{read.rpath.params}()} functions.
#'
#'@return Returns a matrix of mixed trophic impacts.
#'
#'@export 
MTI <- function(Rpath, Rpath.params){
  
  # Useful indices
  N_BIO <- Rpath$NUM_LIVING+Rpath$NUM_DEAD
  
  # Make empty matrices
  AA <- matrix(0, Rpath$NUM_GROUPS, Rpath$NUM_GROUPS)
  colnames(AA) <- Rpath$Group
  rownames(AA) <- Rpath$Group
  DCC <- AA
  
  # Diet matrix removing import row, group column, adding detrital columns
  DCC[1:N_BIO,1:N_BIO] <- cbind(as.matrix(Rpath.params$diet[,-"Group"])[1:(dim(Rpath.params$diet)[1]-1),],
                                matrix(0, nrow=dim(Rpath.params$diet)[1]-1, ncol=Rpath$NUM_DEAD))
  
  # Landings only not discards in DCC (discards never included at all in this??)
  DCC[,(N_BIO+1):Rpath$NUM_GROUPS] <- Rpath$Landings 
  
  # Normalizing and cleaning DCC (diet composition and landings)
  DCC[is.na(DCC)] <- 0
  DCC <- t(t(DCC)/colSums(DCC))
  DCC[is.nan(DCC) | is.infinite(DCC)] <- 0
  
  # Living total flows (and flows FROM detritus via diets, but not TO detritus)
  for (i in 1:N_BIO){
    for (j in 1:Rpath$NUM_LIVING){
      AA[i,j] <- Rpath$Biomass[j] * Rpath$QB[j] * DCC[i,j]
    }
  }
  
  # These lines not needed because flow into detritus in AA matrix is 0'd out below?
  # Keeping in comments for reference.
  #  detflow <- (Rpath$PB*Rpath$Biomass*(1-Rpath$EE)) + Rpath$QB*Rpath$Biomass*Rpath$Unassim
  #  detmat <- detflow * Rpath$DetFate
  #  AA[,(Rpath$NUM_LIVING+1):N_BIO] <- detmat
  # Should Discards be here?  Seems to work but check later?  
  AA[,(N_BIO+1):Rpath$NUM_GROUPS] <- Rpath$Landings + Rpath$Discards 
  fishing_det <- colSums(Rpath$Discards) * Rpath$DetFate[(N_BIO+1):Rpath$NUM_GROUPS,]
  AA[(N_BIO+1):Rpath$NUM_GROUPS, (Rpath$NUM_LIVING+1):N_BIO] <- fishing_det
  
  # "Detritus is not considered a predator" - remove flows to detritus before normalizing
  AA[,(Rpath$NUM_LIVING+1):N_BIO] <- 0
  
  # Normalize FC and and then flip so that predators are in columns    
  FC <- AA/rowSums(AA)
  FC[is.nan(FC) | is.infinite(FC)] <- 0
  FC_T <- t(FC)
  
  # Direct impacts matrix
  MTI <- (DCC - FC_T)
  
  # From Ulanowicz paper, find inverse of ([Identity] - MTI)
  MTI_diag <- diag(Rpath$NUM_GROUPS) - MTI
  
  # Invert the matrix    
  MTI_inv <- MASS::ginv((MTI_diag))
  
  # Now subtract 1 off the diagonal matrix  
  diag(MTI_inv) <- diag(MTI_inv) - 1
  
  # The final(?) piece, hidden in Villy's '03 code - 0 out the detritial Identity 
  for (i in (Rpath$NUM_LIVING+1):N_BIO){
    MTI_inv[i,i] <- 0
  }
  
  colnames(MTI_inv) <- Rpath$Group
  rownames(MTI_inv) <- Rpath$Group
  
  return(MTI_inv)  
}



# MTI <- function(Rpath, Rpath.params, increase = T){
#   #Need to define variables to eliminate check() note about no visible binding
#   Group <- Type <- V1 <- NULL
#   x <- copy(Rpath.params)
#   y <- copy(Rpath)
#   
#   #MTI is a matrix where MTIij = DCij - FCji
#   #Get number of fleets and detrital groups
#   nfleets       <- length(which(x$model$Type == 3))
#   ndetritus     <- length(which(x$model$Type == 2))
#   ngroups       <- length(x$model$Group)
#   fleetnames    <- x$model$Group[which(x$model$Type == 3)]
#   detritusnames <- x$model$Group[which(x$model$Type == 2)]
#   livingnames   <- x$model$Group[which(x$model$Type < 2)]
#   allnames      <- x$model$Group
#   
#   #Set up DCij including detritus and fleet
#   DC <- as.data.table(x$diet)
#   
#   #Remove import and re-standardize remaining DC
#   DC <- DC[Group != 'Import']
#   DC.types  <- x$model[Type < 2, Type]
#   DC.colsum <- DC[, colSums(.SD, na.rm = T), .SDcols = livingnames]
#   DC.check <- round(DC.colsum + DC.types, 3)
#   for(ipred in 1:length(livingnames)){
#     if(DC.check[ipred] != 1){
#       DC[, ipred + 1 := .SD/DC.colsum[ipred], .SDcols = livingnames[ipred]]
#     }
#   }
#   
#   #Add fleet "prey" rows
#   fleetrows <- as.data.table(matrix(rep(0, nfleets * ncol(DC)), nfleets, 
#                                     ncol(DC)))
#   fleetrows[, V1 := fleetnames]
#   DC <- rbindlist(list(DC, fleetrows), use.names = F)
#   
#   #Add detritus and fleet "pred" columns
#   detcols <- as.data.table(matrix(rep(0, ndetritus * ngroups), ngroups, ndetritus))
#   setnames(detcols, paste0('V', seq_along(detritusnames)), detritusnames)
#   DC <- cbind(DC, detcols)
#   
#   #Calculate proportion of catch for fishing "DC"
#   totcatch <- y$Landings + y$Discards
#   totcatch.sum <- colSums(totcatch)
#   for(ifleet in 1:ncol(totcatch)){
#     fleet.prop <- as.data.table(totcatch[, ifleet] / totcatch.sum[ifleet])
#     setnames(fleet.prop, 'V1', fleetnames[ifleet])
#     if(ifleet == 1){
#       fleetcols <- fleet.prop
#     }else{
#       fleetcols <- cbind(fleetcols, fleet.prop)
#     }
#   }
#   DC <- cbind(DC, fleetcols)
#   
#   #Ensure column order is correct
#   DCi.order <- c('Group', DC$Group)
#   setcolorder(DC, DCi.order)
#   
#   #Fix NAs
#   DC[is.na(DC)] <- 0
#   
#   #Remove Group names
#   DC[, Group := NULL]
#   
#   #FCji - Proportion of predation on j from predator i
#   #Calculate consumption
#   bio  <- y$BB
#   BQB  <- bio * y$QB
#   Tij <- DC * BQB[col(as.matrix(DC))]
#   
#   #Add fishery removals
#   Tij <- Tij[, which(!names(Tij) %in% fleetnames), with = F]
#   Tij <- cbind(Tij, totcatch)
#   
#   #Calculate net production
#   Tim <- rowSums(Tij)
#   
#   #Calculate fraction of i's production consumed by pred j
#   FCij <- c()
#   for(iprey in 1:nrow(Tij)){
#     fij <- Tij[iprey, ] / Tim[iprey]
#     FCij <- rbind(FCij, fij)
#   }
#   FCji <- t(FCij)
#   FCji[is.na(FCji)] <- 0
#  
#   #Merge pred and prey
#   net.impact <- as.matrix(DC - FCji)
#   
#   #Create Identity Matrix
#   identity.matrix <- diag(ncol(net.impact))
#   
#   #Calculate all Mixed Trophic Impacts
#   MTI <- MASS::ginv(identity.matrix - net.impact) - identity.matrix
#   
#   
#   return(MTI)
#   
#   
# }
