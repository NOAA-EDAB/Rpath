#Functions for Rpath objects
#Print Rpath
#'@import utils
#'@export
print.Rpath <- function(x, rows = NA, morts = F, ...){
  cat(paste("Rpath model:", attr(x, 'eco.name'), "\n"))
  cat(paste("Model Area:",  attr(x, 'eco.area'), "\n"))
  if(max(x$EE, na.rm = T) > 1){
    unbalanced.groups <- x$Group[which(x$EE > 1)]
    cat("     Status: Unbalanced! \nThe following groups have EE > 1:\n")
    print(unbalanced.groups)
    cat("\n")
  } else {
    cat("     Status: Balanced\n")
  }
  if(morts == F){
    removals <- rowSums(x$Landings) + rowSums(x$Discard)
    out <- data.frame(Group    = x$Group,
                      type     = x$type,
                      TL       = x$TL,
                      Biomass  = x$Biomass,
                      PB       = x$PB,
                      QB       = x$QB,
                      EE       = x$EE,
                      GE       = x$GE,
                      Removals = removals)
  }
  if(morts == T){
    ngroup <- x$NUM_LIVING + x$NUM_DEAD
    out <- data.frame(Group    = x$Group[1:ngroup],
                      type     = x$type [1:ngroup],
                      PB       = x$PB   [1:ngroup])
    #Calculate M0
    M0  <- c(x$PB[1:x$NUM_LIVING] * (1 - x$EE[1:x$NUM_LIVING]), 
             x$EE[(x$NUM_LIVING + 1):ngroup])
    out <- cbind(out, M0)
    #Calculate F mortality
    totcatch <- x$Landings + x$Discards
    Fmort    <- as.data.frame(totcatch / x$Biomass[row(as.matrix(totcatch))])
    setnames(Fmort, x$Group[(ngroup +1):x$NUM_GROUPS], 
                    paste('F.', x$Group[(ngroup +1):x$NUM_GROUPS], sep = ''))
    out  <- cbind(out, Fmort[1:ngroup, ])
    #Calculate M2
    bio  <- x$Biomass[1:x$NUM_LIVING]
    BQB  <- bio * x$QB[1:x$NUM_LIVING]
    diet <- as.data.frame(x$DC)
    nodetrdiet <- diet[1:x$NUM_LIVING, ]
    detrdiet   <- diet[(x$NUM_LIVING +1):ngroup, ]
    newcons    <- nodetrdiet * BQB[col(as.matrix(nodetrdiet))]
    predM      <- newcons / bio[row(as.matrix(newcons))]
    detcons    <- detrdiet * BQB[col(as.matrix(detrdiet))]
    predM      <- rbind(predM, detcons)
    setnames(predM, x$Group[1:x$NUM_LIVING], 
             paste('M2.', x$Group[1:x$NUM_LIVING], sep = ''))
    out <- cbind(out, predM)
  }
  if(is.na(rows)) print(out, nrows = Inf) else head(out, n = rows)
}

#Print Rpath.sim
#'@import utils
#'@export
print.Rsim.output <- function(x, rows = NA, ...){
  cat(paste("Rpath sim results:", attr(x, 'eco.name'),"\n"))
  if(x$crash_year > 0) cat(paste("Run crashed at", x$crash_year, "\n", sep = ''))
  
  gear.zero <- rep(0, x$params$NUM_GEARS)
  start_Catch <- c(x$out_Catch[2, ], gear.zero)
  end_Catch   <- c(x$out_Catch[nrow(x$out_Catch) - 1, ], gear.zero)
  out <- data.frame(Group      = x$params$spname,
                    StartBio   = x$start_state$Biomass,
                    EndBio     = x$end_state$Biomass,
                    BioES      = x$end_state$Biomass / 
                                 x$start_state$Biomass,
                    StartCatch = start_Catch * 12,
                    EndCatch   = end_Catch * 12,
                    CatchES    = (end_Catch * 12) / (start_Catch * 12))
  
  if(is.na(rows)) print(out, nrows = Inf) else head(out, n = rows)
}

#Print Rsim.scenario
#'@import utils
#'@export
print.Rsim.scenario <- function(x, ...){
  cat(paste("Rpath scenario for", attr(x, 'eco.name'), "\n\n"))
  cat("$params contains the parameters from rpath
$forcing contains the forcing parameters
$fishing contains the fishing parameters
$state contains the initial state parameters \n
Use adjust functions to modify
$forcing or $fishing to alter scenario run")
}

#Print Rsim.scenario
#'@import utils
#'@export
print.Rsim.params <- function(x, ...){
  cat(paste("Rsim parameters for", attr(x, 'eco.name'), "\n\n"))
  out <- data.frame(NumGroups   = x$NUM_GROUPS,
                    NumLiving   = x$NUM_LIVING,
                    NumDetritus = x$NUM_DEAD,
                    NumFleets   = x$NUM_GEARS)
  print(out)
  cat("\n$params also includes:\n")
  print(names(x))
}


#Summary for Rpath
#'@export
summary.Rpath <- function(object, ...){
  cat(paste("Rpath model:", attr(object, 'eco.name'),"\n"))
  if(max(object$EE, na.rm = T) > 1){
    unbalanced.groups <- object$Group[which(object$EE > 1)]
    cat("     Status: Unbalanced! \nThe following groups have EE > 1:\n")
    print(unbalanced.groups)
    cat("\n")
  } else {
    cat("     Status: Balanced\n")
  }
  cat("\nSummary Statistics:\n")
  totbiomass <- sum(object$Biomass[which(object$type == 0)], na.rm = T)
  totland    <- sum(object$Landings, na.rm = T)
  out <- data.frame(NumGroups   = object$NUM_GROUPS,
                    NumLiving   = object$NUM_LIVING,
                    NumDetritus = object$NUM_DEAD,
                    NumFleets   = object$NUM_GEARS,
                    TotBiomass   = totbiomass,
                    TotLandings  = totland)
  print(out)
  cat("\nRpath model also includes:\n")
  print(names(object))
}

#Summary for Rpath.sim
#'@export
summary.Rsim.output <- function(object, ...){
  cat(paste("Rsim parameters for:", attr(object, 'eco.name'),"\n"))
  if(object$crash_year > 0) cat(paste("Run crashed at", object$crash_year, "\n", sep = ''))
  cat("\nSummary Statistics:\n")
  totbiomass.start <- sum(object$out_Biomass[1, ], na.rm = T)
  totbiomass.end   <- sum(object$out_Biomass[nrow(object$out_Biomass), ], na.rm = T)
  totcatch.start   <- sum(object$out_Catch[1, ] * 12, na.rm = T)
  totcatch.end     <- sum(object$out_Catch[nrow(object$out_Catch) - 1, ] * 12, na.rm = T)
  out <- data.frame(NumGroups      = length(object$params$spname) - 1,
                    NumLiving      = object$params$NUM_LIVING,
                    NumDetritus    = object$params$NUM_DEAD,
                    NumFleets      = object$params$NUM_GEARS,
                    TotBiomassStart = totbiomass.start,
                    TotBiomassEnd   = totbiomass.end,
                    TotCatchStart   = totcatch.start,
                    TotCatchEnd     = totcatch.end)
  print(out)
  cat("\nRpath sim also includes:\n")
  print(names(object))
}



