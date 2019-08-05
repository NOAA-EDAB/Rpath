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
    removals <- rowSums(x$Catch) + rowSums(x$Discard)
    out <- data.frame(Group    = x$Group,
                      type     = x$type,
                      TL       = x$TL,
                      Biomass  = x$BB,
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
    totcatch <- x$Catch + x$Discards
    Fmort    <- as.data.frame(totcatch / x$BB[row(as.matrix(totcatch))])
    setnames(Fmort, paste('V',  1:x$NUM_GEARS,                     sep = ''), 
                    paste('F.', x$Group[(ngroup +1):x$NUM_GROUPS], sep = ''))
    out  <- cbind(out, Fmort[1:ngroup, ])
    #Calculate M2
    bio  <- x$BB[1:x$NUM_LIVING]
    BQB  <- bio * x$QB[1:x$NUM_LIVING]
    diet <- as.data.frame(x$DC)
    nodetrdiet <- diet[1:x$NUM_LIVING, ]
    detrdiet   <- diet[(x$NUM_LIVING +1):ngroup, ]
    newcons    <- nodetrdiet * BQB[col(as.matrix(nodetrdiet))]
    predM      <- newcons / bio[row(as.matrix(newcons))]
    detcons    <- detrdiet * BQB[col(as.matrix(detrdiet))]
    predM      <- rbind(predM, detcons)
    setnames(predM, paste('V',  1:x$NUM_LIVING,    sep = ''), 
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
  start_CC <- c(x$out_CC[2, ], gear.zero)
  end_CC   <- c(x$out_CC[nrow(x$out_CC) - 1, ], gear.zero)
  out <- data.frame(Group      = x$params$spname,
                    StartBio   = x$start_state$BB,
                    EndBio     = x$end_state$BB,
                    BioES      = x$end_state$BB / 
                                 x$start_state$BB,
                    StartCatch = start_CC * 12,
                    EndCatch   = end_CC * 12,
                    CatchES    = (end_CC * 12) / (start_CC * 12))
  
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
Modify $forcing or $fishing to alter scenario run")
}

#Print Rsim.scenario
#'@import utils
#'@export
print.Rsim.params <- function(x, ...){
  cat(paste("Rsim parameters for", attr(x, 'eco.name'), "\n\n"))
  out <- data.frame(Num.Groups   = x$NUM_GROUPS,
                    Num.Living   = x$NUM_LIVING,
                    Num.Detritus = x$NUM_DEAD,
                    Num.Fleets   = x$NUM_GEARS)
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
  totbiomass <- sum(object$BB[which(object$type == 0)],    na.rm = T)
  totcatch   <- sum(object$Catch, na.rm = T)
  out <- data.frame(Num.Groups   = object$NUM_GROUPS,
                    Num.Living   = object$NUM_LIVING,
                    Num.Detritus = object$NUM_DEAD,
                    Num.Fleets   = object$NUM_GEARS,
                    TotBiomass   = totbiomass,
                    TotCatch     = totcatch)
  print(out)
  cat("\nRpath model also includes:\n")
  print(names(object))
}

#Summary for Rpath.sim
#'@export
summary.Rsim.output <- function(object, ...){
  cat(paste("Rsim parameters for:", attr(object, 'eco.name'),"\n"))
  if(object$CRASH_YEAR > 0) cat(paste("Run crashed at", object$CRASH_YEAR, "\n", sep = ''))
  cat("\nSummary Statistics:\n")
  totbiomass.start <- sum(object$out_BB[1, ],                       na.rm = T)
  totbiomass.end   <- sum(object$out_BB[nrow(object$out_BB), ],          na.rm = T)
  totcatch.start   <- sum(object$out_CC[1, ] * 12,                  na.rm = T)
  totcatch.end     <- sum(object$out_CC[nrow(object$out_CC) - 1, ] * 12, na.rm = T)
  out <- data.frame(Num.Groups        = object$NUM_GROUPS,
                    Num.Living      = object$NUM_LIVING,
                    Num.Detritus    = object$NUM_DEAD,
                    Num.Fleets      = object$NUM_GEARS,
                    TotBiomassStart = totbiomass.start,
                    TotBiomassEnd   = totbiomass.end,
                    TotCatchStart   = totcatch.start,
                    TotCatchEnd     = totcatch.end)
  print(out)
  cat("\nRpath sim also includes:\n")
  print(names(object))
}



