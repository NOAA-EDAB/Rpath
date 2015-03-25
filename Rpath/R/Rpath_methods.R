#Functions for Rpath objects
#Print Rpath
#'@export
print.Rpath <- function(x, rows = NA, morts = F, ...){
  cat(paste("Rpath model:", attr(x, 'eco.name'),"\n"))
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
#'@export
print.Rpath.sim <- function(x, rows = NA, ...){
  cat(paste("Rpath sim results:", attr(x, 'eco.name'),"\n"))
  if(x$CRASH_YEAR > 0) cat(paste("Run crashed at", x$CRASH_YEAR, "\n", sep = ''))
  out <- c()
  for(i in 1:(x$NUM_LIVING + x$NUM_DEAD)){
    sp.out <- data.frame(Group      = x$spname[i],
                         StartBio   = x$out_BB[1, i],
                         EndBio     = x$out_BB[nrow(x$out_CC), i],
                         StartCatch = x$out_CC[1, i],
                         EndCatch   = x$out_CC[nrow(x$out_CC) - 1, i])
    out <- rbind(out, sp.out)
  }
  if(is.na(rows)) print(out, nrows = Inf) else head(out, n = rows)
}
  
#Summary for Rpath
#'@export
summary.Rpath <- function(x, ...){
  cat(paste("Rpath model:", attr(x, 'eco.name'),"\n"))
  if(max(x$EE, na.rm = T) > 1){
    unbalanced.groups <- x$Group[which(x$EE > 1)]
    cat("     Status: Unbalanced! \nThe following groups have EE > 1:\n")
    print(unbalanced.groups)
    cat("\n")
  } else {
    cat("     Status: Balanced\n")
  }
  cat("\nSummary Statistics:\n")
  totbiomass <- sum(x$BB[which(x$type == 0)],    na.rm = T)
  totcatch   <- sum(x$Catch, na.rm = T)
  out <- data.frame(Num.Groups   = x$NUM_GROUPS,
                    Num.Living   = x$NUM_LIVING,
                    Num.Detritus = x$NUM_DEAD,
                    Num.Fleets   = x$NUM_GEARS,
                    TotBiomass   = totbiomass,
                    TotCatch     = totcatch)
  print(out)
  cat("\nRpath model also includes:\n")
  print(names(x))
}

#Summary for Rpath.sim
#'@export
summary.Rpath.sim <- function(x, ...){
  cat(paste("Rpath sim results:", attr(x, 'eco.name'),"\n"))
  if(x$CRASH_YEAR > 0) cat(paste("Run crashed at", x$CRASH_YEAR, "\n", sep = ''))
  cat("\nSummary Statistics:\n")
  totbiomass.start <- sum(x$out_BB[1, ],                  na.rm = T)
  totbiomass.end   <- sum(x$out_BB[nrow(x$out_BB), ],     na.rm = T)
  totcatch.start   <- sum(x$out_CC[1, ],                  na.rm = T)
  totcatch.end     <- sum(x$out_CC[nrow(x$out_CC) - 1, ], na.rm = T)
  out <- data.frame(Num.Groups        = x$NUM_GROUPS,
                    Num.Living      = x$NUM_LIVING,
                    Num.Detritus    = x$NUM_DEAD,
                    Num.Fleets      = x$NUM_GEARS,
                    TotBiomassStart = totbiomass.start,
                    TotBiomassEnd   = totbiomass.end,
                    TotCatchStart   = totcatch.start,
                    TotCatchEnd     = totcatch.end)
  print(out)
  cat("\nRpath sim also includes:\n")
  print(names(x))
}



