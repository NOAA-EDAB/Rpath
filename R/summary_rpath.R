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