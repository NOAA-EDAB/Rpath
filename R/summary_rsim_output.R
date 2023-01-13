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