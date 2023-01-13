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