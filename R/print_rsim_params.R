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