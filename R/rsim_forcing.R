#'Generate Rsim forcing matrix
#'
#'Creates a matrix for forcing functions not related to fishing (\emph{ForcedPrey}
#'\emph{ForcedMort}, \emph{ForcedRecs}, \emph{ForcedSearch}, \emph{ForcedMigrate}, 
#'\emph{ForcedBio}).
#'
#'@inheritParams rsim.fishing
#'
#'@export
#'
rsim.forcing <- function(params, years){
  # Monthly index defaulting to to 1.0, for environmental forcing list
  
  nyrs <- length(years)
  MF <- (matrix(1.0, nyrs * 12, params$NUM_GROUPS + 1))
  year.m <- c()
  for(iyear in 1:nyrs){
    iyear.m <- paste0(years[iyear], '.', 1:12)
    year.m  <- c(year.m, iyear.m)
  }
  rownames(MF) <- year.m
  colnames(MF) <- params$spname
  forcing <- list(ForcedPrey    = MF, 
                  ForcedMort    = MF, 
                  ForcedRecs    = MF, 
                  ForcedSearch  = MF,
                  ForcedActresp = MF,
                  ForcedMigrate = MF * 0,
                  ForcedBio     = MF * -1)
  
  class(forcing) <- "Rsim.forcing"
  return (forcing)
}