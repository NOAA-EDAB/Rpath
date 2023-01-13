#'Generate Rsim state matrix
#'
#'Creates a matrix of state variables used by Rsim.
#'
#'@inheritParams rsim.fishing
#'
#'@export
#'
rsim.state <- function(params){
  state  <- list(Biomass = params$B_BaseRef, 
                 N       = rep(0, params$NUM_GROUPS + 1),
                 Ftime   = rep(1, length(params$B_BaseRef)))
  class(state) <- "Rsim.state"
  return(state)
}