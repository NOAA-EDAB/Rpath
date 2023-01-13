#'Retrieve Rsim starting state values
#'
#'Helper function that will retrieve the starting state values that were used 
#'in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.state} object.
#'
#'@export
#'   
get.rsim.start_state<-function(Rsim.scenario){
  return(Rsim.scenario$start_state)
}