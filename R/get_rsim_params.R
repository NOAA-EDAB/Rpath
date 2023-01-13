#'Retrieve Rsim scenario parameters
#'
#'Helper function that will retrieve parameters that were used in an Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.params} object.
#'
#'@export
#'    
get.rsim.params<-function(Rsim.scenario){
  return(Rsim.scenario$params)
}