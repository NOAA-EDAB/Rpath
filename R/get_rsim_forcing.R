#'Retrieve Rsim forcing parameters
#'
#'Helper function that will retrieve the forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.forcing} object.
#'
#'@export
#'   
get.rsim.forcing<-function(Rsim.scenario){
  return(Rsim.scenario$forcing)
}