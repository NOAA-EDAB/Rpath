#'Retrieve Rsim fishing forcing parameters
#'
#'Helper function that will retrieve the fishing forcing parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.fishing} object.
#'
#'@export
#'  
get.rsim.fishing<-function(Rsim.scenario){
  return(Rsim.scenario$fishing)
}