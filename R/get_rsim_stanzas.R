#'Retrieve Rsim stanza parameters
#'
#'Helper function that will retrieve the stanza parameters that were used in an 
#'Rsim scenario 
#'
#'@family Get functions
#'
#'@inheritParams rsim.run
#'
#'@return Returns an \code{Rsim.stanzas} object.
#'
#'@export
#' 
get.rsim.stanzas<-function(Rsim.scenario){
  return(Rsim.scenario$stanzas)
}