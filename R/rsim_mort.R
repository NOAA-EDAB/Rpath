#'Output mortality on a group
#'
#'Creates a matrix of mortality by predators on a particular prey.
#'
#'@inheritParams rsim.diet 
#'
#'@export
#'
rsim.mort <- function(Rsim.output, group){
  pmat <- Rsim.output$annual_Qlink[, Rsim.output$prey == group]
  colnames(pmat) <- Rsim.output$pred[Rsim.output$prey == group]
  return(pmat)
} 