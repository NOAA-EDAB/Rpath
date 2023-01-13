#'Output consumption by a group
#'
#'Creates a matrix of consumption of prey by a particular predator.
#'
#'@param Rsim.output R object containing the output from \code{rsim.run}.
#'@param group Group from the \code{Rpath} model that is of interest
#'
#'@export
#'
rsim.diet <- function(Rsim.output, group){
  pmat <- Rsim.output$annual_Qlink[, Rsim.output$pred == group]
  colnames(pmat) <- Rsim.output$prey[Rsim.output$pred == group]
  return(pmat)
}  