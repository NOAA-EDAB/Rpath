#Print Rsim.scenario
#'@import utils
#'@export
print.Rsim.scenario <- function(x, ...){
  cat(paste("Rsim scenario for", attr(x, 'eco.name'), "\n\n"))
  cat("$params contains the parameters from rpath
$forcing contains the forcing parameters
$fishing contains the fishing parameters
$state contains the initial state parameters \n
Use adjust functions to modify
$forcing or $fishing to alter scenario run")
}