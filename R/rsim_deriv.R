#'Calculate the derivatives for a time step 
#'
#'Calculates the derivative for a single time step and saves the output 
#'
#'@inheritParams rsim.run
#'@param sim.year Will inherit from apply functions
#'@param sim.month Will inherit from apply functions
#'@param tstep Sub-monthly time step usually set to 0. 
#'
#'@export
#'
rsim.deriv <- function(Rsim.scenario, sim.year = 0, sim.month = 0, tstep = 0){
  scene <- copy(Rsim.scenario)
  rout <- deriv_vector(scene$params,  scene$start_state, 
                       scene$forcing, scene$fishing,
                       scene$stanzas, sim.year, sim.month, tstep)
  
  rtab <- data.frame(scene$params$spname,rout$DerivT,rout$TotGain,rout$TotLoss, 
                     rout$FoodGain, rout$DetritalGain, rout$FishingGain,      
                     rout$UnAssimLoss,rout$ActiveRespLoss,
                     rout$FoodLoss,rout$MzeroLoss,rout$FishingLoss,rout$DetritalLoss)
  colnames(rtab)<-c("Species","DerivT","TotGain","TotLoss","FoodGain","DetritalGain","FishingGain",
                    "UnAssimLoss","ActiveRespLoss","FoodLoss","MzeroLoss","FishingLoss","DetritalLoss")
  return(rtab)
}