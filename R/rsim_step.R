#'Function to add a time step to an Rsim run
#'
#'Runs \code{rsim.run} and appends the output to a previous \code{Rsim.output}.
#'Currently only works with the Adams-Bashforth numerical integration method (AB).
#'
#'@family Rsim functions
#'
#'@inheritParams rsim.run
#'@inheritParams rsim.diet
#'@param year.end The final year for the step interval
#'
#'@return Returns an \code{Rsim.output} object.
#'
#'@export
#'
rsim.step <- function(Rsim.scenario, Rsim.output, method = 'AB',year.end){
  scene    <- copy(Rsim.scenario)
  full.run <- copy(Rsim.output)
  
  # KYA adds run date and some random salt to ensure uniquieness  
    scene$rundate <- paste(Sys.time(),":salt:",runif(1))    
  
  scene.years <- row.names(Rsim.scenario$fishing$ForcedFRate)
  year.start  <- as.numeric(tail(row.names(Rsim.output$annual_Biomass),1)) + 1
  step.start  <- which(scene.years==year.start)
  step.end    <- which(scene.years==year.end)
  
  if(method == 'AB'){
    next.run <- Adams_run(Rsim.scenario$params, full.run$end_state,
                          Rsim.scenario$forcing, Rsim.scenario$fishing,
                          Rsim.scenario$stanzas, step.start, step.end,
                          full.run$dyt)
  }
  #Merge runs
  last <- nrow(full.run$out_Biomass)
  start.month <- last + 1
  end.month   <- step.end * 12
  
  # KYA - Need to subset here and add names, to prevent losing matrix format
  # and labels when only 1 year is run.
  annual_Catch    <- next.run$annual_Catch[step.start:step.end,,drop=FALSE]
  annual_Biomass    <- next.run$annual_Biomass[step.start:step.end,,drop=FALSE]
  annual_QB    <- next.run$annual_QB[step.start:step.end,,drop=FALSE]
  annual_Qlink <- next.run$annual_Qlink[step.start:step.end,,drop=FALSE]
  ylist <- year.start:year.end
  rownames(annual_Catch) <-    ylist 
  rownames(annual_Biomass) <-    ylist 
  rownames(annual_QB) <-    ylist
  rownames(annual_Qlink) <- ylist  
  
  full.run$out_Biomass <- rbind(full.run$out_Biomass, next.run$out_Biomass[start.month:end.month, ])
  full.run$out_Catch <- rbind(full.run$out_Catch, next.run$out_Catch[start.month:end.month, ])
  full.run$out_Gear_Catch <- rbind(full.run$out_Gear_Catch, 
                                next.run$out_Gear_Catch[start.month:end.month, ])
  full.run$annual_Biomass <- rbind(full.run$annual_Biomass, annual_Biomass)
  full.run$annual_Catch <- rbind(full.run$annual_Catch, annual_Catch)
  full.run$annual_QB <- rbind(full.run$annual_QB, annual_QB)
  full.run$annual_Qlink <- rbind(full.run$annual_Qlink, annual_Qlink)
  full.run$end_state <- next.run$end_state
  
  return(full.run)
}
