#####################################################################################
#'Initial function to step through ecosim
#'
#'Runs rsim.run at intervals and combines the run to one output
#'
#'@family Rpath functions
#'
#'@param Rsim.scenario Rpath object created using rsim.scenario.
#'@param Rsim.output Rpath ecosim run created by the rsim.run() function.
#'@param method Which integration algorithim rsim.run will use (Currently only
#'              supports 'AB')
#'@param step.end The end year for the current interval
#'
#'@return Returns an Rsim.output object.
#'@export
rsim.step <- function(Rsim.scenario, Rsim.output, method = 'AB',year.end){
  scene    <- copy(Rsim.scenario)
  full.run <- copy(Rsim.output)
  
  scene.years <- row.names(Rsim.scenario$fishing$FRATE)
  year.start  <- as.numeric(tail(row.names(Rsim.output$annual_BB),1)) + 1
  step.start  <- which(scene.years==year.start)
  step.end    <- which(scene.years==year.end)
  
  if(method == 'AB'){
    next.run <- Adams_run(Rsim.scenario$params, full.run$end_state,
                          Rsim.scenario$forcing, Rsim.scenario$fishing,
                          Rsim.scenario$stanzas, step.start, step.end,
                          full.run$dyt)
  }
  #Merge runs
  last <- nrow(full.run$out_BB)
  start.month <- last + 1
  end.month   <- step.end * 12
  
  # KYA - Need to subset here and add names, to prevent losing matrix format
  # and labels when only 1 year is run.
  annual_CC    <- next.run$annual_CC[step.start:step.end,,drop=FALSE]
  annual_BB    <- next.run$annual_BB[step.start:step.end,,drop=FALSE]
  annual_QB    <- next.run$annual_QB[step.start:step.end,,drop=FALSE]
  annual_Qlink <- next.run$annual_Qlink[step.start:step.end,,drop=FALSE]
  ylist <- year.start:year.end
  rownames(annual_CC) <-    ylist 
  rownames(annual_BB) <-    ylist 
  rownames(annual_QB) <-    ylist
  rownames(annual_Qlink) <- ylist  
  
  full.run$out_BB <- rbind(full.run$out_BB, next.run$out_BB[start.month:end.month, ])
  full.run$out_CC <- rbind(full.run$out_CC, next.run$out_CC[start.month:end.month, ])
  full.run$out_Gear_CC <- rbind(full.run$out_Gear_CC, 
                                next.run$out_Gear_CC[start.month:end.month, ])
  full.run$annual_BB <- rbind(full.run$annual_BB, annual_BB)
  full.run$annual_CC <- rbind(full.run$annual_CC, annual_CC)
  full.run$annual_QB <- rbind(full.run$annual_QB, annual_QB)
  full.run$annual_Qlink <- rbind(full.run$annual_Qlink, annual_Qlink)
  full.run$end_state <- next.run$end_state
  
  return(full.run)
}
