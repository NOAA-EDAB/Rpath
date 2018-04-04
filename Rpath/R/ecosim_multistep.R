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
rsim.step <- function(Rsim.scenario, Rsim.output, method = 'AB', step.end){
  scene    <- copy(Rsim.scenario)
  full.run <- copy(Rsim.output)
  last <- nrow(full.run$out_BB)
  step.start <- last/12 + 1
  if(method == 'AB'){
    next.run <- Adams_run(scene$params, full.run$end_state,
                          scene$forcing, scene$fishing,
                          scene$stanzas, step.start, step.end,
                          full.run$dyt)
    }
  #Merge runs
  start.month <- last + 1
  end.month   <- step.end * 12
  full.run$out_BB <- rbind(full.run$out_BB, next.run$out_BB[start.month:end.month, ])
  full.run$out_CC <- rbind(full.run$out_CC, next.run$out_CC[start.month:end.month, ])
  full.run$out_Gear_CC <- rbind(full.run$out_Gear_CC, 
                                next.run$out_Gear_CC[start.month:end.month, ])
  full.run$annual_BB <- rbind(full.run$annual_BB,
                              next.run$annual_BB[step.start:step.end, ])
  full.run$annual_CC <- rbind(full.run$annual_CC,
                              next.run$annual_CC[step.start:step.end, ])
  full.run$annual_QB <- rbind(full.run$annual_QB,
                              next.run$annual_QB[step.start:step.end, ])
  full.run$annual_Qlink <- rbind(full.run$annual_Qlink,
                                 next.run$annual_Qlink[step.start:step.end, ])
  full.run$end_state <- next.run$end_state
  
  return(full.run)
}
