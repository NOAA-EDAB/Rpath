#####################################################################################
#'Initial function to step through ecosim
#'
#'Runs rsim.run at intervals and combines the run to one output
#'
#'@family Rpath functions
#'
#'@param Rsim.scenario Rpath object created using rsim.scenario.
#'@param method Which integration algorithim rsim.run will use (Currently only
#'              supports 'AB')
#'@param max.year maximum length of the simulation
#'@param interval number of years inbetween simulations
#'
#'@return Returns an Rsim.output object.
#'@export
rsim.step <- function(Rsim.scenario, Rsim.run, method = 'AB', start, end){
  full.run <- Rsim.run
  if(method == 'AB'){
    next.run <- Adams_run(Rsim.scenario$params, Rsim.run$end_state,
                          Rsim.scenario$forcing, Rsim.scenario$fishing,
                          Rsim.scenario$stanzas, start - 1, end,
                          Rsim.run$dyt)
    }
  #Merge runs
  first <- (start - 1) * 12 + 1
  last  <- end * 12 + 1
  full.run$out_BB <- rbind(full.run$out_BB, next.run$out_BB[(first + 1):last, ])
  full.run$out_CC <- rbind(full.run$out_CC, next.run$out_CC[(first + 1):last, ])
  full.run$out_Gear_CC <- rbind(full.run$out_Gear_CC, 
                                next.run$out_Gear_CC[(first + 1):last, ])
  full.run$annual_BB <- rbind(full.run$annual_BB[1:start - 1, ],
                              next.run$annual_BB[start:end, ])
  full.run$annual_CC <- rbind(full.run$annual_CC[1:start - 1, ],
                              next.run$annual_CC[start:end, ])
  full.run$annual_QB <- rbind(full.run$annual_QB[1:start - 1, ],
                              next.run$annual_QB[start:end, ])
  full.run$annual_Qlink <- rbind(full.run$annual_Qlink[1:start - 1, ],
                                 next.run$annual_Qlink[start:end, ])
  full.run$end_state <- next.run$end_state
  
  return(full.run)
}
