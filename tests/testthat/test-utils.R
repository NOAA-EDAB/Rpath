
source("test-constants.R")


randomNumber <- function(seed) {
  set.seed(seed)
  rval <- runif(1,min=-JITTER_AMOUNT_PCT/100,max=JITTER_AMOUNT_PCT/100)[1]
  return(rval)
}


# This is an increment function
inc <- function(value) eval.parent(substitute(value <- value + 1))

readDataFile <- function(filename,fill=TRUE,sep="") {
  readRDS(filename)
  # read.table(filename,fill=fill,sep=sep)
}

writeDataFile <- function(object,filename) {
  saveRDS(object,file=filename)
  # write.table(object,file=filename)
}

printStatsScenario <- function(msg,scenario) {
  print("- - -")
  print(paste0(msg,", forcing$ForcedBio:     ",sum(scenario$forcing$ForcedBio)))
  print(paste0(msg,", forcing$ForcedMigrate: ",sum(scenario$forcing$ForcedMigrate)))
  print(paste0(msg,", start_state$Biomass:   ",sum(scenario$start_state$Biomass)))
}

printStatsSimulation <- function(msg,sim) {
  print("- - -")
  print(paste0(msg,", sim$out_Biomass:    ",sum(sim$out_Biomass)))
  print(paste0(msg,", sim$out_Catch:      ",sum(sim$out_Catch)))
  print(paste0(msg,", sim$out_Gear_Catch: ",sum(sim$out_Gear_Catch)))
}
