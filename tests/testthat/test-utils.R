source("test-constants.R")

#' Random number generator
#' 
#' This function returns a random number from a uniform distribution between
#' a min and max value.
#'
#' @param pctToJitter : the percent to jitter (value=0.5 means 50%)
#' @param positiveOnly : boolean that signifies random number returned is limited to non-negative numbers
#'
#' @return Returns the random value
#' 
randomNumber <- function(pctToJitter,positiveOnly) {
  minJitter <- -pctToJitter
  maxJitter <- pctToJitter
  if (positiveOnly) {
    minJitter <- 0
  }
  rval <- runif(1,min=minJitter,max=maxJitter) # [1]
  # print(paste0("rval: ",rval))
  return(rval)
}

#' An increment function
#' 
#' This function increments a number by 1
#'
#' @param value : the number to increment by 1
#'
#' @return Returns the incremented value
#' 
inc <- function(value) eval.parent(substitute(value <- value + 1))

#' Read the passed data file
#' 
#' This function reads a data file (either a data table or an rds file)
#'
#' @param filename : the file to read in
#' @param fill : boolean to fill blanks with 0's
#' @param sep : the separator used in the file
#'
#' @return Returns the read in file
#' 
readDataFile <- function(filename,fill=TRUE,sep="") {
  readRDS(filename)
  # read.table(filename,fill=fill,sep=sep)
}

#' Write out the passed data file
#' 
#' This function writes out a data file (either a data table or an rds file)
#'
#' @param object : the object to write out (i.e., serialize)
#' @param filename : the name of the output file to write
#'
#' @return None
#' 
writeDataFile <- function(object,filename) {
  saveRDS(object,file=filename,version=2)
  # write.table(object,file=filename)
}

#' Print scenario statistics
#' 
#' Print some statistics from the current scenario
#'
#' @param msg : message header
#' @param scenario : object resulting from current Rpath run
#'
#' @return None
#' 
printStatsScenario <- function(msg,scenario) {
  print("- - -")
  print(paste0(msg,", forcing$ForcedBio:     ",sum(scenario$forcing$ForcedBio)))
  print(paste0(msg,", forcing$ForcedMigrate: ",sum(scenario$forcing$ForcedMigrate)))
  print(paste0(msg,", start_state$Biomass:   ",sum(scenario$start_state$Biomass)))
}

#' Print simulation statistics
#' 
#' Print some statistics from the current simulation
#'
#' @param msg : message header
#' @param sim : object resulting from last Rpath simulation run
#'
#' @return None
#' 
printStatsSimulation <- function(msg,sim) {
  print("- - -")
  print(paste0(msg,", sim$out_Biomass:    ",sum(sim$out_Biomass)))
  print(paste0(msg,", sim$out_Catch:      ",sum(sim$out_Catch)))
  print(paste0(msg,", sim$out_Gear_Catch: ",sum(sim$out_Gear_Catch)))
}

#' Modify a scenario$fishing matrix
#' 
#' Modifies the appropriate columns from a scenario$fishing matrix.
#'
#' @param species : A vector of fish species
#' @param fleets : A vector of fleets
#' @param typeData : The type of forced data (i.e., Forced Effort, Forced FRate, Forced Catch)
#' @param forcingData : The forced data matrix
#'
#' @return Returns a matrix that's been updated withe the forced data
#' 
modifyFishingMatrix <- function(modNum,species,fleets,typeData,forcingData,model,randomNumberType) {
  group <- model$Group
  pb <- model$PB
  ForcedMatrix <- forcingData
  speciesOrFleets <- c()
  const1 <- 1
  const2 <- 1
  upperLimit <- 0.5
  usePBValue <- FALSE
  scaleFactorPB <- 0.0001 # Needed this to add enough randomness to the plots, else they'd be fairly smooth
  if (typeData == "Forced Effort") {
    speciesOrFleets <- fleets
  } else if (typeData == "Forced FRate" || typeData == "Forced Catch") {
    speciesOrFleets <- species
    const1 <- 0
    usePBValue <- TRUE
  } else {
    print(paste0("Error: Found invalid typeData of: ",typeData))
    return(ForcedMatrix)
  }
  index <- match(speciesOrFleets,group)
  if (usePBValue) {
    const2 <- scaleFactorPB*pb[index]
  }
  
  j <- 0
  for (i in 1:length(speciesOrFleets)) {
    item <- speciesOrFleets[i]
    vectorData           <- ForcedMatrix[,item]
    newVectorWithJitter <- c()
    for (value in vectorData) {
      j <- j + 1
      # Not sure why I need the [1] index here, this should always be just a single value but sometimes it's a list
      jitteredValue <- (const1+const2*randomNumber(upperLimit,randomNumberType))[1]
      newVectorWithJitter = append(newVectorWithJitter,jitteredValue)
    }
    ForcedMatrix[,item]  <- newVectorWithJitter
  }
  return(ForcedMatrix) 
}

#' Modify a scenario$forcing matrix
#' 
#' Modifies the appropriate columns from a scenario$forcing matrix.
#'
#' @param species : A vector of fish species
#' @param modifyType : The type of modification (i.e., Jittered or Stepped)
#' @param typeData : The type of forced data (i.e., Forced Bio, Forced Migrate)
#' @param forcingData : The forced data matrix
#' @param scenario : The REcosystem_scenario object 
#'
#' @return Returns a matrix that's been updated withe the forced data
#' 
modifyForcingMatrix <- function (modNum,species,modifyType,typeData,forcingData,scenario,randomNumberType) {
  ForcedMatrix <- forcingData
  numMonths <- nrow(ForcedMatrix)
  scaleFactors <- c(100,50,200,300,10000)
  scaleFactors <- c(.6,.6,.6,.6,.6)
  if ((typeData == FORCED_BIOMASS) || (typeData == FORCED_MIGRATION)) {
    for (i in 1:length(species)) {
      aSpecies <- species[[i]]
      speciesBiomass <- scenario$start_state$Biomass[aSpecies]
      startValue <- speciesBiomass
      if (typeData == FORCED_MIGRATION) {
        # rval <- (randomNumber(FORCED_MIGRATION_BIOMASS_PCT)+FORCED_MIGRATION_BIOMASS_PCT)/2.0
        rval <- randomNumber(FORCED_MIGRATION_BIOMASS_PCT,randomNumberType)
        startValue <- rval
        # print(paste0("species: ",aSpecies,", biomass: ",speciesBiomass,", randomNum: ",rval,", startValue: ",startValue))
      }
      if (modifyType == JITTERED) {
        # print(paste0("start value: ",i,", ",startValue))
        ForcedMatrix[,aSpecies] <- createJitterVectorFromValue(typeData, startValue, numMonths, 
                                                               "Months","Biomass (mt/km\U00B2)",
                                                               paste0(typeData,' with ',modifyType,' Noise - ',aSpecies),
                                                               POSITIVE_AND_NEGATIVE)
      } else {
        stepType <- ((i-1)%%3)+1 # Only current step types are 1, 2, or 3
        ForcedMatrix[,aSpecies] <- stepifyBiomass(typeData, startValue, numMonths, stepType, "Months","Biomass (mt/km\U00B2)",
                                                  paste0(typeData,' with ',modifyType,' Noise - ',aSpecies), scaleFactors[i])
      }
    }
  }
  return(ForcedMatrix)
}
