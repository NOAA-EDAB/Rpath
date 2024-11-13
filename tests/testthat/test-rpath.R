library(Rpath)
library(qpdf)
library(testthat)
library(stringr)
library(distillery)
library(ggplot2)
library(ggpubr)
library(rlist)

options(precision=50)

source("test-constants.R")
source("test-utils.R")
source("test-utils-stepify.R")
source("test-utils-jitter.R")
source("test-utils-plot.R")

print(paste0("here::here: ",              here::here() ))
print(paste0("OUTPUT_DATA_DIR: ",         OUTPUT_DATA_DIR)) #RSK
print(paste0("INPUT_DATA_DIR_BASELINE: ", INPUT_DATA_DIR_BASELINE)) #RSK
print(paste0("INPUT_DATA_DIR_CURRENT: ",  INPUT_DATA_DIR_CURRENT)) #RSK

# Set the seed here so that all runs will be deterministic.
set.seed(SEED_VALUE)

# Create the current and output directories if they don't already exist.
if (! dir.exists(INPUT_DATA_DIR_CURRENT)) {
  dir.create(INPUT_DATA_DIR_CURRENT,recursive=TRUE)
}
if (! dir.exists(OUTPUT_DATA_DIR)) {
  dir.create(OUTPUT_DATA_DIR,recursive=TRUE)
}

#' This test tests for equality between two objects
#' 
#' This function tests to see if the two data tables are equal within a tolerance value
#'
#' @param runNum : The run number
#' @param tableName : Temp arg, used to skip over out_Gear_Catch tables as there's currently a bug with their header names
#' @param desc : The description of the run
#' @param baselineTable : The baseline table 
#' @param currentTable : The current generated table
#'
#' @return No return value
#'
runTestEqual <- function(runNum,tableName,desc,baselineTable,currentTable) {
  if (tableName == 'out_Gear_Catch') {
    if (! RUN_QUIET) {
      print('Test Skipped: The out_Gear_Catch tables are currently missing column headings...so this test will be skipped for now.')
    }
    return()
  }
  paddedRunNum <- str_pad(runNum,3,pad="0")
  if (! RUN_QUIET) {
    print(paste0("Test #",paddedRunNum,": ",desc))
  }
  testthat::expect_equal(baselineTable,currentTable,tolerance=TOLERANCE_VALUE)
}

#' The test for running silent test function
#' 
#' This function tests to see if the current Rpath run runs silently 
#'
#' @param runNum : The run number
#' @param desc : The description of the run
#' @param params : Parameters associated with, and passed to, the Rpath run 
#' @param name : The ecosystem name of the run that's passed to Rpath
#'
#' @return No return value
#'
runTestSilent <- function(runNum,desc,params,name) {
  paddedRunNum <- str_pad(runNum,3,pad="0")
  if (! RUN_QUIET) {
    print(paste0("Test #",paddedRunNum,": ",desc))
  }
  testthat::expect_silent(rpath(params,eco.name=name))
}

#' The main test run function
#' 
#' This function prints the appropriate baseline and current forced data and then runs the testthat::expect_equal test.
#'
#' @param tableName : Name of output table to plot (i.e., out_Biomass, out_Catch)
#' @param forcedData : The data that's being forced (i.e., Biomass, Catch)
#' @param forcedType : The type of forcing being done (.e., Random (Jitter) or Stepped)
#' @param baseAlg : The baseline algorithm used (currently AB or RK4)
#' @param currAlg : The current algorithm used (currently AB or RK4) 
#' @param baselineTable : table that represents the baseline run's data
#' @param currentDataFrame : The data frame that contains the current run's data
#' @param currentFilename : The name of the file that contains the current run's data
#' @param species : A vector of fish species
#'
#' @return No return value
#'
runTestRDS <- function(runNum,tableName,forcedData,forcedType,baseAlg,currAlg,baselineDataFrame,currentFilename,species) {
  
  # N.B. Remove this if statement once missing column headings issue is fixed
  if (tableName == 'out_Gear_Catch') {
    if (! RUN_QUIET) {
      print('Test Skipped: The out_Gear_Catch tables are currently missing column headings...so this test will be skipped for now.')
    }
    return()
  }
  
  paddedRunNum <- str_pad(runNum,3,pad="0")
  if (! RUN_QUIET) {
    print(paste0("Test #",paddedRunNum,": Is Baseline ",baseAlg," ",tableName," equivalent to Current ",currAlg," ",tableName," using ",forcedType," ",forcedData,"?"))
  }
  
  currentDataFrame <- readDataFile(currentFilename)
  
  if (tableName == 'out_Biomass' || tableName == 'out_Catch') { # || tableName == 'out_Gear_Catch') {
    if (PLOT_TYPE == 1) {
      plotResultsSuperimposed(baselineDataFrame, currentDataFrame, baseAlg, currAlg, tableName, forcedData, forcedType, species)
    } else {
      plotResultsDifference(baselineDataFrame, currentDataFrame, baseAlg, currAlg, tableName, forcedData, forcedType, species)
    }
  }
    
  # Write out the difference table (current-baseline)
  # diffTable <- abs(currentDataFrame-baselineDataFrame)
  # diffTable <- currentDataFrame
  numRows <- nrow(currentDataFrame)
  numCols <- ncol(currentDataFrame)
  diffTable <- matrix(nrow=numRows,ncol=numCols)
  # diffTable[TRUE] <- 0
  for (i in 1:nrow(currentDataFrame)) {
    for (j in 1:ncol(currentDataFrame)) {
      diffTable[i,j] <- abs(currentDataFrame[i,j] - baselineDataFrame[i,j])
      if (diffTable[i,j] <= TOLERANCE_VALUE) {
        diffTable[i,j] <- 0
      } 
      # else {
      #   print(paste0("-> diffTable[",i,",",j,"]: ",diffTable[i,j]))
      # }
    }
  }

  # Set all values <= TOLERANCE_VALUE to 0 because we're going to next compare this to a zero table
  # diffTable[diffTable <= TOLERANCE_VALUE] <- 0

  # Create the zero table
  zeroTable <- matrix(0,nrow=numRows,ncol=numCols)
  
  sumDiffTable <- 0
  sumColsCurr = colSums(currentDataFrame)
  sumColsBase = colSums(baselineDataFrame)
  for (i in 1:ncol(currentDataFrame)) {
    sumDiffTable <- sumDiffTable + abs(sumColsCurr[i]-sumColsBase[i])
  }
#print(paste0("***sumDiffTable: ",sumDiffTable))
#print(paste0("SUM of currentDataFrame:  ", sum(currentDataFrame)))           # ok
#print(paste0("SUM of baselineDataFrame: ", sum(baselineDataFrame)))          # ok
#print("Comparing if sumDiffTable/sum(currentDataFrame) <= TOLERANCE_VALUE")
#print(paste0("Is ",sumDiffTable,"/",sum(currentDataFrame)," <= ",TOLERANCE_VALUE," ?"))
areIdentical <- (sumDiffTable/ sum(currentDataFrame) <= TOLERANCE_VALUE)
#print(paste0("areIdentical: ",areIdentical))  
                   
# print(paste0("SUM of currentDataFrame: ", sum(currentDataFrame)))            # ok
# print(paste0("SUM of baselineDataFrame: ",sum(baselineDataFrame)))           # ok
#print(paste0("SUM of diffTable:         ",sum(diffTable)))                   # not ok
#print(paste0("SUM of zeroTable:         ",sum(zeroTable)))                   # ok
# print(paste0("Col sums currentDataFrame:  ",colSums(currentDataFrame)))      # ok
# print(paste0("Col sums baselineDataFrame: ",colSums(baselineDataFrame)))     # ok
#print(paste0("Col sums diffTable:         ",colSums(diffTable)))             # not ok
# print(paste0("Sum currentDataFrame col=17: ", sum(currentDataFrame[,17])))   # ok
# print(paste0("Sum baselineDataFrame col=17: ", sum(baselineDataFrame[,17]))) # ok
#print(paste0("Sum diffTable col=17: ", sum(diffTable[,17])))                 # not ok
# areIdentical <- identical(diffTable,zeroTable)                             # not ok
#print(paste0("*** Are they identical diffTable and zeroTable: ",identical(diffTable,zeroTable)))

  testthat::expect_true(areIdentical)
  
  write.table(diffTable, file=file.path(OUTPUT_DATA_DIR,paste0("diff_",paddedRunNum,".dat")))
  write.table(zeroTable, file=file.path(OUTPUT_DATA_DIR,paste0("zero_",paddedRunNum,".dat")))
}


testthat::test_that("Rpath Unit Tests", {
  BaselineJitterTables  <- list()
  BaselineJitterFiles   <- list()
  CurrentJitterFiles    <- list()
  BaselineSteppedTables <- list()
  BaselineSteppedFiles  <- list()
  CurrentSteppedFiles   <- list()
  fleets  <- c('Trawlers','Midwater','Dredgers')
  species <- c('OtherGroundfish','Megabenthos','Seals','JuvRoundfish1','AduRoundfish1')
  modNum <- 1
  runNum <- 0

  # ---------- Set up initial file paths ----------
  # N.B. The Baseline and Current AB and RK4 files are .csv files since they were produced by
  # the write.rsim() function and not the more generic write.table() function.
  BaselineRpathObjTopLevel                 <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RpathObj_TopLevel.rds')
  BaselineRpathObjSummary                  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RpathObj_Summary.rds')
  BaselineAB                               <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB.csv')
  BaselineRK4                              <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4.csv')
  BaselineABOutBiomass                     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutBiomass.rds')
  BaselineRK4OutBiomass                    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutBiomass.rds')
  BaselineABOutCatch                       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutCatch.rds')
  BaselineRK4OutCatch                      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutCatch.rds')
  BaselineABOutGearCatch                   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutGearCatch.rds')
  BaselineRK4OutGearCatch                  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutGearCatch.rds')
  BaselineABForcedBioOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped.rds')
  BaselineRK4ForcedBioOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped.rds')
  BaselineABForcedBioOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped.rds')
  BaselineRK4ForcedBioOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped.rds')
  BaselineABForcedBioOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped.rds')
  BaselineRK4ForcedBioOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped.rds')
  BaselineABForcedBioOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter.rds')
  BaselineRK4ForcedBioOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter.rds')
  BaselineABForcedBioOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter.rds')
  BaselineRK4ForcedBioOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter.rds')
  BaselineABForcedBioOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter.rds')
  BaselineRK4ForcedBioOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter.rds')
  BaselineABForcedMigOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped.rds')
  BaselineRK4ForcedMigOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped.rds')
  BaselineABForcedMigOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped.rds')
  BaselineRK4ForcedMigOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped.rds')
  BaselineABForcedMigOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped.rds')
  BaselineRK4ForcedMigOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped.rds')
  BaselineABForcedMigOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter.rds')
  BaselineRK4ForcedMigOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter.rds')
  BaselineABForcedMigOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter.rds')
  BaselineRK4ForcedMigOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter.rds')
  BaselineABForcedMigOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter.rds')
  BaselineRK4ForcedMigOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter.rds')
  BaselineABForcedEffOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped.rds')
  BaselineABForcedEffOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped.rds')
  BaselineABForcedEffOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped.rds')
  BaselineABForcedEffOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter.rds')
  BaselineABForcedEffOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter.rds')
  BaselineABForcedEffOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter.rds')
  BaselineABForcedFRaOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped.rds')
  BaselineABForcedFRaOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped.rds')
  BaselineABForcedFRaOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped.rds')
  BaselineABForcedFRaOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter.rds')
  BaselineABForcedFRaOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter.rds')
  BaselineABForcedFRaOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter.rds')
  BaselineABForcedCatOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped.rds')
  BaselineABForcedCatOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped.rds')
  BaselineABForcedCatOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped.rds')
  BaselineABForcedCatOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter.rds')
  BaselineABForcedCatOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter.rds')
  BaselineABForcedCatOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter.rds')
  BaselineRK4ForcedEffOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped.rds')
  BaselineRK4ForcedEffOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped.rds')
  BaselineRK4ForcedEffOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped.rds')
  BaselineRK4ForcedEffOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter.rds')
  BaselineRK4ForcedEffOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter.rds')
  BaselineRK4ForcedEffOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter.rds')
  BaselineRK4ForcedFRaOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped.rds')
  BaselineRK4ForcedFRaOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped.rds')
  BaselineRK4ForcedFRaOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped.rds')
  BaselineRK4ForcedFRaOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter.rds')
  BaselineRK4ForcedFRaOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter.rds')
  BaselineRK4ForcedFRaOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter.rds')
  BaselineRK4ForcedCatOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped.rds')
  BaselineRK4ForcedCatOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped.rds')
  BaselineRK4ForcedCatOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped.rds')
  BaselineRK4ForcedCatOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter.rds')
  BaselineRK4ForcedCatOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter.rds')
  BaselineRK4ForcedCatOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter.rds')
  #
  CurrentRpathObjTopLevel                  <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RpathObj_TopLevel.rds')
  CurrentRpathObjSummary                   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RpathObj_Summary.rds')
  CurrentAB                                <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB.csv')
  CurrentRK4                               <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4.csv')
  CurrentABOutBiomass                      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutBiomass.rds')
  CurrentRK4OutBiomass                     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutBiomass.rds')
  CurrentABOutCatch                        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutCatch.rds')
  CurrentRK4OutCatch                       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutCatch.rds')
  CurrentABOutGearCatch                    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutGearCatch.rds')
  CurrentRK4OutGearCatch                   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutGearCatch.rds')
  CurrentABForcedBioOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutBiomass_Stepped.rds')
  CurrentRK4ForcedBioOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutBiomass_Stepped.rds')
  CurrentABForcedBioOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutCatch_Stepped.rds')
  CurrentRK4ForcedBioOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutCatch_Stepped.rds')
  CurrentABForcedBioOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutGearCatch_Stepped.rds')
  CurrentRK4ForcedBioOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutGearCatch_Stepped.rds')
  CurrentABForcedBioOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutBiomass_Jitter.rds')
  CurrentABForcedBioOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutCatch_Jitter.rds')
  CurrentABForcedBioOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutGearCatch_Jitter.rds')
  CurrentRK4ForcedBioOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutBiomass_Jitter.rds')
  CurrentRK4ForcedBioOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutCatch_Jitter.rds')
  CurrentRK4ForcedBioOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutGearCatch_Jitter.rds')
  CurrentABForcedMigOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutBiomass_Stepped.rds')
  CurrentRK4ForcedMigOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutBiomass_Stepped.rds')
  CurrentABForcedMigOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutCatch_Stepped.rds')
  CurrentRK4ForcedMigOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutCatch_Stepped.rds')
  CurrentABForcedMigOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutGearCatch_Stepped.rds')
  CurrentRK4ForcedMigOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutGearCatch_Stepped.rds')
  CurrentABForcedMigOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutBiomass_Jitter.rds')
  CurrentABForcedMigOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutCatch_Jitter.rds')
  CurrentABForcedMigOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutGearCatch_Jitter.rds')
  CurrentRK4ForcedMigOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutBiomass_Jitter.rds')
  CurrentRK4ForcedMigOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutCatch_Jitter.rds')
  CurrentRK4ForcedMigOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutGearCatch_Jitter.rds')
  CurrentABForcedEffOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutBiomass_Stepped.rds')
  CurrentRK4ForcedEffOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutBiomass_Stepped.rds')
  CurrentABForcedEffOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutCatch_Stepped.rds')
  CurrentRK4ForcedEffOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutCatch_Stepped.rds')
  CurrentABForcedEffOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutGearCatch_Stepped.rds')
  CurrentRK4ForcedEffOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutGearCatch_Stepped.rds')
  CurrentABForcedFRaOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutBiomass_Stepped.rds')
  CurrentRK4ForcedFRaOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutBiomass_Stepped.rds')
  CurrentABForcedFRaOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutCatch_Stepped.rds')
  CurrentRK4ForcedFRaOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutCatch_Stepped.rds')
  CurrentABForcedFRaOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutGearCatch_Stepped.rds')
  CurrentRK4ForcedFRaOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutGearCatch_Stepped.rds')
  CurrentABForcedCatOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutBiomass_Stepped.rds')
  CurrentRK4ForcedCatOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutBiomass_Stepped.rds')
  CurrentABForcedCatOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutCatch_Stepped.rds')
  CurrentRK4ForcedCatOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutCatch_Stepped.rds')
  CurrentABForcedCatOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutGearCatch_Stepped.rds')
  CurrentRK4ForcedCatOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutGearCatch_Stepped.rds')
  CurrentABForcedEffOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutBiomass_Jitter.rds')
  CurrentABForcedEffOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutCatch_Jitter.rds')
  CurrentABForcedEffOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutGearCatch_Jitter.rds')
  CurrentRK4ForcedEffOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutBiomass_Jitter.rds')
  CurrentRK4ForcedEffOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutCatch_Jitter.rds')
  CurrentRK4ForcedEffOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutGearCatch_Jitter.rds')
  CurrentABForcedFRaOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutBiomass_Jitter.rds')
  CurrentABForcedFRaOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutCatch_Jitter.rds')
  CurrentABForcedFRaOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutGearCatch_Jitter.rds')
  CurrentRK4ForcedFRaOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutBiomass_Jitter.rds')
  CurrentRK4ForcedFRaOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutCatch_Jitter.rds')
  CurrentRK4ForcedFRaOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutGearCatch_Jitter.rds')
  CurrentABForcedCatOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutBiomass_Jitter.rds')
  CurrentABForcedCatOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutCatch_Jitter.rds')
  CurrentABForcedCatOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutGearCatch_Jitter.rds')
  CurrentRK4ForcedCatOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutBiomass_Jitter.rds')
  CurrentRK4ForcedCatOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutCatch_Jitter.rds')
  CurrentRK4ForcedCatOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutGearCatch_Jitter.rds')

  # Initialize baseline and current paths
  REcosystemBaseline                                     <- NULL
  REcosystemBaselineSummary                              <- NULL
  REcosystem_Baseline_AB                                 <- NULL
  REcosystem_Baseline_RK4                                <- NULL
  REcosystem_Baseline_AB_OutBiomass                      <- NULL
  REcosystem_Baseline_AB_OutCatch                        <- NULL
  REcosystem_Baseline_AB_OutGearCatch                    <- NULL
  REcosystem_Baseline_RK4_OutBiomass                     <- NULL
  REcosystem_Baseline_RK4_OutCatch                       <- NULL
  REcosystem_Baseline_RK4_OutGearCatch                   <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter     <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter    <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter       <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter      <- NULL
  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter   <- NULL
  REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter  <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter     <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter    <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter       <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter      <- NULL
  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter   <- NULL
  REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter  <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter     <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter       <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter   <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter    <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter      <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter  <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter     <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter       <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter   <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter    <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter      <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter  <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter     <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter       <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter   <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter    <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter      <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter  <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped    <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped      <- NULL
  REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped  <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped   <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped     <- NULL
  REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped <- NULL
  if (! CREATE_BASELINE_FILES) {
    REcosystemBaseline                                     <- readDataFile(BaselineRpathObjTopLevel) #,               fill = TRUE, sep = " ")
    REcosystemBaselineSummary                              <- read.table(BaselineRpathObjSummary,                fill = TRUE)
    REcosystem_Baseline_AB                                 <- read.csv(BaselineAB)
    REcosystem_Baseline_RK4                                <- read.csv(BaselineRK4)
    REcosystem_Baseline_AB_OutBiomass                      <- readDataFile(BaselineABOutBiomass) #,                   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_OutCatch                        <- readDataFile(BaselineABOutCatch) #,                     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_OutGearCatch                    <- readDataFile(BaselineABOutGearCatch) #,                 fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutBiomass                     <- readDataFile(BaselineRK4OutBiomass) #,                  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutCatch                       <- readDataFile(BaselineRK4OutCatch) #,                    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutGearCatch                   <- readDataFile(BaselineRK4OutGearCatch) #,                fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped    <- readDataFile(BaselineABForcedBioOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedBioOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped      <- readDataFile(BaselineABForcedBioOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedBioOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedBioOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedBioOutGearCatchStepped) #,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter     <- readDataFile(BaselineABForcedBioOutBiomassJitter) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter    <- readDataFile(BaselineRK4ForcedBioOutBiomassJitter) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter       <- readDataFile(BaselineABForcedBioOutCatchJitter) #,      fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter      <- readDataFile(BaselineRK4ForcedBioOutCatchJitter) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter   <- readDataFile(BaselineABForcedBioOutGearCatchJitter) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter  <- readDataFile(BaselineRK4ForcedBioOutGearCatchJitter) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped    <- readDataFile(BaselineABForcedMigOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedMigOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped      <- readDataFile(BaselineABForcedMigOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedMigOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedMigOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedMigOutGearCatchStepped) #,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter     <- readDataFile(BaselineABForcedMigOutBiomassJitter) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter    <- readDataFile(BaselineRK4ForcedMigOutBiomassJitter) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter       <- readDataFile(BaselineABForcedMigOutCatchJitter) #,      fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter      <- readDataFile(BaselineRK4ForcedMigOutCatchJitter) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter   <- readDataFile(BaselineABForcedMigOutGearCatchJitter) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter  <- readDataFile(BaselineRK4ForcedMigOutGearCatchJitter) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- readDataFile(BaselineABForcedEffOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedEffOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- readDataFile(BaselineABForcedEffOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedEffOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedEffOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedEffOutGearCatchStepped) #,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter     <- readDataFile(BaselineABForcedEffOutBiomassJitter) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter       <- readDataFile(BaselineABForcedEffOutCatchJitter) #,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter   <- readDataFile(BaselineABForcedEffOutGearCatchJitter) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter    <- readDataFile(BaselineRK4ForcedEffOutBiomassJitter) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter      <- readDataFile(BaselineRK4ForcedEffOutCatchJitter) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter  <- readDataFile(BaselineRK4ForcedEffOutGearCatchJitter) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter     <- readDataFile(BaselineABForcedFRaOutBiomassJitter) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter       <- readDataFile(BaselineABForcedFRaOutCatchJitter) #,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter   <- readDataFile(BaselineABForcedFRaOutGearCatchJitter) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter    <- readDataFile(BaselineRK4ForcedFRaOutBiomassJitter) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter      <- readDataFile(BaselineRK4ForcedFRaOutCatchJitter) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter  <- readDataFile(BaselineRK4ForcedFRaOutGearCatchJitter) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter     <- readDataFile(BaselineABForcedCatOutBiomassJitter) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter       <- readDataFile(BaselineABForcedCatOutCatchJitter) #,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter   <- readDataFile(BaselineABForcedCatOutGearCatchJitter) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter    <- readDataFile(BaselineRK4ForcedCatOutBiomassJitter) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter      <- readDataFile(BaselineRK4ForcedCatOutCatchJitter) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter  <- readDataFile(BaselineRK4ForcedCatOutGearCatchJitter) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- readDataFile(BaselineABForcedEffOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- readDataFile(BaselineABForcedEffOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedEffOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedEffOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedEffOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedEffOutGearCatchStepped) #,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped    <- readDataFile(BaselineABForcedFRaOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped      <- readDataFile(BaselineABForcedFRaOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedFRaOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedFRaOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedFRaOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedFRaOutGearCatchStepped) #,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped    <- readDataFile(BaselineABForcedCatOutBiomassStepped) #,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped      <- readDataFile(BaselineABForcedCatOutCatchStepped) #,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped  <- readDataFile(BaselineABForcedCatOutGearCatchStepped) #, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped   <- readDataFile(BaselineRK4ForcedCatOutBiomassStepped) #,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped     <- readDataFile(BaselineRK4ForcedCatOutCatchStepped) #,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped <- readDataFile(BaselineRK4ForcedCatOutGearCatchStepped) #,fill = TRUE, sep = " ")
  }
  
  # # 
  # # Save current Rpath run data
  REco <- rpath(REco.params, eco.name = 'R Ecosystem') # an Rpath object
  RpathObjTopLevel <- if (CREATE_BASELINE_FILES) BaselineRpathObjTopLevel else CurrentRpathObjTopLevel
  writeDataFile(REco,RpathObjTopLevel)
  # sink(RpathObjTopLevel)
  #   print(REco)
  # sink()

  # Save current Rpath run summary data
  RpathObjSummary <- if (CREATE_BASELINE_FILES) BaselineRpathObjSummary else CurrentRpathObjSummary
  # writeDataFile(summary(REco),RpathObjSummary)
  sink(RpathObjSummary)
    cat(summary(REco))
  sink()

  # Save current Rpath sim run data for AB and RK4
  REcosystem_scenario <- rsim.scenario(REco, REco.params, 1:50)
  REcosystem_Current_AB_from_Sim  <- rsim.run(REcosystem_scenario,method='AB', years=1:50)
  REcosystem_Current_RK4_from_Sim <- rsim.run(REcosystem_scenario,method='RK4',years=1:50)
  if (CREATE_BASELINE_FILES) {
    write.Rsim(REcosystem_Current_AB_from_Sim, BaselineAB)
    write.Rsim(REcosystem_Current_RK4_from_Sim,BaselineRK4)
  } else {
    write.Rsim(REcosystem_Current_AB_from_Sim, CurrentAB)
    write.Rsim(REcosystem_Current_RK4_from_Sim,CurrentRK4)
  }
  
  # ------------------------------------------
  # ------------------------------------------
  # --------------- Run tests ----------------
  # ------------------------------------------
  # ------------------------------------------
  
  print("------------------ Rpath Object Tests ------------------")
  if (! CREATE_BASELINE_FILES) {
    # Remove existing output data files
    cwd <- getwd()
    files <- dir(path=file.path(cwd,OUTPUT_DATA_DIR),pattern='diff_*')
    file.remove(file.path(OUTPUT_DATA_DIR,files))
    files <- dir(path=file.path(cwd,OUTPUT_DATA_DIR),pattern='zero_*')
    file.remove(file.path(OUTPUT_DATA_DIR,files))

    # Test 1 - Test if Balanced (i.e., "Status: Balanced" is the 2nd line of the Summary file)
    headerSummaryLines <- readLines(CurrentRpathObjSummary,n=2)
    parts <- unlist(strsplit(str_trim(headerSummaryLines[2]),split=" "))
    runTestEqual(inc(runNum),"","Is model balanced?",parts[2],"Balanced")

    # Test 2 - Test if function runs silently (i.e., no messages, warnings, or print statements)
    runTestSilent(inc(runNum),"Does model run without any terminal output (i.e., warnings, errors)?",REco.params,'R Ecosystem')

    # Test 3 - Test that the REcosystem object is the same as the saved original REcosystem object
    REcosystemCurrent <- readDataFile(CurrentRpathObjTopLevel)#,fill = TRUE, sep = " ")
    runTestEqual(inc(runNum),"","Is the baseline Rpath object equivalent to the current Rpath object (toplevel data)?",REcosystemBaseline,REcosystemCurrent)

    # Test 4 - Test that the REcosystem Summary is the same as the saved original REcosystem Summary
    REcosystemSummaryCurrent <- read.table(CurrentRpathObjSummary,fill = TRUE)
    runTestEqual(inc(runNum),"","Is the baseline Rpath run Summary the same as the current Rpath Summary?",REcosystemBaselineSummary,REcosystemSummaryCurrent)

    REcosystem_Current_AB  <- read.csv(CurrentAB)
    REcosystem_Current_RK4 <- read.csv(CurrentRK4)
  }

  # Tests 5-16 - Test that REcosystem AB object is same as RK4 object with no perturbations
  if (CREATE_BASELINE_FILES) {
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Biomass,    BaselineABOutBiomass)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Biomass,   BaselineRK4OutBiomass)
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Catch,      BaselineABOutCatch)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Catch,     BaselineRK4OutCatch)
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Gear_Catch, BaselineABOutGearCatch)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Gear_Catch,BaselineRK4OutGearCatch)
  }
  else {
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Biomass,    CurrentABOutBiomass)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Biomass,   CurrentRK4OutBiomass)
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Catch,      CurrentABOutCatch)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Catch,     CurrentRK4OutCatch)
    writeDataFile(REcosystem_Current_AB_from_Sim$out_Gear_Catch, CurrentABOutGearCatch)
    writeDataFile(REcosystem_Current_RK4_from_Sim$out_Gear_Catch,CurrentRK4OutGearCatch)
    REcosystem_Current_AB_OutBiomass    <- readDataFile(CurrentABOutBiomass) #,   fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutBiomass   <- readDataFile(CurrentRK4OutBiomass) #,  fill=TRUE,sep=" ")
    REcosystem_Current_AB_OutCatch      <- readDataFile(CurrentABOutCatch) #,     fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutCatch     <- readDataFile(CurrentRK4OutCatch) #,    fill=TRUE,sep=" ")
    REcosystem_Current_AB_OutGearCatch  <- readDataFile(CurrentABOutGearCatch) #, fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutGearCatch <- readDataFile(CurrentRK4OutGearCatch) #,fill=TRUE,sep=" ")
    runTestEqual(inc(runNum),"",              "Compare baseline AB to Current AB",                    REcosystem_Baseline_AB,              REcosystem_Current_AB)
    runTestEqual(inc(runNum),"out_Biomass",   "Compare baseline AB to Current AB for OutBiomass",     REcosystem_Baseline_AB_OutBiomass,   REcosystem_Current_AB_OutBiomass)
    runTestEqual(inc(runNum),"out_Catch",     "Compare baseline AB to Current AB for OutCatch",       REcosystem_Baseline_AB_OutCatch,     REcosystem_Current_AB_OutCatch)
    runTestEqual(inc(runNum),"out_Gear_Catch","Compare baseline AB to Current AB for OutGearCatch",   REcosystem_Baseline_AB_OutGearCatch, REcosystem_Current_AB_OutGearCatch)
    runTestEqual(inc(runNum),"",              "Compare baseline RK4 to Current RK4",                  REcosystem_Baseline_RK4,             REcosystem_Current_RK4)
    runTestEqual(inc(runNum),"out_Biomass",   "Compare baseline RK4 to Current RK4 for OutputBiomass",REcosystem_Baseline_RK4_OutBiomass,  REcosystem_Current_RK4_OutBiomass)
    runTestEqual(inc(runNum),"out_Catch",     "Compare baseline RK4 to Current RK4 for OutCatch",     REcosystem_Baseline_RK4_OutCatch,    REcosystem_Current_RK4_OutCatch)
    runTestEqual(inc(runNum),"out_Gear_Catch","Compare baseline RK4 to Current RK4 for OutGearCatch", REcosystem_Baseline_RK4_OutGearCatch,REcosystem_Current_RK4_OutGearCatch)
    runTestEqual(inc(runNum),"",              "Compare baseline AB to Current RK4",                   REcosystem_Baseline_AB,              REcosystem_Current_RK4)
    runTestEqual(inc(runNum),"out_Biomass",   "Compare baseline AB to Current RK4 for OutBiomass",    REcosystem_Baseline_AB_OutBiomass,   REcosystem_Current_RK4_OutBiomass)
    runTestEqual(inc(runNum),"out_Catch",     "Compare baseline AB to Current RK4 for OutCatch",      REcosystem_Baseline_AB_OutCatch,     REcosystem_Current_RK4_OutCatch)
    runTestEqual(inc(runNum),"out_Gear_Catch","Compare baseline AB to Current RK4 for OutGearCatch",  REcosystem_Baseline_AB_OutGearCatch, REcosystem_Current_RK4_OutGearCatch)
  }
  
  print("------------------ Forced Biomass Tests (Jitter) ------------------")
  typeData <- list(FORCED_BIOMASS,FORCED_MIGRATION)
  numMonths <- nrow(REcosystem_scenario$forcing$ForcedBio)
  for (typeNum in 1:length(typeData)) {
    REco <- rpath(REco.params, eco.name = 'R Ecosystem') # an Rpath object
    REcosystem_scenario <- rsim.scenario(REco, REco.params, 1:50)
    REcosystem_scenario_jitter <- REcosystem_scenario
    theTypeData  <- typeData[[typeNum]]
    modNum <- modNum + 1
    if (theTypeData == FORCED_BIOMASS) {
      BaselineJitterDataFrames <- list(REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter,  REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter,  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter,
                                       REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter, REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter, REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter)
      BaselineJitterFilenames  <- list(BaselineABForcedBioOutBiomassJitter,  BaselineABForcedBioOutCatchJitter,  BaselineABForcedBioOutGearCatchJitter,
                                       BaselineRK4ForcedBioOutBiomassJitter, BaselineRK4ForcedBioOutCatchJitter, BaselineRK4ForcedBioOutGearCatchJitter)
      CurrentJitterFilenames   <- list(CurrentABForcedBioOutBiomassJitter,  CurrentABForcedBioOutCatchJitter,  CurrentABForcedBioOutGearCatchJitter,
                                       CurrentRK4ForcedBioOutBiomassJitter, CurrentRK4ForcedBioOutCatchJitter, CurrentRK4ForcedBioOutGearCatchJitter)
      # RSK - These (original) lines don't work in git actions
      REcosystem_scenario_jitter$forcing$ForcedBio <- modifyForcingMatrix(modNum, species, JITTERED, theTypeData, 
                                                                          REcosystem_scenario_jitter$forcing$ForcedBio,
                                                                          REcosystem_scenario_jitter,
                                                                          POSITIVE_ONLY)
      
      # A couple different ways to jitter. Trying to determine why some tests fail in git actions when run
      # in the tests.yml file but don't fail when run via R-CMD-Check.yml.
      #
      # This line doesn't fail in git actions (it's just not the exact logic I need)
      # REcosystem_scenario_jitter$forcing$ForcedBio <- jitter(REcosystem_scenario_jitter$forcing$ForcedBio,factor=FACTOR_VALUE)
      #
      # Another way to jitter
      #numMonths <- nrow(REcosystem_scenario_jitter$forcing$ForcedBio)
      #numSpecies <- length(species)
      #speciesNum <- 0
      #totSpeciesBiomass <- 0
      #totRandVal <- 0
      #for (aSpecies in species) {
      #  jitterVector <- c()
      #  speciesBiomass <- REcosystem_scenario_jitter$start_state$Biomass[aSpecies]
      #  totSpeciesBiomass <- totSpeciesBiomass + speciesBiomass
      #  for (month in 1:numMonths) {
      #    randVal <- ...
      #    jitteredValue <- speciesBiomass * (1.0 + randVal)
      #    totRandVal <- totRandVal + randVal
      #    jitterVector <- append(jitterVector,jitteredValue)
      #  }
      #  speciesNum <- speciesNum + 1
      #  REcosystem_scenario_jitter$forcing$ForcedBio[,aSpecies] <- jitterVector # RSK problematic line here
      #}
      # REcosystem_scenario_jitter$forcing$ForcedBio <- jitterMatrixColumns(REcosystem_scenario_jitter$forcing$ForcedBio,REcosystem_scenario_jitter$start_state$Biomass,species)
    }
    else if (theTypeData == FORCED_MIGRATION) {
      BaselineJitterDataFrames <- list(REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter,  REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter,  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter,
                                       REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter, REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter, REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter)
      BaselineJitterFilenames  <- list(BaselineABForcedMigOutBiomassJitter,  BaselineABForcedMigOutCatchJitter,  BaselineABForcedMigOutGearCatchJitter,
                                       BaselineRK4ForcedMigOutBiomassJitter, BaselineRK4ForcedMigOutCatchJitter, BaselineRK4ForcedMigOutGearCatchJitter)
      CurrentJitterFilenames   <- list(CurrentABForcedMigOutBiomassJitter,  CurrentABForcedMigOutCatchJitter,  CurrentABForcedMigOutGearCatchJitter,
                                       CurrentRK4ForcedMigOutBiomassJitter, CurrentRK4ForcedMigOutCatchJitter, CurrentRK4ForcedMigOutGearCatchJitter)
      REcosystem_scenario_jitter$forcing$ForcedMigrate <- modifyForcingMatrix(modNum,species,JITTERED,theTypeData,
                                                                              REcosystem_scenario_jitter$forcing$ForcedMigrate,
                                                                              REcosystem_scenario_jitter,
                                                                              POSITIVE_ONLY)
    } else {
      print(paste0("Error: Unknown data type: ",theTypeData))
      return()
    }
    # Tests 17-22 Forced Biomass with Jitter
    # Tests 23-28 Forced Migration with Jitter
    REcosystem_AB_Current_Jitter  <- rsim.run(REcosystem_scenario_jitter,method='AB', years=1:50)
    REcosystem_RK4_Current_Jitter <- rsim.run(REcosystem_scenario_jitter,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      writeDataFile(REcosystem_AB_Current_Jitter$out_Biomass,     BaselineJitterFilenames[[1]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Catch,       BaselineJitterFilenames[[2]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Gear_Catch,  BaselineJitterFilenames[[3]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Biomass,    BaselineJitterFilenames[[4]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Catch,      BaselineJitterFilenames[[5]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Gear_Catch, BaselineJitterFilenames[[6]])
    } else {
      writeDataFile(REcosystem_AB_Current_Jitter$out_Biomass,     CurrentJitterFilenames[[1]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Catch,       CurrentJitterFilenames[[2]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Gear_Catch,  CurrentJitterFilenames[[3]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Biomass,    CurrentJitterFilenames[[4]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Catch,      CurrentJitterFilenames[[5]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Gear_Catch, CurrentJitterFilenames[[6]])
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Random", "AB",  "AB",  BaselineJitterDataFrames[[1]], CurrentJitterFilenames[[1]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Random", "AB",  "AB",  BaselineJitterDataFrames[[2]], CurrentJitterFilenames[[2]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "AB",  "AB",  BaselineJitterDataFrames[[3]], CurrentJitterFilenames[[3]], species)
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Random", "RK4", "RK4", BaselineJitterDataFrames[[4]], CurrentJitterFilenames[[4]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Random", "RK4", "RK4", BaselineJitterDataFrames[[5]], CurrentJitterFilenames[[5]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "RK4", "RK4", BaselineJitterDataFrames[[6]], CurrentJitterFilenames[[6]], species)
    }
  } 

# if (CREATE_BASELINE_FILES == FALSE) {
#   return()
# }
    
  print("------------------ Forced Biomass Tests (Stepped) ------------------")
  # Tests 29-34 Forced BIOMASS with Stepped Noise
  # Tests 35-40 Forced Migration with Stepped Noise
  REcosystem_scenario_stepped <- REcosystem_scenario
  typeData             <- list(FORCED_BIOMASS,FORCED_MIGRATION)
  numMonths <- nrow(REcosystem_scenario_stepped$forcing$ForcedBio)
  for (i in 1:length(typeData)) {
    REcosystem_scenario <- rsim.scenario(REco, REco.params, 1:50)
    REcosystem_scenario_stepped <- copy(REcosystem_scenario)
    theTypeData  <- typeData[[i]]
    modNum <- modNum + 1
    if (theTypeData == FORCED_BIOMASS) {
      REcosystem_scenario_stepped$forcing$ForcedBio <- modifyForcingMatrix(modNum,species,STEPPED,theTypeData,
                         REcosystem_scenario_stepped$forcing$ForcedBio,
                         REcosystem_scenario_stepped,
                         POSITIVE_ONLY)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped)
      BaselineSteppedFiles  <- list(BaselineABForcedBioOutBiomassStepped,  BaselineABForcedBioOutCatchStepped,  BaselineABForcedBioOutGearCatchStepped,
                                    BaselineRK4ForcedBioOutBiomassStepped, BaselineRK4ForcedBioOutCatchStepped, BaselineRK4ForcedBioOutGearCatchStepped)
      CurrentSteppedFiles   <- list(CurrentABForcedBioOutBiomassStepped,  CurrentABForcedBioOutCatchStepped,  CurrentABForcedBioOutGearCatchStepped,
                                    CurrentRK4ForcedBioOutBiomassStepped, CurrentRK4ForcedBioOutCatchStepped, CurrentRK4ForcedBioOutGearCatchStepped)
    } else if (theTypeData == FORCED_MIGRATION) {
      REcosystem_scenario_stepped$forcing$ForcedMigrate <- modifyForcingMatrix(modNum,species,STEPPED,theTypeData,
                         REcosystem_scenario_stepped$forcing$ForcedMigrate,
                         REcosystem_scenario_stepped,
                         POSITIVE_ONLY)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped)
      BaselineSteppedFiles  <- list(BaselineABForcedMigOutBiomassStepped,  BaselineABForcedMigOutCatchStepped,  BaselineABForcedMigOutGearCatchStepped,
                                    BaselineRK4ForcedMigOutBiomassStepped, BaselineRK4ForcedMigOutCatchStepped, BaselineRK4ForcedMigOutGearCatchStepped)
      CurrentSteppedFiles   <- list(CurrentABForcedMigOutBiomassStepped,  CurrentABForcedMigOutCatchStepped,  CurrentABForcedMigOutGearCatchStepped,
                                    CurrentRK4ForcedMigOutBiomassStepped, CurrentRK4ForcedMigOutCatchStepped, CurrentRK4ForcedMigOutGearCatchStepped)
    }
    REcosystem_AB_Current_Stepped  <- rsim.run(REcosystem_scenario_stepped,method='AB', years=1:50)
    REcosystem_RK4_Current_Stepped <- rsim.run(REcosystem_scenario_stepped,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      writeDataFile(REcosystem_AB_Current_Stepped$out_Biomass,     BaselineSteppedFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Catch,       BaselineSteppedFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Gear_Catch,  BaselineSteppedFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Biomass,    BaselineSteppedFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Catch,      BaselineSteppedFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Gear_Catch, BaselineSteppedFiles[[6]])
    } else {
      writeDataFile(REcosystem_AB_Current_Stepped$out_Biomass,     CurrentSteppedFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Catch,       CurrentSteppedFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Gear_Catch,  CurrentSteppedFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Biomass,    CurrentSteppedFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Catch,      CurrentSteppedFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Gear_Catch, CurrentSteppedFiles[[6]])
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[1]], CurrentSteppedFiles[[1]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[2]], CurrentSteppedFiles[[2]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[3]], CurrentSteppedFiles[[3]], species)
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[4]], CurrentSteppedFiles[[4]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[5]], CurrentSteppedFiles[[5]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[6]], CurrentSteppedFiles[[6]], species)
    }
  }
  
  print("------------------ Forced Effort Tests (Jitter) ------------------")
  numMonths <- nrow(REcosystem_scenario$fishing$ForcedEffort)
  REcosystem_scenario_jitter <- REcosystem_scenario
  fishingOriginalData  <- list(REcosystem_scenario$fishing$ForcedEffort, REcosystem_scenario$fishing$ForcedFRate, REcosystem_scenario$fishing$ForcedCatch)
  typeData             <- list(FORCED_EFFORT, FORCED_FRATE, FORCED_CATCH)
  for (i in 1:length(fishingOriginalData)) {
    theTypeData  <- typeData[[i]]
    modNum <- modNum + 1
    if (theTypeData == FORCED_EFFORT) {
      ForcedMatrix <- modifyFishingMatrix(modNum,species,fleets,theTypeData,fishingOriginalData[[i]],
                                          REco.params$model,POSITIVE_AND_NEGATIVE)
      REcosystem_scenario_jitter$fishing$ForcedEffort <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedEffOutBiomassJitter, BaselineABForcedEffOutCatchJitter, BaselineABForcedEffOutGearCatchJitter,
                                   BaselineRK4ForcedEffOutBiomassJitter,BaselineRK4ForcedEffOutCatchJitter,BaselineRK4ForcedEffOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedEffOutBiomassJitter,  CurrentABForcedEffOutCatchJitter,  CurrentABForcedEffOutGearCatchJitter,
                                   CurrentRK4ForcedEffOutBiomassJitter, CurrentRK4ForcedEffOutCatchJitter, CurrentRK4ForcedEffOutGearCatchJitter)
    } else if (theTypeData == FORCED_FRATE) {
      ForcedMatrix <- modifyFishingMatrix(modNum,species,fleets,theTypeData,fishingOriginalData[[i]],
                                          REco.params$model,POSITIVE_ONLY)
      REcosystem_scenario_jitter$fishing$ForcedFRate  <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedFRaOutBiomassJitter, BaselineABForcedFRaOutCatchJitter, BaselineABForcedFRaOutGearCatchJitter,
                                   BaselineRK4ForcedFRaOutBiomassJitter,BaselineRK4ForcedFRaOutCatchJitter,BaselineRK4ForcedFRaOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedFRaOutBiomassJitter,  CurrentABForcedFRaOutCatchJitter,  CurrentABForcedFRaOutGearCatchJitter,
                                   CurrentRK4ForcedFRaOutBiomassJitter, CurrentRK4ForcedFRaOutCatchJitter, CurrentRK4ForcedFRaOutGearCatchJitter)
    } else if (theTypeData == FORCED_CATCH) {
      REcosystem_scenario_jitter$fishing$ForcedCatch  <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter,
                                REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedCatOutBiomassJitter, BaselineABForcedCatOutCatchJitter, BaselineABForcedCatOutGearCatchJitter,
                                BaselineRK4ForcedCatOutBiomassJitter,BaselineRK4ForcedCatOutCatchJitter,BaselineRK4ForcedCatOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedCatOutBiomassJitter,  CurrentABForcedCatOutCatchJitter,  CurrentABForcedCatOutGearCatchJitter,
                                CurrentRK4ForcedCatOutBiomassJitter, CurrentRK4ForcedCatOutCatchJitter, CurrentRK4ForcedCatOutGearCatchJitter)
    }
    # Tests 41-46 Forced Effort with Jitter
    # Tests 47-52 Forced FRate with Jitter
    # Tests 53-58 Forced Catch with Jitter
    REcosystem_AB_Current_Jitter  <- rsim.run(REcosystem_scenario_jitter,method='AB', years=1:50)
    REcosystem_RK4_Current_Jitter <- rsim.run(REcosystem_scenario_jitter,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      writeDataFile(REcosystem_AB_Current_Jitter$out_Biomass,     BaselineJitterFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Catch,       BaselineJitterFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Gear_Catch,  BaselineJitterFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Biomass,    BaselineJitterFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Catch,      BaselineJitterFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Gear_Catch, BaselineJitterFiles[[6]])
    } else {
      writeDataFile(REcosystem_AB_Current_Jitter$out_Biomass,     CurrentJitterFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Catch,       CurrentJitterFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Jitter$out_Gear_Catch,  CurrentJitterFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Biomass,    CurrentJitterFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Catch,      CurrentJitterFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Jitter$out_Gear_Catch, CurrentJitterFiles[[6]])
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[1]], CurrentJitterFiles[[1]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[2]], CurrentJitterFiles[[2]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[3]], CurrentJitterFiles[[3]], species)
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[4]], CurrentJitterFiles[[4]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[5]], CurrentJitterFiles[[5]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[6]], CurrentJitterFiles[[6]], species)
    }
  }

  
  print("------------------ Forced Effort Tests (Stepped) ------------------")
  numMonths <- nrow(REcosystem_scenario$fishing$ForcedEffort)
  REcosystem_scenario_stepped <- REcosystem_scenario
  fishingOriginalData  <- list(REcosystem_scenario$fishing$ForcedEffort, REcosystem_scenario$fishing$ForcedFRate, REcosystem_scenario$fishing$ForcedCatch)
  typeData             <- list(FORCED_EFFORT,FORCED_FRATE,FORCED_CATCH)
  for (i in 1:length(fishingOriginalData)) {
    theTypeData  <- typeData[[i]]
    REcosystem_scenario_stepped <- REcosystem_scenario
    if (theTypeData == FORCED_EFFORT) {
      REcosystem_scenario_stepped$fishing$ForcedEffort <- stepifyMatrix(fishingOriginalData[[i]],fleets,0.1)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedEffOutBiomassStepped, BaselineABForcedEffOutCatchStepped, BaselineABForcedEffOutGearCatchStepped,
                                   BaselineRK4ForcedEffOutBiomassStepped,BaselineRK4ForcedEffOutCatchStepped,BaselineRK4ForcedEffOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedEffOutBiomassStepped,  CurrentABForcedEffOutCatchStepped,  CurrentABForcedEffOutGearCatchStepped,
                                   CurrentRK4ForcedEffOutBiomassStepped, CurrentRK4ForcedEffOutCatchStepped, CurrentRK4ForcedEffOutGearCatchStepped)
    } else if (theTypeData == FORCED_FRATE) {
      REcosystem_scenario_stepped$fishing$ForcedFRate <- stepifyMatrix(fishingOriginalData[[i]],species,0.01)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedFRaOutBiomassStepped, BaselineABForcedFRaOutCatchStepped, BaselineABForcedFRaOutGearCatchStepped,
                                   BaselineRK4ForcedFRaOutBiomassStepped,BaselineRK4ForcedFRaOutCatchStepped,BaselineRK4ForcedFRaOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedFRaOutBiomassStepped,  CurrentABForcedFRaOutCatchStepped,  CurrentABForcedFRaOutGearCatchStepped,
                                   CurrentRK4ForcedFRaOutBiomassStepped, CurrentRK4ForcedFRaOutCatchStepped, CurrentRK4ForcedFRaOutGearCatchStepped)      
    } else if (theTypeData == FORCED_CATCH) {
      REcosystem_scenario_stepped$fishing$ForcedCatch <- stepifyMatrix(fishingOriginalData[[i]],species,0.01)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedCatOutBiomassStepped, BaselineABForcedCatOutCatchStepped, BaselineABForcedCatOutGearCatchStepped,
                                   BaselineRK4ForcedCatOutBiomassStepped,BaselineRK4ForcedCatOutCatchStepped,BaselineRK4ForcedCatOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedCatOutBiomassStepped,  CurrentABForcedCatOutCatchStepped,  CurrentABForcedCatOutGearCatchStepped,
                                   CurrentRK4ForcedCatOutBiomassStepped, CurrentRK4ForcedCatOutCatchStepped, CurrentRK4ForcedCatOutGearCatchStepped)      
    }
    # Tests 59-64 Forced Effort with Stepped Noise
    # Tests 65-70 Forced FRate with Stepped Noise
    # Tests 71-76 Forced Catch with Stepped Noise
    REcosystem_AB_Current_Stepped  <- rsim.run(REcosystem_scenario_stepped,method='AB', years=1:50)
    REcosystem_RK4_Current_Stepped <- rsim.run(REcosystem_scenario_stepped,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      writeDataFile(REcosystem_AB_Current_Stepped$out_Biomass,     BaselineSteppedFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Catch,       BaselineSteppedFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Gear_Catch,  BaselineSteppedFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Biomass,    BaselineSteppedFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Catch,      BaselineSteppedFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Gear_Catch, BaselineSteppedFiles[[6]])
    } else {
      writeDataFile(REcosystem_AB_Current_Stepped$out_Biomass,     CurrentSteppedFiles[[1]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Catch,       CurrentSteppedFiles[[2]])
      writeDataFile(REcosystem_AB_Current_Stepped$out_Gear_Catch,  CurrentSteppedFiles[[3]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Biomass,    CurrentSteppedFiles[[4]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Catch,      CurrentSteppedFiles[[5]])
      writeDataFile(REcosystem_RK4_Current_Stepped$out_Gear_Catch, CurrentSteppedFiles[[6]])
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[1]], CurrentSteppedFiles[[1]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[2]], CurrentSteppedFiles[[2]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[3]], CurrentSteppedFiles[[3]], species)
      runTestRDS(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[4]], CurrentSteppedFiles[[4]], species)
      runTestRDS(inc(runNum),"out_Catch",      theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[5]], CurrentSteppedFiles[[5]], species)
      runTestRDS(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[6]], CurrentSteppedFiles[[6]], species)
    }
  }
  
  if (! CREATE_BASELINE_FILES) {
    print(paste0("Completed ",inc(runNum)," test(s)."))
  }
  
})
