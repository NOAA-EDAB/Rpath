library(Rpath)
library(qpdf)
library(testthat)
library(stringr)
library(distillery)
library(ggplot2)
library(ggpubr)
library(rlist)
# library(plotly)
# library(htmlwidgets)

data(package="Rpath")

# ---- Modify this toggle to TRUE to generate the baseline files. ----
# ---- Reset it back to FALSE to run the tests. ----------------------
# CREATE_BASELINE_FILES <- TRUE
CREATE_BASELINE_FILES <- FALSE

NUMBER_OF_STEPS <- 5 # should be an odd multiple of nrows=600 (i.e., 5,15,30)
FACTOR <- 5
SEED   <- 1
SEED_OFFSET <- 1000
TOLERANCE <- 1e-5
RUN_QUIET <- TRUE
YLIMIT_DIFFERENCE_PLOTS <- 0.05
PLOT_TYPE <- 1 # 1 = Baseline and Current superimposed, 2 = difference of (Current-Baseline)
PLOT_SHOW <- 1 # 1 - All Plots, 2 = Only plots reflecting test errors # Not sure if can be implemented
INPUT_DATA_DIR_BASELINE  <- 'data/input/baseline'
INPUT_DATA_DIR_CURRENT   <- here::here('tests/testthat/data/input')
INPUT_DATA_DIR_BASELINE  <- INPUT_DATA_DIR_CURRENT # RSK
OUTPUT_DATA_DIR          <- here::here('tests/testthat/data/output')

print(paste0("OUTPUT_DATA_DIR: ",         OUTPUT_DATA_DIR)) #RSK
print(paste0("INPUT_DATA_DIR_BASELINE: ", INPUT_DATA_DIR_BASELINE)) #RSK
print(paste0("INPUT_DATA_DIR_CURRENT: ",  INPUT_DATA_DIR_CURRENT)) #RSK

# Create the current and output directories if they don't already exist.
if (! dir.exists(INPUT_DATA_DIR_CURRENT)) {
  dir.create(INPUT_DATA_DIR_CURRENT,recursive=TRUE)
}
if (! dir.exists(OUTPUT_DATA_DIR)) {
  dir.create(OUTPUT_DATA_DIR,recursive=TRUE)
}

#' Stepify Effort
#' 
#' This function stepifies a vector of effort data. Stepification consists of modifying (i.e., multiplying) the
#' vector data in question by a series of constant values (i.e., steps). Plotting these steps
#' shows a stair-stepped function of various forms based upon the passed type. Type 1 = mountain
#' shape, Type 2 = valley shape, Type 3 = mixed shape.
#'
#' @param effort : vector effort data
#' @param type : type of "stepification"
#'
#' @return Returns the "stepped" vector data
#' 
stepifyVector <- function(vectorData,type=1) {
  steps     <- getSteps(type)
  numMonths <- length(vectorData)
  incMonths <- numMonths/ NUMBER_OF_STEPS;
  inc <- 0
  for (i in 1:NUMBER_OF_STEPS) {
    start <- inc + 1
    start <- if (start == 0) 1 else start
    end   <- inc+incMonths
    vectorData[start:end] <- vectorData[start:end] * steps[i]
    inc   <- inc + incMonths
  } 
  return(vectorData)
}

#' Stepify a Matrix
#' 
#' This function sequences through the appropriately named columns and calls stepifyVector
#' to stepify those columns.
#' 
#' @family Rpath functions
#'
#' @param forcedEffortMatrix : matrix in which specific columns will be "stepified"
#' @param memberNames : vector of either species or fleet names
#'
#' @return Returns the updated forced effort matrix
#'
stepifyMatrix <- function(forcedEffortMatrix,memberNames) {
  ForcedMatrix <- forcedEffortMatrix
  for (i in 1:length(memberNames)) {
    memberName <- memberNames[i]
    memberCol <- ForcedMatrix[,memberName]
    memberStepped <- stepifyVector(memberCol,i)
    ForcedMatrix[,memberName] <- memberStepped
    # plot(memberStepped,type='l',lwd=5,xlab="Months",ylab="Effort",main=paste0("Forced Effort with Stepped Noise - ",memberName))
  }
  return(ForcedMatrix)
}

#' Get stepification steps
#' 
#' This function returns the step offsets in order to give a function a stair-stepped pattern.
#'
#' @param type : The type of stair-stepped pattern. There are currently 3 types: 1 (mountain pattern), 2 (valley pattern), and 3 (mix of mountain and valley pattern)
#'
#' @return Returns a vector of stair-stepped offsets
#'
getSteps <- function(type=1) {
  STEPS_WIDTH <- 0.1 # approx. +/- 10% of initial value
  steps = c()
  tog <- -1
  min <- 1.0
# max <- if (type == 2 || type == 3) 0.5 else 2.0
  max <- if (type == 2 || type == 3) (1.0-STEPS_WIDTH) else (1.0+STEPS_WIDTH)
  if (is.even(NUMBER_OF_STEPS)) {
     NUMBER_OF_STEPS <- NUMBER_OF_STEPS + 1
  }
  halfway   <- NUMBER_OF_STEPS %/% 2
  increment <- (max-min)/halfway

  for (i in 1:NUMBER_OF_STEPS) {
    if (i <= halfway+1) {
      inc <- increment
      start <- min
      ii <- i-1
    } else {
      inc <- -increment
      start <- max-inc
      ii <- i-halfway
    }
    tog <- if (type == 3) -tog else 1.0
    steps[i] <- abs(start + tog*ii*inc)
  }
  return(steps)
}

#' Stepify Biomass
#'
#' Takes as input a biomass value and stepifies the value by dividing the number of months
#' by the number of steps desired and multiplying each interval by a different step value.
#'
#' @param numMonths : number of months in model
#' @param biomass : biomass value to stepify across months
#' @param type : type of stepification desired (1, 2, or 3)
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot
#' @param title : title for plot
#'
#' @return Returns the vector of stepified segments
#' 
stepifyBiomass <- function(numMonths,biomass,type=1,xlabel,ylabel,title) {
  steps <- getSteps(type)
  incMonths <- numMonths / NUMBER_OF_STEPS;
  parts <- c()
  for (i in 1:NUMBER_OF_STEPS) {
    part  <- replicate(incMonths, steps[i]*biomass)
    parts <- c(parts,part)
  }
  # plot(parts,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)
  return(parts)
}

#' Add Jitter (i.e., random noise)
#' 
#' Adds random noise to the specified matrix.
#'
#' @param matrix : data matrix to be jittered
#' @param factor : the R jitter factor (the larger the number the greater the jitter amount)
#' @param seedOffset : offset to the seed value, useful for getting "groups" of seed values
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot 
#' @param title : main title for plot
#'
#' @return Returns the jittered matrix
#' 
addJitter <- function(matrix,seedOffset,xlabel,ylabel,title) {
  set.seed(seedOffset) # *SEED)
  # From jitter() doc: If amount == 0, jitter returns factor * z/50, where
  # z = max(x0) - min(x), aka the range. So if factor=5 and amount=0, jitter()
  # returns a random value within a tenth of the range.
  jitteredMatrix <- jitter(matrix,factor=FACTOR,amount=0)
  if (xlabel != '' && ylabel != '' & title != '') {
    # plot(jitteredMatrix,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)    
  } 
  return(jitteredMatrix)
}

#' Create a jittered vector
#' 
#' Used for jittering a matrix column. Can't use replicate because need to pass a different seed for every value for reproducibility.
#'
#' @param value : value to add jitter to
#' @param numElements : number of elements in vector
#' @param seedOffset : offset to the main seed value
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot
#' @param title : main title for plot
#'
#' @return Returns the column-jittered matrix
#' 
createJitterVectorFromValue <- function(value,numElements,seedOffset,xlabel,ylabel,title) {
  jitterVector <- c()
  for (i in 1:numElements) {
    jitterVector <- append(jitterVector,addJitter(value,seedOffset+i,'','',''))
    # currentSeed  <- seedOffset*SEED + i
  }
# plot(jitterVector,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)
  return(jitterVector)
}

#' Plot test results superimposed on the same plot
#' 
#' Plot specific columns from that passed in matrices and plot them as groups of plots per page.
#' The plots will be of the baseline run and the current run plot superimposed on the same plot.
#'
#' @param BaseData : Baseline data to compare current data to
#' @param CurrData : Current data to compare against pre-run baseline data
#' @param baseAlg : The baseline algorithm used (currently AB or RK4)
#' @param currAlg : The current algorithm used (currently AB or RK4)
#' @param tableName : Table name used for the plot's title
#' @param forcedData : The data that's being forced (i.e., Biomass, Catch)
#' @param forcedType : The type of forcing being done (.e., Random (Jitter) or Stair-Stepped)
#' @param species : A vector of fish species
#'
#' @return Returns the final, combined plot
#' 
plotResultsSuperimposed <- function(BaseData,CurrData,baseAlg,currAlg,tableName,forcedData,forcedType,species) {
  plots  <- list()
  group  <- species
  yLabel <- "Biomass (mt/km²)"
  currDf <- data.frame()
  baseDf <- data.frame()

  for (member in group) {
    xvalues     <- c(1:length(CurrData[,member]))
    numMonths   <- length(CurrData[,member])
    currYvalues <- CurrData[,member]
    baseYvalues <- BaseData[,member]
    currDf    <- data.frame(Legend=replicate(numMonths,paste0('Current ', currAlg)), xvalues,currYvalues)
    baseDf    <- data.frame(Legend=replicate(numMonths,paste0('Baseline ',baseAlg)), xvalues,baseYvalues)

    aPlot <- ggplot() + 
      geom_line(data=baseDf, aes(x=xvalues, y=baseYvalues, color=Legend)) +
      geom_line(data=currDf, aes(x=xvalues, y=currYvalues, color=Legend)) +
      labs(x="Months",y=yLabel,
           title=paste0('Sim Run (',forcedData,' w/ ',forcedType,' Noise) - ',member),
           subtitle=paste0("Dataset: ",tableName)) +
      scale_color_manual(name="Legend:",values=c('red','darkblue')) +
      theme( plot.title = element_text(hjust=0.5,size=7,face="bold"),
             plot.subtitle = element_text(hjust=0.5,size=7),
             axis.text = element_text(size=8),
             axis.title = element_text(size=8),
             legend.text = element_text(size=8),
             legend.title = element_text(size=10),
             legend.position='bottom',
             legend.spacing.y = unit(0.0,'cm'),
             legend.background = element_rect(fill='#f7f7f7'),
             # legend.box.background = element_rect(color = 'black'),
             plot.background = element_rect(color='black',fill=NA,linewidth=1)
           )
    plots <- list.append(plots,aPlot)
    
  }
  combinedPlot <- ggarrange(plotlist=plots,nrow=3,ncol=2)
  # annotate_figure(combinedPlot, top = text_grob("Sample main title here", color = "red", face = "bold", size = 14))
   # saveWidget(ggplotly(combinedPlot), file = "Rplots.html");
   # print(ggplotly(combinedPlot))
   print(combinedPlot)
}

#' Plot the difference of the two runs
#' 
#' Plot specific columns from that passed in matrices and plot them as groups of plots per page.
#' The plots will be of the current run - baseline run. So if there's a perfect match, the plot should
#' be all zeros.
#'
#' @param BaseData : Baseline data to compare current data to
#' @param CurrData : Current data to compare against pre-run baseline data
#' @param baseAlg : The baseline algorithm used (currently AB or RK4)
#' @param currAlg : The current algorithm used (currently AB or RK4)
#' @param tableName : Table name used for the plot's title
#' @param forcedData : The data that's being forced (i.e., Biomass, Catch)
#' @param forcedType : The type of forcing being done (.e., Random (Jitter) or Stair-Stepped)
#' @param species : A vector of fish species
#'
#' @return Returns the final, combined plot
#' 
plotResultsDifference <- function(BaseData,CurrData,baseAlg,currAlg,tableName,forcedData,forcedType,species) {
  plots  <- list()
  group  <- species
  yLabel <- "Biomass (mt/km²)"
  diffDf <- data.frame()
  
  for (member in group) {
    xvalues     <- c(1:length(CurrData[,member]))
    numMonths   <- length(CurrData[,member])
    currYvalues <- CurrData[,member]
    baseYvalues <- BaseData[,member]
    diffDf <- data.frame(Legend=replicate(numMonths,paste0('Current(',currAlg,')-Baseline(',baseAlg,') ')), xvalues,currYvalues-baseYvalues)
    aPlot  <- ggplot() + 
      geom_line(data=diffDf, aes(x=xvalues, y=currYvalues-baseYvalues, color=Legend)) +
      labs(x="Months",y=yLabel,
           title=paste0('Sim Run using ',forcedData,' with ',forcedType,' Noise - ',member),
           subtitle=paste0("Dataset: ",tableName)) +
      scale_color_manual(name="Legend:",values=c('darkblue')) +
      coord_cartesian(ylim = c(-YLIMIT_DIFFERENCE_PLOTS, YLIMIT_DIFFERENCE_PLOTS)) + # RSK
      theme( plot.title = element_text(hjust=0.5,size=7,face="bold"),
             plot.subtitle = element_text(hjust=0.5,size=7),
             axis.text = element_text(size=8),
             axis.title = element_text(size=8),
             legend.text = element_text(size=8),
             legend.title = element_text(size=10),
             legend.position='bottom',
             legend.spacing.y = unit(0.0,'cm'),
             legend.background = element_rect(fill='#f7f7f7'),
             # legend.box.background = element_rect(color = 'black'),
             plot.background = element_rect(color='black',fill=NA,linewidth=1)
      )
    plots <- list.append(plots,aPlot)
  }
  combinedPlot <- ggarrange(plotlist=plots,ncol=1)
  # annotate_figure(combinedPlot, top = text_grob("Sample main title here", color = "red", face = "bold", size = 14))
  print(combinedPlot)
} 

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
  testthat::expect_equal(baselineTable,currentTable,tolerance=TOLERANCE)
}

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
#' @param outputTable : table that represents the current run's data
#' @param outputFile : The output file to contain the current run's data
#' @param species : A vector of fish species
#'
#' @return No return value
#' 
runTest <- function(runNum,tableName,forcedData,forcedType,baseAlg,currAlg,baselineTable,outputTable,outputFile,species) {
 
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
  if (tableName == 'out_Biomass' || tableName == 'out_Catch') { # || tableName == 'out_Gear_Catch') {
    if (PLOT_TYPE == 1) {
      plotResultsSuperimposed(baselineTable, outputTable, baseAlg, currAlg, tableName, forcedData, forcedType, species)
    } else {
      plotResultsDifference(baselineTable, outputTable, baseAlg, currAlg, tableName, forcedData, forcedType, species)
    }
  }

  write.table(outputTable, file=outputFile) 
  inputTable <- read.table(outputFile, fill = TRUE,sep = " ")
  retv <- testthat::expect_equal(baselineTable,inputTable,tolerance=TOLERANCE)

  # Write out the difference table (current-baseline)
  diffTable <- (outputTable-baselineTable)
  write.table(diffTable, file=file.path(OUTPUT_DATA_DIR,paste0("test_",paddedRunNum,".dat")))
}

#' Modify a scene$fishing matrix
#' 
#' Modifies the appropriate columns from a scene$fishing matrix.
#'
#' @param species : A vector of fish species
#' @param fleets : A vector of fleets
#' @param typeData : The type of forced data (i.e., Forced Effort, Forced FRate, Forced Catch)
#' @param forcingData : The forced data matrix
#'
#' @return Returns a matrix that's been updated withe the forced data
#' 
modifyFishingMatrix <- function(modNum,species,fleets,typeData,forcingData) {
  ForcedMatrix <- forcingData
  items <- c()
  if (typeData == "Forced Effort") {
    items <- fleets
  } else if (typeData == "Forced FRate" || typeData == "Forced Catch") {
    items <- species 
  } else {
    print(paste0("Error: Found invalid typeData of: ",typeData))
    return(ForcedMatrix)
  }

  for (i in 1:length(items)) {
    item <- items[i]
    matrixData           <- ForcedMatrix[,item]
    matrixDataWithJitter <- addJitter(matrixData,modNum*SEED_OFFSET*SEED+i,"Months","Effort",paste0(typeData," with Random Noise - ",item))
    ForcedMatrix[,item]  <- matrixDataWithJitter
  }
  return(ForcedMatrix)
}

#' Modify a scene$forcing matrix
#' 
#' Modifies the appropriate columns from a scene$forcing matrix.
#'
#' @param species : A vector of fish species
#' @param modifyType : The type of modification (i.e., Jittered or Stepped)
#' @param typeData : The type of forced data (i.e., Forced Bio, Forced Migrate)
#' @param forcingData : The forced data matrix
#' @param scene : The REcosystem_scene object 
#'
#' @return Returns a matrix that's been updated withe the forced data
#' 
modifyForcingMatrix <- function (modNum,species,modifyType,typeData,forcingData,scene) {
  ForcedMatrix <- forcingData
  numMonths <- nrow(ForcedMatrix)
  if (typeData == "Forced Bio" || typeData == "Forced Migrate") {
    for (i in 1:length(species)) {
      aSpecies <- species[[i]]
      speciesBiomass <- scene$start_state$Biomass[aSpecies]
      if (modifyType == 'Jittered') {
        ForcedMatrix[,aSpecies] <- createJitterVectorFromValue(speciesBiomass,numMonths,modNum*i*SEED_OFFSET, "Months","Biomass (mt/km²)",paste0(typeData,' with ',modifyType,' Noise - ',aSpecies))
      } else {
        stepType <- ((i-1)%%3)+1 # Only current step types are 1, 2, or 3
        ForcedMatrix[,aSpecies] <- stepifyBiomass(numMonths,speciesBiomass,stepType,"Months","Biomass (mt/km²)",paste0(typeData,' with ',modifyType,' Noise - ',aSpecies))
      }
    }
  }

  return(ForcedMatrix)
}

# This is an increment function
inc <- function(value) eval.parent(substitute(value <- value + 1))


testthat::test_that("Rpath Unit Tests", {
  BaselineJitterTables  <- list()
  BaselineJitterFiles   <- list()
  CurrentJitterFiles    <- list()
  BaselineSteppedTables <- list()
  BaselineSteppedFiles  <- list()
  CurrentSteppedFiles   <- list()
  fleets  <- c('Trawlers','Midwater','Dredgers')
  species <- c('OtherGroundfish','Megabenthos','Seals','JuvRoundfish1','AduRoundfish1')
  originalWorkingDir <- getwd();
  modNum <- 1
  runNum <- 0
  
  # ---------- Set up initial file paths ----------
  # N.B. The Baseline and Current AB and RK4 files are .csv files since they were produce by
  # the write.rsim() function and not the more generic write.table() function.
  BaselineRpathObjTopLevel                 <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RpathObj_TopLevel.dat')
  BaselineRpathObjSummary                  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RpathObj_Summary.dat')
  BaselineAB                               <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB.csv')
  BaselineRK4                              <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4.csv')
  BaselineABOutBiomass                     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutBiomass.dat')
  BaselineRK4OutBiomass                    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutBiomass.dat')
  BaselineABOutCatch                       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutCatch.dat')
  BaselineRK4OutCatch                      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutCatch.dat')
  BaselineABOutGearCatch                   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_OutGearCatch.dat')
  BaselineRK4OutGearCatch                  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_OutGearCatch.dat')
  BaselineABForcedBioOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped.dat')
  BaselineRK4ForcedBioOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped.dat')
  BaselineABForcedBioOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped.dat')
  BaselineRK4ForcedBioOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped.dat')
  BaselineABForcedBioOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped.dat')
  BaselineRK4ForcedBioOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped.dat')
  BaselineABForcedBioOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter.dat')
  BaselineRK4ForcedBioOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter.dat')
  BaselineABForcedBioOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter.dat')
  BaselineRK4ForcedBioOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter.dat')
  BaselineABForcedBioOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter.dat')
  BaselineRK4ForcedBioOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter.dat')
  BaselineABForcedMigOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped.dat')
  BaselineRK4ForcedMigOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped.dat')
  BaselineABForcedMigOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped.dat')
  BaselineRK4ForcedMigOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped.dat')
  BaselineABForcedMigOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped.dat')
  BaselineRK4ForcedMigOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped.dat')
  BaselineABForcedMigOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter.dat')
  BaselineRK4ForcedMigOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter.dat')
  BaselineABForcedMigOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter.dat')
  BaselineRK4ForcedMigOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter.dat')
  BaselineABForcedMigOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter.dat')
  BaselineRK4ForcedMigOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter.dat')
  BaselineABForcedEffOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped.dat')
  BaselineABForcedEffOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped.dat')
  BaselineABForcedEffOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped.dat')
  BaselineABForcedEffOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter.dat')
  BaselineABForcedEffOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter.dat')
  BaselineABForcedEffOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter.dat')
  BaselineABForcedFRaOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped.dat')
  BaselineABForcedFRaOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped.dat')
  BaselineABForcedFRaOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped.dat')
  BaselineABForcedFRaOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter.dat')
  BaselineABForcedFRaOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter.dat')
  BaselineABForcedFRaOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter.dat')
  BaselineABForcedCatOutBiomassStepped     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped.dat')
  BaselineABForcedCatOutCatchStepped       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped.dat')
  BaselineABForcedCatOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped.dat')
  BaselineABForcedCatOutBiomassJitter      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter.dat')
  BaselineABForcedCatOutCatchJitter        <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter.dat')
  BaselineABForcedCatOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter.dat')
  BaselineRK4ForcedEffOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped.dat')
  BaselineRK4ForcedEffOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped.dat')
  BaselineRK4ForcedEffOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped.dat')
  BaselineRK4ForcedEffOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter.dat')
  BaselineRK4ForcedEffOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter.dat')
  BaselineRK4ForcedEffOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter.dat')
  BaselineRK4ForcedFRaOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped.dat')
  BaselineRK4ForcedFRaOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped.dat')
  BaselineRK4ForcedFRaOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped.dat')
  BaselineRK4ForcedFRaOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter.dat')
  BaselineRK4ForcedFRaOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter.dat')
  BaselineRK4ForcedFRaOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter.dat')
  BaselineRK4ForcedCatOutBiomassStepped    <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped.dat')
  BaselineRK4ForcedCatOutCatchStepped      <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped.dat')
  BaselineRK4ForcedCatOutGearCatchStepped  <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped.dat')
  BaselineRK4ForcedCatOutBiomassJitter     <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter.dat')
  BaselineRK4ForcedCatOutCatchJitter       <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter.dat')
  BaselineRK4ForcedCatOutGearCatchJitter   <- file.path(INPUT_DATA_DIR_BASELINE,'REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter.dat')
  #
  CurrentRpathObjTopLevel                  <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RpathObj_TopLevel.dat')
  CurrentRpathObjSummary                   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RpathObj_Summary.dat')
  CurrentAB                                <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB.csv')
  CurrentRK4                               <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4.csv')
  CurrentABOutBiomass                      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutBiomass.dat')
  CurrentRK4OutBiomass                     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutBiomass.dat')
  CurrentABOutCatch                        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutCatch.dat')
  CurrentRK4OutCatch                       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutCatch.dat')
  CurrentABOutGearCatch                    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_OutGearCatch.dat')
  CurrentRK4OutGearCatch                   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_OutGearCatch.dat')
  CurrentABForcedBioOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutBiomass_Stepped.dat')
  CurrentRK4ForcedBioOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutBiomass_Stepped.dat')
  CurrentABForcedBioOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutCatch_Stepped.dat')
  CurrentRK4ForcedBioOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutCatch_Stepped.dat')
  CurrentABForcedBioOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutGearCatch_Stepped.dat')
  CurrentRK4ForcedBioOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutGearCatch_Stepped.dat')
  CurrentABForcedBioOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutBiomass_Jitter.dat')
  CurrentABForcedBioOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutCatch_Jitter.dat')
  CurrentABForcedBioOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedBio_OutGearCatch_Jitter.dat')
  CurrentRK4ForcedBioOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutBiomass_Jitter.dat')
  CurrentRK4ForcedBioOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutCatch_Jitter.dat')
  CurrentRK4ForcedBioOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedBio_OutGearCatch_Jitter.dat')
  CurrentABForcedMigOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutBiomass_Stepped.dat')
  CurrentRK4ForcedMigOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutBiomass_Stepped.dat')
  CurrentABForcedMigOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutCatch_Stepped.dat')
  CurrentRK4ForcedMigOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutCatch_Stepped.dat')
  CurrentABForcedMigOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutGearCatch_Stepped.dat')
  CurrentRK4ForcedMigOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutGearCatch_Stepped.dat')
  CurrentABForcedMigOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutBiomass_Jitter.dat')
  CurrentABForcedMigOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutCatch_Jitter.dat')
  CurrentABForcedMigOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedMig_OutGearCatch_Jitter.dat')
  CurrentRK4ForcedMigOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutBiomass_Jitter.dat')
  CurrentRK4ForcedMigOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutCatch_Jitter.dat')
  CurrentRK4ForcedMigOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedMig_OutGearCatch_Jitter.dat')
  CurrentABForcedEffOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutBiomass_Stepped.dat')
  CurrentRK4ForcedEffOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutBiomass_Stepped.dat')
  CurrentABForcedEffOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutCatch_Stepped.dat')
  CurrentRK4ForcedEffOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutCatch_Stepped.dat')
  CurrentABForcedEffOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutGearCatch_Stepped.dat')
  CurrentRK4ForcedEffOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutGearCatch_Stepped.dat')
  CurrentABForcedFRaOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutBiomass_Stepped.dat')
  CurrentRK4ForcedFRaOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutBiomass_Stepped.dat')
  CurrentABForcedFRaOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutCatch_Stepped.dat')
  CurrentRK4ForcedFRaOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutCatch_Stepped.dat')
  CurrentABForcedFRaOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutGearCatch_Stepped.dat')
  CurrentRK4ForcedFRaOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutGearCatch_Stepped.dat')
  CurrentABForcedCatOutBiomassStepped      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutBiomass_Stepped.dat')
  CurrentRK4ForcedCatOutBiomassStepped     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutBiomass_Stepped.dat')
  CurrentABForcedCatOutCatchStepped        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutCatch_Stepped.dat')
  CurrentRK4ForcedCatOutCatchStepped       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutCatch_Stepped.dat')
  CurrentABForcedCatOutGearCatchStepped    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutGearCatch_Stepped.dat')
  CurrentRK4ForcedCatOutGearCatchStepped   <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutGearCatch_Stepped.dat')
  CurrentABForcedEffOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutBiomass_Jitter.dat')
  CurrentABForcedEffOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutCatch_Jitter.dat')
  CurrentABForcedEffOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedEff_OutGearCatch_Jitter.dat')
  CurrentRK4ForcedEffOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutBiomass_Jitter.dat')
  CurrentRK4ForcedEffOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutCatch_Jitter.dat')
  CurrentRK4ForcedEffOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedEff_OutGearCatch_Jitter.dat')
  CurrentABForcedFRaOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutBiomass_Jitter.dat')
  CurrentABForcedFRaOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutCatch_Jitter.dat')
  CurrentABForcedFRaOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedFRa_OutGearCatch_Jitter.dat')
  CurrentRK4ForcedFRaOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutBiomass_Jitter.dat')
  CurrentRK4ForcedFRaOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutCatch_Jitter.dat')
  CurrentRK4ForcedFRaOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedFRa_OutGearCatch_Jitter.dat')
  CurrentABForcedCatOutBiomassJitter       <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutBiomass_Jitter.dat')
  CurrentABForcedCatOutCatchJitter         <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutCatch_Jitter.dat')
  CurrentABForcedCatOutGearCatchJitter     <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_AB_ForcedCat_OutGearCatch_Jitter.dat')
  CurrentRK4ForcedCatOutBiomassJitter      <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutBiomass_Jitter.dat')
  CurrentRK4ForcedCatOutCatchJitter        <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutCatch_Jitter.dat')
  CurrentRK4ForcedCatOutGearCatchJitter    <- file.path(INPUT_DATA_DIR_CURRENT,'REcosystem_Current_RK4_ForcedCat_OutGearCatch_Jitter.dat')

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
    REcosystemBaseline                                     <- read.table(BaselineRpathObjTopLevel,               fill = TRUE, sep = " ")
    REcosystemBaselineSummary                              <- read.table(BaselineRpathObjSummary,                fill = TRUE)
    REcosystem_Baseline_AB                                 <- read.csv(BaselineAB)
    REcosystem_Baseline_RK4                                <- read.csv(BaselineRK4)
    REcosystem_Baseline_AB_OutBiomass                      <- read.table(BaselineABOutBiomass,                   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_OutCatch                        <- read.table(BaselineABOutCatch,                     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_OutGearCatch                    <- read.table(BaselineABOutGearCatch,                 fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutBiomass                     <- read.table(BaselineRK4OutBiomass,                  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutCatch                       <- read.table(BaselineRK4OutCatch,                    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_OutGearCatch                   <- read.table(BaselineRK4OutGearCatch,                fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped    <- read.table(BaselineABForcedBioOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedBioOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped      <- read.table(BaselineABForcedBioOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped     <- read.table(BaselineRK4ForcedBioOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped  <- read.table(BaselineABForcedBioOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedBioOutGearCatchStepped,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter     <- read.table(BaselineABForcedBioOutBiomassJitter,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter    <- read.table(BaselineRK4ForcedBioOutBiomassJitter,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter       <- read.table(BaselineABForcedBioOutCatchJitter,      fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter      <- read.table(BaselineRK4ForcedBioOutCatchJitter,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter   <- read.table(BaselineABForcedBioOutGearCatchJitter,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter  <- read.table(BaselineRK4ForcedBioOutGearCatchJitter, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped    <- read.table(BaselineABForcedMigOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedMigOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped      <- read.table(BaselineABForcedMigOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped     <- read.table(BaselineRK4ForcedMigOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped  <- read.table(BaselineABForcedMigOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedMigOutGearCatchStepped,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter     <- read.table(BaselineABForcedMigOutBiomassJitter,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter    <- read.table(BaselineRK4ForcedMigOutBiomassJitter,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter       <- read.table(BaselineABForcedMigOutCatchJitter,      fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter      <- read.table(BaselineRK4ForcedMigOutCatchJitter,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter   <- read.table(BaselineABForcedMigOutGearCatchJitter,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter  <- read.table(BaselineRK4ForcedMigOutGearCatchJitter, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- read.table(BaselineABForcedEffOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedEffOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- read.table(BaselineABForcedEffOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- read.table(BaselineRK4ForcedEffOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- read.table(BaselineABForcedEffOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedEffOutGearCatchStepped,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter     <- read.table(BaselineABForcedEffOutBiomassJitter,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter       <- read.table(BaselineABForcedEffOutCatchJitter,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter   <- read.table(BaselineABForcedEffOutGearCatchJitter,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter    <- read.table(BaselineRK4ForcedEffOutBiomassJitter,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter      <- read.table(BaselineRK4ForcedEffOutCatchJitter,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter  <- read.table(BaselineRK4ForcedEffOutGearCatchJitter, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter     <- read.table(BaselineABForcedFRaOutBiomassJitter,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter       <- read.table(BaselineABForcedFRaOutCatchJitter,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter   <- read.table(BaselineABForcedFRaOutGearCatchJitter,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter    <- read.table(BaselineRK4ForcedFRaOutBiomassJitter,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter      <- read.table(BaselineRK4ForcedFRaOutCatchJitter,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter  <- read.table(BaselineRK4ForcedFRaOutGearCatchJitter, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter     <- read.table(BaselineABForcedCatOutBiomassJitter,    fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter       <- read.table(BaselineABForcedCatOutCatchJitter,      fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter   <- read.table(BaselineABForcedCatOutGearCatchJitter,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter    <- read.table(BaselineRK4ForcedCatOutBiomassJitter,   fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter      <- read.table(BaselineRK4ForcedCatOutCatchJitter,     fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter  <- read.table(BaselineRK4ForcedCatOutGearCatchJitter, fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped    <- read.table(BaselineABForcedEffOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped      <- read.table(BaselineABForcedEffOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped  <- read.table(BaselineABForcedEffOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedEffOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped     <- read.table(BaselineRK4ForcedEffOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedEffOutGearCatchStepped,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped    <- read.table(BaselineABForcedFRaOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped      <- read.table(BaselineABForcedFRaOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped  <- read.table(BaselineABForcedFRaOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedFRaOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped     <- read.table(BaselineRK4ForcedFRaOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedFRaOutGearCatchStepped,fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped    <- read.table(BaselineABForcedCatOutBiomassStepped,   fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped      <- read.table(BaselineABForcedCatOutCatchStepped,     fill = TRUE, sep = " ")
    REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped  <- read.table(BaselineABForcedCatOutGearCatchStepped, fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped   <- read.table(BaselineRK4ForcedCatOutBiomassStepped,  fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped     <- read.table(BaselineRK4ForcedCatOutCatchStepped,    fill = TRUE, sep = " ")
    REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped <- read.table(BaselineRK4ForcedCatOutGearCatchStepped,fill = TRUE, sep = " ")
  }
  
  # 
  # Save current Rpath run data
  REco <- rpath(REco.params, eco.name = 'R Ecosystem') # an Rpath object
  RpathObjTopLevel <- if (CREATE_BASELINE_FILES) BaselineRpathObjTopLevel else CurrentRpathObjTopLevel
  sink(RpathObjTopLevel)
    # capture.output(print(REco),file=RpathObjTopLevel)
    print(REco)
  sink()
  
  # Save current Rpath run summary data
  setwd(originalWorkingDir)
  RpathObjSummary <- if (CREATE_BASELINE_FILES) BaselineRpathObjSummary else CurrentRpathObjSummary
  sink(RpathObjSummary)
    cat(summary(REco))
  sink()

  # Save current Rpath sim run data for AB and RK4
  setwd(originalWorkingDir)
  REcosystem_scene <- rsim.scenario(REco, REco.params, 1:50)
  REcosystem_Current_AB_from_Sim  <- rsim.run(REcosystem_scene,method='AB', years=1:50)
  REcosystem_Current_RK4_from_Sim <- rsim.run(REcosystem_scene,method='RK4',years=1:50)
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
    files <- dir(path=file.path(cwd,OUTPUT_DATA_DIR),pattern='test_*')
    file.remove(file.path(OUTPUT_DATA_DIR,files))

    # Test 1 - Test if Balanced (i.e., "Status: Balanced" is the 2nd line of the Summary file)
    headerSummaryLines <- readLines(CurrentRpathObjSummary,n=2)
    parts <- unlist(strsplit(str_trim(headerSummaryLines[2]),split=" "))
    runTestEqual(inc(runNum),"","Is model balanced?",parts[2],"Balanced")

    # Test 2 - Test if function runs silently (i.e., no messages, warnings, or print statements)
    runTestSilent(inc(runNum),"Does model run without any terminal output (i.e., warnings, errors)?",REco.params,'R Ecosystem')
    
    # Test 3 - Test that the REcosystem object is the same as the saved original REcosystem object
    REcosystemCurrent <- read.table(CurrentRpathObjTopLevel,fill = TRUE, sep = " ")
    runTestEqual(inc(runNum),"","Is the baseline Rpath object equivalent to the current Rpath object (toplevel data)?",REcosystemBaseline,REcosystemCurrent)

    # Test 4 - Test that the REcosystem Summary is the same as the saved original REcosystem Summary
    REcosystemSummaryCurrent <- read.table(CurrentRpathObjSummary,fill = TRUE)
    runTestEqual(inc(runNum),"","Is the baseline Rpath run Summary the same as the current Rpath Summary?",REcosystemBaselineSummary,REcosystemSummaryCurrent)

    # Tests 5-16 - Test that REcosystem AB object is same as RK4 object with no perturbations
    REcosystem_Current_AB  <- read.csv(CurrentAB)
    REcosystem_Current_RK4 <- read.csv(CurrentRK4)
  }

  if (CREATE_BASELINE_FILES) { 
    write.table(REcosystem_Current_AB_from_Sim$out_Biomass,    file=BaselineABOutBiomass)
    write.table(REcosystem_Current_RK4_from_Sim$out_Biomass,   file=BaselineRK4OutBiomass)
    write.table(REcosystem_Current_AB_from_Sim$out_Catch,      file=BaselineABOutCatch)
    write.table(REcosystem_Current_RK4_from_Sim$out_Catch,     file=BaselineRK4OutCatch)
    write.table(REcosystem_Current_AB_from_Sim$out_Gear_Catch, file=BaselineABOutGearCatch)
    write.table(REcosystem_Current_RK4_from_Sim$out_Gear_Catch,file=BaselineRK4OutGearCatch)  
  } else {
    write.table(REcosystem_Current_AB_from_Sim$out_Biomass,    file=CurrentABOutBiomass)
    write.table(REcosystem_Current_RK4_from_Sim$out_Biomass,   file=CurrentRK4OutBiomass)
    write.table(REcosystem_Current_AB_from_Sim$out_Catch,      file=CurrentABOutCatch)
    write.table(REcosystem_Current_RK4_from_Sim$out_Catch,     file=CurrentRK4OutCatch)
    write.table(REcosystem_Current_AB_from_Sim$out_Gear_Catch, file=CurrentABOutGearCatch)
    write.table(REcosystem_Current_RK4_from_Sim$out_Gear_Catch,file=CurrentRK4OutGearCatch)
    REcosystem_Current_AB_OutBiomass    <- read.table(CurrentABOutBiomass,   fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutBiomass   <- read.table(CurrentRK4OutBiomass,  fill=TRUE,sep=" ")
    REcosystem_Current_AB_OutCatch      <- read.table(CurrentABOutCatch,     fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutCatch     <- read.table(CurrentRK4OutCatch,    fill=TRUE,sep=" ")
    REcosystem_Current_AB_OutGearCatch  <- read.table(CurrentABOutGearCatch, fill=TRUE,sep=" ")
    REcosystem_Current_RK4_OutGearCatch <- read.table(CurrentRK4OutGearCatch,fill=TRUE,sep=" ")
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
  setwd(originalWorkingDir)
  forcingOriginalData  <- list(REcosystem_scene$forcing$ForcedBio, REcosystem_scene$forcing$ForcedMigrate)
  typeData             <- list('Forced Bio','Forced Migrate')
  numMonths <- nrow(REcosystem_scene$forcing$ForcedBio)
  REcosystem_scene_jitter <- REcosystem_scene
  for (i in 1:length(forcingOriginalData)) {
    theTypeData  <- typeData[[i]]
    ForcedMatrix <- modifyForcingMatrix(modNum,species,'Jittered',theTypeData,forcingOriginalData[[i]],REcosystem_scene_jitter)
    modNum <- modNum + 1
    if (theTypeData == 'Forced Bio') {
      REcosystem_scene_jitter$forcing$ForcedBio <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedBio_OutBiomass_Jitter,  REcosystem_Baseline_AB_ForcedBio_OutCatch_Jitter,  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Jitter, REcosystem_Baseline_RK4_ForcedBio_OutCatch_Jitter, REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Jitter)
      BaselineJitterFiles  <- list(BaselineABForcedBioOutBiomassJitter,  BaselineABForcedBioOutCatchJitter,  BaselineABForcedBioOutGearCatchJitter,
                                   BaselineRK4ForcedBioOutBiomassJitter, BaselineRK4ForcedBioOutCatchJitter, BaselineRK4ForcedBioOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedBioOutBiomassJitter,  CurrentABForcedBioOutCatchJitter,  CurrentABForcedBioOutGearCatchJitter,
                                   CurrentRK4ForcedBioOutBiomassJitter, CurrentRK4ForcedBioOutCatchJitter, CurrentRK4ForcedBioOutGearCatchJitter)
    } else if (theTypeData == 'Forced Migrate') {
      REcosystem_scene_jitter$forcing$ForcedMigrate <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedMig_OutBiomass_Jitter,  REcosystem_Baseline_AB_ForcedMig_OutCatch_Jitter,  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Jitter, REcosystem_Baseline_RK4_ForcedMig_OutCatch_Jitter, REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Jitter)
      BaselineJitterFiles  <- list(BaselineABForcedMigOutBiomassJitter,  BaselineABForcedMigOutCatchJitter,  BaselineABForcedMigOutGearCatchJitter,
                                   BaselineRK4ForcedMigOutBiomassJitter, BaselineRK4ForcedMigOutCatchJitter, BaselineRK4ForcedMigOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedMigOutBiomassJitter,  CurrentABForcedMigOutCatchJitter,  CurrentABForcedMigOutGearCatchJitter,
                                   CurrentRK4ForcedMigOutBiomassJitter, CurrentRK4ForcedMigOutCatchJitter, CurrentRK4ForcedMigOutGearCatchJitter)
    }
    
    REcosystem_AB_Current_Jitter  <- rsim.run(REcosystem_scene_jitter,method='AB', years=1:50)
    REcosystem_RK4_Current_Jitter <- rsim.run(REcosystem_scene_jitter,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      write.table(REcosystem_AB_Current_Jitter$out_Biomass,     file=BaselineJitterFiles[[1]])
      write.table(REcosystem_AB_Current_Jitter$out_Catch,       file=BaselineJitterFiles[[2]])
      write.table(REcosystem_AB_Current_Jitter$out_Gear_Catch,  file=BaselineJitterFiles[[3]])
      write.table(REcosystem_RK4_Current_Jitter$out_Biomass,    file=BaselineJitterFiles[[4]])
      write.table(REcosystem_RK4_Current_Jitter$out_Catch,      file=BaselineJitterFiles[[5]])
      write.table(REcosystem_RK4_Current_Jitter$out_Gear_Catch, file=BaselineJitterFiles[[6]])
    } else {
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Random", "AB",  "AB",  BaselineJitterTables[[1]], REcosystem_AB_Current_Jitter$out_Biomass,     CurrentJitterFiles[[1]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Random", "AB",  "AB",  BaselineJitterTables[[2]], REcosystem_AB_Current_Jitter$out_Catch,       CurrentJitterFiles[[2]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "AB",  "AB",  BaselineJitterTables[[3]], REcosystem_AB_Current_Jitter$out_Gear_Catch,  CurrentJitterFiles[[3]], species)
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Random", "RK4", "RK4", BaselineJitterTables[[4]], REcosystem_RK4_Current_Jitter$out_Biomass,    CurrentJitterFiles[[4]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Random", "RK4", "RK4", BaselineJitterTables[[5]], REcosystem_RK4_Current_Jitter$out_Catch,      CurrentJitterFiles[[5]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "RK4", "RK4", BaselineJitterTables[[6]], REcosystem_RK4_Current_Jitter$out_Gear_Catch, CurrentJitterFiles[[6]], species)
    }
  }

  print("------------------ Forced Biomass Tests (Stepped) ------------------")
  setwd(originalWorkingDir)
  REcosystem_scene_stepped <- REcosystem_scene
  forcingOriginalData  <- list(REcosystem_scene_stepped$forcing$ForcedBio, REcosystem_scene_stepped$forcing$ForcedMigrate)
  typeData             <- list('Forced Bio','Forced Migrate')
  numMonths <- nrow(REcosystem_scene_stepped$forcing$ForcedBio)
  for (i in 1:length(forcingOriginalData)) {
    theTypeData  <- typeData[[i]]
    ForcedMatrix <- modifyForcingMatrix(modNum,species,'Stepped',theTypeData,forcingOriginalData[[i]],REcosystem_scene_stepped)
    modNum <- modNum + 1
    if (theTypeData == 'Forced Bio') {
      REcosystem_scene_stepped$forcing$ForcedBio <- ForcedMatrix
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedBio_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedBio_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedBio_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedBio_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedBio_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedBio_OutGearCatch_Stepped)
      BaselineSteppedFiles  <- list(BaselineABForcedBioOutBiomassStepped,  BaselineABForcedBioOutCatchStepped,  BaselineABForcedBioOutGearCatchStepped,
                                    BaselineRK4ForcedBioOutBiomassStepped, BaselineRK4ForcedBioOutCatchStepped, BaselineRK4ForcedBioOutGearCatchStepped)
      CurrentSteppedFiles   <- list(CurrentABForcedBioOutBiomassStepped,  CurrentABForcedBioOutCatchStepped,  CurrentABForcedBioOutGearCatchStepped,
                                    CurrentRK4ForcedBioOutBiomassStepped, CurrentRK4ForcedBioOutCatchStepped, CurrentRK4ForcedBioOutGearCatchStepped)
    } else if (theTypeData == 'Forced Migrate') {
      REcosystem_scene_stepped$forcing$ForcedMigrate <- ForcedMatrix
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedMig_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedMig_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedMig_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedMig_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedMig_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedMig_OutGearCatch_Stepped)
      BaselineSteppedFiles  <- list(BaselineABForcedMigOutBiomassStepped,  BaselineABForcedMigOutCatchStepped,  BaselineABForcedMigOutGearCatchStepped,
                                    BaselineRK4ForcedMigOutBiomassStepped, BaselineRK4ForcedMigOutCatchStepped, BaselineRK4ForcedMigOutGearCatchStepped)
      CurrentSteppedFiles   <- list(CurrentABForcedMigOutBiomassStepped,  CurrentABForcedMigOutCatchStepped,  CurrentABForcedMigOutGearCatchStepped,
                                    CurrentRK4ForcedMigOutBiomassStepped, CurrentRK4ForcedMigOutCatchStepped, CurrentRK4ForcedMigOutGearCatchStepped)
    }
    REcosystem_AB_Current_Stepped  <- rsim.run(REcosystem_scene_stepped,method='AB', years=1:50)
    REcosystem_RK4_Current_Stepped <- rsim.run(REcosystem_scene_stepped,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      write.table(REcosystem_AB_Current_Stepped$out_Biomass,     file=BaselineSteppedFiles[[1]])
      write.table(REcosystem_AB_Current_Stepped$out_Catch,       file=BaselineSteppedFiles[[2]])
      write.table(REcosystem_AB_Current_Stepped$out_Gear_Catch,  file=BaselineSteppedFiles[[3]])
      write.table(REcosystem_RK4_Current_Stepped$out_Biomass,    file=BaselineSteppedFiles[[4]])
      write.table(REcosystem_RK4_Current_Stepped$out_Catch,      file=BaselineSteppedFiles[[5]])
      write.table(REcosystem_RK4_Current_Stepped$out_Gear_Catch, file=BaselineSteppedFiles[[6]])
    } else {
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[1]], REcosystem_AB_Current_Stepped$out_Biomass,    CurrentSteppedFiles[[1]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[2]], REcosystem_AB_Current_Stepped$out_Catch,      CurrentSteppedFiles[[2]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "AB", "AB", BaselineSteppedTables[[3]], REcosystem_AB_Current_Stepped$out_Gear_Catch, CurrentSteppedFiles[[3]], species)
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[4]], REcosystem_RK4_Current_Stepped$out_Biomass,   CurrentSteppedFiles[[4]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[5]], REcosystem_RK4_Current_Stepped$out_Catch,     CurrentSteppedFiles[[5]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "RK4","RK4",BaselineSteppedTables[[6]], REcosystem_RK4_Current_Stepped$out_Gear_Catch,CurrentSteppedFiles[[6]], species)
    }
  }
  
  print("------------------ Forced Effort Tests (Jitter) ------------------")
  setwd(originalWorkingDir)
  numMonths <- nrow(REcosystem_scene$fishing$ForcedEffort)
  REcosystem_scene_jitter <- REcosystem_scene
  fishingOriginalData  <- list(REcosystem_scene$fishing$ForcedEffort, REcosystem_scene$fishing$ForcedFRate, REcosystem_scene$fishing$ForcedCatch)
  typeData             <- list('Forced Effort','Forced FRate','Forced Catch')
  for (i in 1:length(fishingOriginalData)) {
    theTypeData  <- typeData[[i]]
    ForcedMatrix <- modifyFishingMatrix(modNum,species,fleets,theTypeData,fishingOriginalData[[i]])
    modNum <- modNum + 1
    if (theTypeData == 'Forced Effort') {
      REcosystem_scene_jitter$fishing$ForcedEffort <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedEff_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedEff_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedEff_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedEffOutBiomassJitter, BaselineABForcedEffOutCatchJitter, BaselineABForcedEffOutGearCatchJitter,
                                   BaselineRK4ForcedEffOutBiomassJitter,BaselineRK4ForcedEffOutCatchJitter,BaselineRK4ForcedEffOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedEffOutBiomassJitter,  CurrentABForcedEffOutCatchJitter,  CurrentABForcedEffOutGearCatchJitter,
                                   CurrentRK4ForcedEffOutBiomassJitter, CurrentRK4ForcedEffOutCatchJitter, CurrentRK4ForcedEffOutGearCatchJitter)
    } else if (theTypeData == 'Forced FRate') {
      REcosystem_scene_jitter$fishing$ForcedFRate  <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedFRa_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Jitter,
                                   REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedFRaOutBiomassJitter, BaselineABForcedFRaOutCatchJitter, BaselineABForcedFRaOutGearCatchJitter,
                                   BaselineRK4ForcedFRaOutBiomassJitter,BaselineRK4ForcedFRaOutCatchJitter,BaselineRK4ForcedFRaOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedFRaOutBiomassJitter,  CurrentABForcedFRaOutCatchJitter,  CurrentABForcedFRaOutGearCatchJitter,
                                   CurrentRK4ForcedFRaOutBiomassJitter, CurrentRK4ForcedFRaOutCatchJitter, CurrentRK4ForcedFRaOutGearCatchJitter)
    } else if (theTypeData == 'Forced Catch') {
      REcosystem_scene_jitter$fishing$ForcedCatch  <- ForcedMatrix
      BaselineJitterTables <- list(REcosystem_Baseline_AB_ForcedCat_OutBiomass_Jitter, REcosystem_Baseline_AB_ForcedCat_OutCatch_Jitter, REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Jitter,
                                REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Jitter,REcosystem_Baseline_RK4_ForcedCat_OutCatch_Jitter,REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Jitter) 
      BaselineJitterFiles  <- list(BaselineABForcedCatOutBiomassJitter, BaselineABForcedCatOutCatchJitter, BaselineABForcedCatOutGearCatchJitter,
                                BaselineRK4ForcedCatOutBiomassJitter,BaselineRK4ForcedCatOutCatchJitter,BaselineRK4ForcedCatOutGearCatchJitter)
      CurrentJitterFiles   <- list(CurrentABForcedCatOutBiomassJitter,  CurrentABForcedCatOutCatchJitter,  CurrentABForcedCatOutGearCatchJitter,
                                CurrentRK4ForcedCatOutBiomassJitter, CurrentRK4ForcedCatOutCatchJitter, CurrentRK4ForcedCatOutGearCatchJitter)
    }
    REcosystem_AB_Current_Jitter  <- rsim.run(REcosystem_scene_jitter,method='AB', years=1:50)
    REcosystem_RK4_Current_Jitter <- rsim.run(REcosystem_scene_jitter,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      write.table(REcosystem_AB_Current_Jitter$out_Biomass,     file=BaselineJitterFiles[[1]])
      write.table(REcosystem_AB_Current_Jitter$out_Catch,       file=BaselineJitterFiles[[2]])
      write.table(REcosystem_AB_Current_Jitter$out_Gear_Catch,  file=BaselineJitterFiles[[3]])
      write.table(REcosystem_RK4_Current_Jitter$out_Biomass,    file=BaselineJitterFiles[[4]])
      write.table(REcosystem_RK4_Current_Jitter$out_Catch,      file=BaselineJitterFiles[[5]])
      write.table(REcosystem_RK4_Current_Jitter$out_Gear_Catch, file=BaselineJitterFiles[[6]])
    } else {
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[1]], REcosystem_AB_Current_Jitter$out_Biomass,    CurrentJitterFiles[[1]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[2]], REcosystem_AB_Current_Jitter$out_Catch,      CurrentJitterFiles[[2]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "AB", "AB",  BaselineJitterTables[[3]], REcosystem_AB_Current_Jitter$out_Gear_Catch, CurrentJitterFiles[[3]], species)
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[4]], REcosystem_RK4_Current_Jitter$out_Biomass,   CurrentJitterFiles[[4]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[5]], REcosystem_RK4_Current_Jitter$out_Catch,     CurrentJitterFiles[[5]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Random", "RK4","RK4", BaselineJitterTables[[6]], REcosystem_RK4_Current_Jitter$out_Gear_Catch,CurrentJitterFiles[[6]], species)
    }
  }

  
  print("------------------ Forced Effort Tests (Stepped) ------------------")
  setwd(originalWorkingDir)
  numMonths <- nrow(REcosystem_scene$fishing$ForcedEffort)
  REcosystem_scene_stepped <- REcosystem_scene
  fishingOriginalData  <- list(REcosystem_scene$fishing$ForcedEffort, REcosystem_scene$fishing$ForcedFRate, REcosystem_scene$fishing$ForcedCatch)
  typeData             <- list('Forced Effort','Forced FRate','Forced Catch')
  for (i in 1:length(fishingOriginalData)) {
    theTypeData  <- typeData[[i]]
    REcosystem_scene_stepped <- REcosystem_scene
    if (theTypeData == 'Forced Effort') {
      REcosystem_scene_stepped$fishing$ForcedEffort <- stepifyMatrix(fishingOriginalData[[i]],fleets)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedEff_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedEff_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedEff_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedEff_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedEff_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedEff_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedEffOutBiomassStepped, BaselineABForcedEffOutCatchStepped, BaselineABForcedEffOutGearCatchStepped,
                                   BaselineRK4ForcedEffOutBiomassStepped,BaselineRK4ForcedEffOutCatchStepped,BaselineRK4ForcedEffOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedEffOutBiomassStepped,  CurrentABForcedEffOutCatchStepped,  CurrentABForcedEffOutGearCatchStepped,
                                   CurrentRK4ForcedEffOutBiomassStepped, CurrentRK4ForcedEffOutCatchStepped, CurrentRK4ForcedEffOutGearCatchStepped)
    } else if (theTypeData == 'Forced FRate') {
      REcosystem_scene_stepped$fishing$ForcedFRate <- stepifyMatrix(fishingOriginalData[[i]],species)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedFRa_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedFRa_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedFRa_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedFRa_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedFRa_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedFRa_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedFRaOutBiomassStepped, BaselineABForcedFRaOutCatchStepped, BaselineABForcedFRaOutGearCatchStepped,
                                   BaselineRK4ForcedFRaOutBiomassStepped,BaselineRK4ForcedFRaOutCatchStepped,BaselineRK4ForcedFRaOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedFRaOutBiomassStepped,  CurrentABForcedFRaOutCatchStepped,  CurrentABForcedFRaOutGearCatchStepped,
                                   CurrentRK4ForcedFRaOutBiomassStepped, CurrentRK4ForcedFRaOutCatchStepped, CurrentRK4ForcedFRaOutGearCatchStepped)      
    } else if (theTypeData == 'Forced Catch') {
      REcosystem_scene_stepped$fishing$ForcedCatch <- stepifyMatrix(fishingOriginalData[[i]],species)
      BaselineSteppedTables <- list(REcosystem_Baseline_AB_ForcedCat_OutBiomass_Stepped,  REcosystem_Baseline_AB_ForcedCat_OutCatch_Stepped,  REcosystem_Baseline_AB_ForcedCat_OutGearCatch_Stepped,
                                    REcosystem_Baseline_RK4_ForcedCat_OutBiomass_Stepped, REcosystem_Baseline_RK4_ForcedCat_OutCatch_Stepped, REcosystem_Baseline_RK4_ForcedCat_OutGearCatch_Stepped)
      BaselineSteppedFiles <- list(BaselineABForcedCatOutBiomassStepped, BaselineABForcedCatOutCatchStepped, BaselineABForcedCatOutGearCatchStepped,
                                   BaselineRK4ForcedCatOutBiomassStepped,BaselineRK4ForcedCatOutCatchStepped,BaselineRK4ForcedCatOutGearCatchStepped)
      CurrentSteppedFiles  <- list(CurrentABForcedCatOutBiomassStepped,  CurrentABForcedCatOutCatchStepped,  CurrentABForcedCatOutGearCatchStepped,
                                   CurrentRK4ForcedCatOutBiomassStepped, CurrentRK4ForcedCatOutCatchStepped, CurrentRK4ForcedCatOutGearCatchStepped)      
    }
    REcosystem_AB_Current_Stepped  <- rsim.run(REcosystem_scene_stepped,method='AB', years=1:50)
    REcosystem_RK4_Current_Stepped <- rsim.run(REcosystem_scene_stepped,method='RK4',years=1:50)
    if (CREATE_BASELINE_FILES) {
      write.table(REcosystem_AB_Current_Stepped$out_Biomass,     file=BaselineSteppedFiles[[1]])
      write.table(REcosystem_AB_Current_Stepped$out_Catch,       file=BaselineSteppedFiles[[2]])
      write.table(REcosystem_AB_Current_Stepped$out_Gear_Catch,  file=BaselineSteppedFiles[[3]])
      write.table(REcosystem_RK4_Current_Stepped$out_Biomass,    file=BaselineSteppedFiles[[4]])
      write.table(REcosystem_RK4_Current_Stepped$out_Catch,      file=BaselineSteppedFiles[[5]])
      write.table(REcosystem_RK4_Current_Stepped$out_Gear_Catch, file=BaselineSteppedFiles[[6]])
    } else {
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[1]], REcosystem_AB_Current_Stepped$out_Biomass,     CurrentSteppedFiles[[1]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[2]], REcosystem_AB_Current_Stepped$out_Catch,       CurrentSteppedFiles[[2]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "AB",  "AB",  BaselineSteppedTables[[3]], REcosystem_AB_Current_Stepped$out_Gear_Catch,  CurrentSteppedFiles[[3]], species)
      runTest(inc(runNum),"out_Biomass",    theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[4]], REcosystem_RK4_Current_Stepped$out_Biomass,    CurrentSteppedFiles[[4]], species)
      runTest(inc(runNum),"out_Catch",      theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[5]], REcosystem_RK4_Current_Stepped$out_Catch,      CurrentSteppedFiles[[5]], species)
      runTest(inc(runNum),"out_Gear_Catch", theTypeData, "Stepped", "RK4", "RK4", BaselineSteppedTables[[6]], REcosystem_RK4_Current_Stepped$out_Gear_Catch, CurrentSteppedFiles[[6]], species)
    }
  }
  
  if (! CREATE_BASELINE_FILES) {
    print(paste0("Completed ",inc(runNum)," test(s)."))
  }
  
})
