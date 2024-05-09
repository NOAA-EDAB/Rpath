source("test-constants.R")

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
stepifyVector <- function(vectorData,type=1, stepFactor) {
  steps     <- getSteps(type)
  numMonths <- length(vectorData)
  incMonths <- numMonths/ NUMBER_OF_STEPS;
  inc <- 0
  for (i in 1:NUMBER_OF_STEPS) {
    start <- inc + 1
    start <- if (start == 0) 1 else start
    end   <- inc+incMonths
    vectorData[start:end] <- vectorData[start:end] * steps[i] * stepFactor
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
stepifyMatrix <- function(forcedEffortMatrix,memberNames,stepFactor) {
  ForcedMatrix <- forcedEffortMatrix
  for (i in 1:length(memberNames)) {
    memberName <- memberNames[i]
    memberCol <- ForcedMatrix[,memberName]
    memberStepped <- stepifyVector(memberCol,i,stepFactor)
    ForcedMatrix[,memberName] <- memberStepped
    # plot(memberStepped,type='l',lwd=5,xlab="Months",ylab="Effort",main=paste0("Forced Effort with Stepped Noise - ",memberName))
  }
  return(ForcedMatrix)
}


#' Stepify Biomass
#'
#' Takes as input a biomass value and stepifies the value by dividing the number of months
#' by the number of steps desired and multiplying each interval by a different step value.
#'
#' @param value : value to stepify across months
#' @param numMonths : number of months in model
#' @param type : type of stepification desired (1, 2, or 3)
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot
#' @param title : title for plot
#'
#' @return Returns the vector of stepified segments
#' 
stepifyBiomass <- function(typeData, value,numMonths,type=1,xlabel,ylabel,title, scaleFactor) {
  steps <- getSteps(type)
  incMonths <- numMonths / NUMBER_OF_STEPS;
  parts <- c()
  if (typeData == FORCED_MIGRATION) {
    value <- scaleFactor*value
  }
  for (i in 1:NUMBER_OF_STEPS) {
    part  <- replicate(incMonths, steps[i]*value)
    parts <- c(parts,part)
  }
  # plot(parts,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)
  return(parts)
}