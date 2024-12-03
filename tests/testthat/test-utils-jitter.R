source("test-constants.R")

#' Add Jitter (i.e., random noise)
#' 
#' Adds random noise to the specified matrix.
#'
#' @param matrix : data matrix to be jittered
#' @param factor : the R jitter factor (the larger the number the greater the jitter amount)
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot 
#' @param title : main title for plot
#'
#' @return Returns the jittered matrix
#' 
addJitter <- function(matrix,xlabel,ylabel,title) {
  # From jitter() doc: If amount == 0, jitter returns factor * z/50, where
  # z = max(x0) - min(x), aka the range. So if factor=5 and amount=0, jitter()
  # returns a random value within a tenth of the range.
  
  jitteredMatrix <- jitter(matrix,factor=FACTOR_VALUE,amount=NULL)
  if (xlabel != '' && ylabel != '' & title != '') {
    # plot(jitteredMatrix,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)    
  } 
  return(jitteredMatrix)
}

#' Create a jittered vector
#' 
#' Used for jittering a matrix column.
#'
#' @param typeData : the type of data to be jittered
#' @param value : value to add jitter to
#' @param numElements : number of elements in vector
#' @param xlabel : x axis label for plot
#' @param ylabel : y axis label for plot
#' @param title : main title for plot
#' @param title : the type of random number (positive only or any)
#'
#' @return Returns the column-jittered matrix
#' 
createJitterVectorFromValue <- function(typeData,value,numElements,xlabel,ylabel,title,randomNumberType) {
  jitterVector <- c()
  migrationScaleFactor <- 1
  if (typeData == FORCED_MIGRATION) {
      migrationScaleFactor = FORCED_MIGRATION_SCALE_FACTOR_JITTER # = 1000
  }
  for (i in 1:numElements) {
    randVal <- randomNumber(migrationScaleFactor*JITTER_AMOUNT_PCT,randomNumberType)
    jitteredValue <- value * (1.0 + randVal)
    jitterVector <- append(jitterVector,jitteredValue)
  }
  # plot(jitterVector,type='l',lwd=5,xlab=xlabel,ylab=ylabel,main=title)
  return(jitterVector)
}