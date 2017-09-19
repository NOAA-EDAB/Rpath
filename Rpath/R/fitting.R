
#'@export

read.rsim.fitting.catch <- function(insim,flist){
  rsim <- insim
  
  if (is.null(rsim$fitting)){rsim$fitting <- list()}  
  
  CATCH <- NULL
  
  for (f in flist){
      cdat  <- read.csv(f)
      ccdat <- cdat[!is.na(cdat$AMOUNT) & cdat$YEAR %in% years,]
      sdat  <- aggregate(as.numeric(ccdat$AMOUNT)/as.numeric(ccdat$SCALE),list(ccdat$YEAR,ccdat$SPECIES),"sum")
      sd    <- 0.1*sdat$x
      wt   <- rep(1,length(sdat$x))
      CATCH <- rbind(CATCH, cbind(sdat,sd,wt))
  }
  rsim$fitting$catch <- CATCH
  colnames(rsim$fitting$catch) <- c("year","species","obs","sd","wt")
  return(rsim)
  
}