
#'@export

read.rsim.fitting.catch <- function(insim,years,flist){
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
  
  # APPLY FISHING TO FITTING
  rsim$fishing$EFFORT[]<-0
  colnames(rsim$fishing$CATCH)<-rsim$params$spname[1:(rsim$params$NUM_BIO+1)]
  rownames(rsim$fishing$CATCH)<-seq(min(years), length.out=dim(rsim$fishing$CATCH)[1])
  rsim$fishing$CATCH[matrix(c(as.character(rsim$fitting$catch$year),as.character(rsim$fitting$catch$species)),
                           length(rsim$fitting$catch$year),2)] <- rsim$fitting$catch$obs 
  
  
  return(rsim)
  
}