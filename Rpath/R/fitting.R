###########################################################
#'@export
rsim.fit.plot.catch <- function(scene,result,pred){
  OBJ <- rsim.fit.obj(scene,result)
  qdat <- OBJ$catch[OBJ$catch$species==pred,]
  mn   <- qdat$obs
  up   <- mn + 1.96*qdat$sd
  dn   <- mn - 1.96*qdat$sd 
  plot(rownames(result$annual_CC),result$annual_CC[,pred],type="l",
       ylim=c(0,max(up,result$annual_CC[,pred])),xlab="",ylab=paste(pred,"catch"))
  points(qdat$year,mn)
  segments(qdat$year,y0=up,y1=dn)  
}

###########################################################
#'@export
rsim.fit.obj <- function(scene,result){

  ONEHALFLOGTWOPI <- 0.5*log(2*pi) #0.918938533204672
  epsilon <- 1e-36
  
  OBJ <- list() 
  
# Catch compared (assumes all catch is clean, absolute values)
  est <- result$annual_CC[matrix(c(as.character(scene$fitting$catch$year),as.character(scene$fitting$catch$species)),
                              ncol=2)] + epsilon
  obs <- scene$fitting$catch$obs + epsilon
  sd  <- scene$fitting$catch$sd  + epsilon
  sdlog  <- sqrt(log(1.0+sd*sd/(obs*obs)))
  sdiff  <- (log(obs)-log(est))/sdlog
  fit    <- scene$fitting$catch$wt * (log(sdlog) + ONEHALFLOGTWOPI + 0.5*sdiff*sdiff)
  OBJ$catch <- cbind(scene$fitting$catch,est,sdiff,fit)  

# Final summation and return
  OBJ$tot <- sum(OBJ$catch$fit) #OBJ$survey$fit, OBJ$diet$fit

  return(OBJ)    
}

#########################################################
#'@export
read.rsim.fitting.catch <- function(insim,years,flist){
  rsim <- insim
  
  #if (is.null(rsim$fitting)){
     rsim$fitting <- list(years=years, NY=length(years) )
  #}  
  
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
  colnames(rsim$fishing$CATCH)<-rsim$params$spname[1:(rsim$params$NUM_BIO+1)]
  rownames(rsim$fishing$CATCH)<-seq(min(years), length.out=dim(rsim$fishing$CATCH)[1])
  rsim$fishing$CATCH[matrix(c(as.character(rsim$fitting$catch$year),as.character(rsim$fitting$catch$species)),
                           length(rsim$fitting$catch$year),2)] <- rsim$fitting$catch$obs 
  return(rsim)
}