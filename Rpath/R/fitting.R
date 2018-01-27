###########################################################
#'@export
rsim.fit.plot <- function(scene,result,pred){
  OBJ <- rsim.fit.obj(scene,result)
  
  par(mfrow=c(1,2))
  # Plot Catch
    qdat <- OBJ$catch[OBJ$catch$species==pred,]
    mn   <- qdat$obs
    up   <- mn + 1.96*qdat$sd
    dn   <- mn - 1.96*qdat$sd 
    plot(rownames(result$annual_CC),result$annual_CC[,pred],type="l",
       ylim=c(0,max(up,result$annual_CC[,pred])),xlab="",ylab=paste(pred,"catch"))
    points(qdat$year,mn)
    segments(qdat$year,y0=up,y1=dn)
    
  # Plot Biomass
    #qdat <- OBJ$survey[OBJ$survey$species==pred,]
    #mn   <- qdat$obs/qdat$survey_q
    #up   <- mn + 1.96*qdat$sd/qdat$survey_q
    #dn   <- mn - 1.96*qdat$sd/qdat$survey_q 
    plot(rownames(result$annual_BB),result$annual_BB[,pred],type="l",
         ylim=c(0,max(up,result$annual_BB[,pred])),xlab="",ylab=paste(pred,"biomass"))
    #points(qdat$year,mn)
    #segments(qdat$year,y0=up,y1=dn)    
    
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
read.rsim.fitting.catch <- function(insim,flist){

  years <- as.numeric(rownames(insim$fishing$CATCH))
  
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
  #colnames(rsim$fishing$CATCH)<-rsim$params$spname[1:(rsim$params$NUM_BIO+1)]
  #rownames(rsim$fishing$CATCH)<-seq(min(years), length.out=dim(rsim$fishing$CATCH)[1])
  rsim$fishing$CATCH[matrix(c(as.character(rsim$fitting$catch$year),as.character(rsim$fitting$catch$species)),
                           length(rsim$fitting$catch$year),2)] <- rsim$fitting$catch$obs 
  return(rsim)
}

#########################################################
#'@export
read.rsim.fitting.biomass <- function(insim,flist){

  years <- as.numeric(rownames(insim$fishing$CATCH))
  
  rsim <- insim
  if (is.null(rsim$fitting)){rsim$fitting <- list()}  
  
  BIO <- NULL
  for (f in flist){
    cdat  <- read.csv(f)
    ccdat <-cdat[!is.na(cdat$AMOUNT)& cdat$YEAR %in% years ,]        
    wt   <- rep(1,length(ccdat$YEAR))
    type <- as.character(rep("absolute",length(ccdat$YEAR)))
    BIO <- rbind(BIO, data.frame(as.numeric(ccdat$YEAR),as.character(ccdat$SPECIES),
                                     as.numeric(ccdat$AMOUNT/ccdat$SCALE),as.numeric(ccdat$SD),as.numeric(wt),as.character(type),
                                     stringsAsFactors = F))
  }
  
  rsim.fitting.biomass <- BIO
  colnames(rsim$fitting$biomass) <- c("year","species","obs","sd","wt","type")  
    
  return(rsim)
}

#########################################################
#'@export  
read.rsim.fitting.diet <- function(insim,flist){
  
  years <- as.numeric(rownames(insim$fishing$CATCH))
  
  rsim <- insim
  if (is.null(rsim$fitting)){rsim$fitting <- list()}  

  QDAT <- NULL
  for (f in flist){
    rdat <- read.csv(f)
    ddat <- rdat[rdat$year %in% years,]
    QDAT <- rbind(QDAT, data.frame(ddat$year,ddat$pred,ddat$cperwMean,ddat$cperwSD,rep(1,length(ddat[,1]))))
  }

  rsim$fitting$ration <- QDAT
  colnames(rsim$fitting$ration)<-c("year","species","obs","sd","wt")  
  
  rsim$fitting$diet <- read_diet_alphas(insim,flist)
  
  return(rsim)
}
  
#########################################################

read_diet_alphas<-function(insim,flist){
  # Implements log-likelihood from Ainsworth et al. 2010 (Ecol. App. 20: 2188-2202)
  # LogLiklihood = lgamma(alpha0) - sum(lgamma(alpha)) + sum( (alpha-1)*P )
  # lgamma(alpha0) is split over other summed lines to make summing easier
  odat <- NULL
  for (f in flist){  
    ddat <- read.csv(f)
    fpy  <- 1+which(colnames(ddat)=="alpha0")
    lpy  <- length(colnames(ddat))
    pred <- ddat$pred
    prey <- as.character(t(matrix(colnames(ddat)[fpy:lpy],length(colnames(ddat)[fpy:lpy]),length(pred))))
    year <- ddat$year
    almat  <- ddat[,fpy:lpy]
    alpha0 <- rowSums(almat)
    alprop <- as.numeric(unlist(almat/alpha0))
    log_a0part <- alprop*lgamma(alpha0)
    
    alpha     <- as.numeric(unlist(ddat[,fpy:lpy]))
    log_alpha <- lgamma(alpha) 
    log_diff  <- log_a0part-log_alpha
    alphaM1   <- alpha-1
    wt <- rep(1,length(alphaM1))     
    odat<-rbind(odat,data.frame(year,pred,prey,alpha,log_alpha,log_a0part,log_diff,alphaM1,alprop,wt))
  }
  dalph  <- odat[order(odat$pred,odat$year),]
  linkprey <- insim$params$spname[insim$params$PreyFrom+1]
  linkpred <- insim$params$spname[insim$params$PreyTo+1]
  simlink <- rep(NA,length(dalph[,1]))
  for (i in 1:(length(dalph[,1]))){
    simlink[i]<-which(linkprey==dalph$prey[i] & linkpred==dalph$pred[i] )     
  }
  return(data.frame(dalph,simlink))
}  

#########################################################

