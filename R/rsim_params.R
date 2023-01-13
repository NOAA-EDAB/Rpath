#'Initial set up for Rsim module of Rpath
#'
#'Converts the outputs from Rpath into rates for use in Rsim.
#'
#'@family Rsim functions
#'
#'@inheritParams rsim.scenario
#'@param mscramble WILL REMOVE
#'@param mhandle WILL REMOVE
#'@param preyswitch WILL REMOVE - Adjust with adjust.scenario
#'@param scrambleselfwt Value of 1 indicates no overlap while 0 indicates complete overlap.
#'@param handleselfwt Value of 1 indicates no overlap while 0 indicates complete overlap.
#'@param steps_yr Number of time steps per year.
#'@param steps_m Number of time steps per month.
#'
#'@return Returns an \code{Rsim.params} object that is passed to the \code{rsim.run} 
#'    function via the \code{rsim.scenario} function.
#'
#'@export
#'
rsim.params <- function(Rpath, mscramble = 2, mhandle = 1000, preyswitch = 1, 
                        scrambleselfwt = 0, handleselfwt = 0, 
                        steps_yr = 12, steps_m = 1){
  
  nliving <- Rpath$NUM_LIVING
  ndead   <- Rpath$NUM_DEAD
  
  simpar <- c()
  
  simpar$NUM_GROUPS <- Rpath$NUM_GROUPS
  simpar$NUM_LIVING <- nliving
  simpar$NUM_DEAD   <- ndead
  simpar$NUM_GEARS  <- Rpath$NUM_GEARS
  simpar$NUM_BIO    <- simpar$NUM_LIVING + simpar$NUM_DEAD
  simpar$spname     <- c("Outside", Rpath$Group)
  simpar$spnum      <- 0:length(Rpath$Biomass) 
  
  #Energetics for Living and Dead Groups
  #Reference biomass for calculating YY
  simpar$B_BaseRef <- c(1.0, Rpath$Biomass) 
  #Mzero proportional to (1-EE)
  simpar$MzeroMort <- c(0.0, Rpath$PB * (1.0 - Rpath$EE)) 
  #Unassimilated is the proportion of CONSUMPTION that goes to detritus.  
  simpar$UnassimRespFrac <- c(0.0, Rpath$Unassim);
  #Active respiration is proportion of CONSUMPTION that goes to "heat"
  #Passive respiration/ VonB adjustment is left out here
  simpar$ActiveRespFrac <-  c(0.0, ifelse(Rpath$QB > 0, 
                                          1.0 - (Rpath$PB / Rpath$QB) - Rpath$Unassim, 
                                          0.0))
  #Ftime related parameters
  simpar$FtimeAdj   <- rep(0.0, length(simpar$B_BaseRef))
  simpar$FtimeQBOpt <-   c(1.0, ifelse(Rpath$type==1,Rpath$PB,Rpath$QB))
  simpar$PBopt      <-   c(1.0, Rpath$PB)            
  
  #Fishing Effort defaults to 0 for non-gear, 1 for gear
  #KYA EFFORT REMOVED FROM PARAMS July 2015 (was used for gear targeting)
  #simpar$fish_Effort <- ifelse(simpar$spnum <= nliving + ndead,
  #                             0.0,
  #                             1.0) 
  
  #NoIntegrate
  simpar$NoIntegrate <- ifelse(simpar$MzeroMort * simpar$B_BaseRef > 
                                 2 * steps_yr * steps_m, 
                               0, 
                               simpar$spnum)  
  
  #Pred/Prey defaults
  simpar$HandleSelf   <- rep(handleselfwt,   Rpath$NUM_GROUPS + 1)
  simpar$ScrambleSelf <- rep(scrambleselfwt, Rpath$NUM_GROUPS + 1)
  
  #primary production links
  #primTo   <- ifelse(Rpath$PB>0 & Rpath$QB<=0, 1:length(Rpath$PB),0 )
  primTo   <- ifelse(Rpath$type > 0 & Rpath$type <= 1, 
                     1:length(Rpath$PB),
                     0)
  primFrom <- rep(0, length(Rpath$PB))
  primQ    <- Rpath$PB * Rpath$Biomass
  #Change production to consumption for mixotrophs
  mixotrophs <- which(Rpath$type > 0 & Rpath$type < 1)
  primQ[mixotrophs] <- primQ[mixotrophs] / Rpath$GE[mixotrophs] * 
    Rpath$type[mixotrophs] 
  
  #Import links
  # importTo <- ifelse(Rpath$DC[nrow(Rpath$DC), ] > 0,
  #                    1:ncol(Rpath$DC),
  #                    0)
  # importFrom <- rep(0, ncol(Rpath$DC))
  # importQ <- Rpath$DC[nrow(Rpath$DC), ] * Rpath$QB[1:nliving] * Rpath$Biomass[1:nliving]
  
  #Predator/prey links
  preyfrom  <- row(Rpath$DC)
  preyto    <- col(Rpath$DC)	
  predpreyQ <- Rpath$DC[1:(nliving + ndead + 1), ] * 
    t(matrix(Rpath$QB[1:nliving] * Rpath$Biomass[1:nliving],
             nliving, nliving + ndead + 1))
  
  #combined
  simpar$PreyFrom <- c(primFrom[primTo > 0], preyfrom [predpreyQ > 0]) 
  #importFrom[importTo > 0])
  #Changed import prey number to 0
  simpar$PreyFrom[which(simpar$PreyFrom == nrow(Rpath$DC))] <- 0
  simpar$PreyTo   <- c(primTo  [primTo > 0], preyto   [predpreyQ > 0]) 
  #importTo  [importTo > 0])
  simpar$QQ       <- c(primQ   [primTo > 0], predpreyQ[predpreyQ > 0])
  #importQ   [importTo > 0])             	
  
  numpredprey <- length(simpar$QQ)
  
  simpar$DD <- rep(mhandle,   numpredprey)
  simpar$VV <- rep(mscramble, numpredprey)
  
  #NOTE:  Original in C didn't set handleswitch for primary production groups.  Error?
  #probably not when group 0 biomass doesn't change from 1.
  simpar$HandleSwitch <- rep(preyswitch, numpredprey)
  
  #scramble combined prey pools
  Btmp <- simpar$B_BaseRef
  py   <- simpar$PreyFrom + 1.0
  pd   <- simpar$PreyTo + 1.0
  VV   <- simpar$VV * simpar$QQ / Btmp[py]
  AA   <- (2.0 * simpar$QQ * VV) / (VV * Btmp[pd] * Btmp[py] - simpar$QQ * Btmp[pd])
  simpar$PredPredWeight <- AA * Btmp[pd] 
  simpar$PreyPreyWeight <- AA * Btmp[py] 
  
  PredTotWeight <- rep(0, length(simpar$B_BaseRef))
  PreyTotWeight <- rep(0, length(simpar$B_BaseRef))
  
  for(links in 1:numpredprey){
    PredTotWeight[py[links]] <- PredTotWeight[py[links]] + simpar$PredPredWeight[links]
    PreyTotWeight[pd[links]] <- PreyTotWeight[pd[links]] + simpar$PreyPreyWeight[links]    
  }  
  #PredTotWeight[]   <- as.numeric(tapply(simpar$PredPredWeight,py,sum))
  #PreyTotWeight[]   <- as.numeric(tapply(simpar$PreyPreyWeight,pd,sum))
  
  simpar$PredPredWeight <- simpar$PredPredWeight/PredTotWeight[py]
  simpar$PreyPreyWeight <- simpar$PreyPreyWeight/PreyTotWeight[pd]
  
  simpar$NumPredPreyLinks <- numpredprey
  simpar$PreyFrom       <- c(0, simpar$PreyFrom); names(simpar$PreyFrom) <- NULL
  simpar$PreyTo         <- c(0, simpar$PreyTo)  ; names(simpar$PreyTo) <- NULL
  simpar$QQ             <- c(0, simpar$QQ)      ; names(simpar$QQ) <- NULL
  simpar$DD             <- c(mhandle, simpar$DD); names(simpar$DD) <- NULL
  simpar$VV             <- c(mscramble, simpar$VV); names(simpar$VV) <- NULL 
  simpar$HandleSwitch   <- c(0, simpar$HandleSwitch); names(simpar$HandleSwitch) <- NULL 
  simpar$PredPredWeight <- c(0, simpar$PredPredWeight); names(simpar$PredPredWeight) <- NULL
  simpar$PreyPreyWeight <- c(0, simpar$PreyPreyWeight); names(simpar$PreyPreyWeight) <- NULL
  
  #landing links
  fishfrom    <- row(as.matrix(Rpath$Landings))
  fishthrough <- col(as.matrix(Rpath$Landings)) + (nliving + ndead)
  fishcatch   <- Rpath$Landings
  fishto      <- fishfrom * 0
  
  if(sum(fishcatch) > 0){
    simpar$FishFrom    <- fishfrom   [fishcatch > 0]
    simpar$FishThrough <- fishthrough[fishcatch > 0]
    simpar$FishQ       <- fishcatch  [fishcatch > 0] / simpar$B_BaseRef[simpar$FishFrom + 1]  
    simpar$FishTo      <- fishto     [fishcatch > 0]
  }
  #discard links
  
  for(d in 1:Rpath$NUM_DEAD){
    detfate <- Rpath$DetFate[(nliving + ndead + 1):Rpath$NUM_GROUPS, d]
    detmat  <- t(matrix(detfate, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
    
    fishfrom    <-  row(as.matrix(Rpath$Discards))                      
    fishthrough <-  col(as.matrix(Rpath$Discards)) + (nliving + ndead)
    fishto      <-  t(matrix(nliving + d, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
    fishcatch   <-  Rpath$Discards * detmat
    if(sum(fishcatch) > 0){
      simpar$FishFrom    <- c(simpar$FishFrom,    fishfrom   [fishcatch > 0])
      simpar$FishThrough <- c(simpar$FishThrough, fishthrough[fishcatch > 0])
      ffrom <- fishfrom[fishcatch > 0]
      simpar$FishQ       <- c(simpar$FishQ,  fishcatch[fishcatch > 0] / simpar$B_BaseRef[ffrom + 1])  
      simpar$FishTo      <- c(simpar$FishTo, fishto   [fishcatch > 0])
    }
  } 
  
  simpar$NumFishingLinks <- length(simpar$FishFrom)  
  simpar$FishFrom        <- c(0, simpar$FishFrom)    ; names(simpar$FishFrom) <- NULL
  simpar$FishThrough     <- c(0, simpar$FishThrough) ; names(simpar$FishThrough) <- NULL
  simpar$FishQ           <- c(0, simpar$FishQ)       ; names(simpar$FishQ) <- NULL
  simpar$FishTo          <- c(0, simpar$FishTo)      ; names(simpar$FishTo) <- NULL
  
  # SET DETRITAL FLOW
  detfrac <- Rpath$DetFate[1:(nliving + ndead), ]
  detfrom <- row(as.matrix(detfrac))
  detto   <- col(as.matrix(detfrac)) + nliving
  
  detout <- 1 - rowSums(as.matrix(Rpath$DetFate[1:(nliving + ndead), ]))
  dofrom <- 1:length(detout)
  doto   <- rep(0, length(detout))
  
  simpar$DetFrac <- c(0, detfrac[detfrac > 0], detout[detout > 0]); names(simpar$DetFrac) <- NULL
  simpar$DetFrom <- c(0, detfrom[detfrac > 0], dofrom[detout > 0]); names(simpar$DetFrom) <- NULL
  simpar$DetTo   <- c(0, detto  [detfrac > 0], doto  [detout > 0]); names(simpar$DetTo) <- NULL
  simpar$NumDetLinks <- length(simpar$DetFrac) - 1
  
  # STATE VARIABLE DEFAULT 
  #simpar$state_BB    <- simpar$B_BaseRef
  #simpar$state_Ftime <- rep(1, length(Rpath$BB) + 1)
  simpar$BURN_YEARS  <- -1
  simpar$COUPLED     <-  1
  simpar$RK4_STEPS   <- 4.0
  simpar$SENSE_LIMIT <- c(1e-4, 1e4)
  
  # KYA April 2020 Adding names to species-length vectors
  gnames <- as.character(c("Outside",Rpath$Group))
  names(simpar$spname)<-names(simpar$spnum)<-names(simpar$B_BaseRef) <- gnames
  names(simpar$MzeroMort)<-names(simpar$UnassimRespFrac)<-names(simpar$ActiveRespFrac) <- gnames   
  names(simpar$FtimeAdj)<-names(simpar$FtimeQBOpt)<-names(simpar$PBopt) <- gnames            
  names(simpar$NoIntegrate)<-names(simpar$HandleSelf)<-names(simpar$ScrambleSelf) <- gnames
  
  class(simpar) <- "Rsim.params"
  return(simpar)
}