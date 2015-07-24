#'@export
ecosim.random.params <- function(Rpath){
  #
  # TODO Add Randomness to simpar below based on matrix Rpath$pedigree
  # ## MZERO USED AS AN EXAMPLE BELOW
  
  #Old path_to_rates--------------------------------------------------------------------
  MSCRAMBLE      <- 2.0
  MHANDLE        <- 1000.0
  PREYSWITCH     <- 1.0
  # For SelfWts 1.0 = no overlap, 0.0 = complete overlap
  ScrambleSelfWt <- 1.0
  HandleSelfWt   <- 1.0
  
  simpar <- c()
  
  simpar$NUM_GROUPS <- Rpath$NUM_GROUPS
  simpar$NUM_LIVING <- Rpath$NUM_LIVING
  simpar$NUM_DEAD   <- Rpath$NUM_DEAD
  simpar$NUM_GEARS  <- Rpath$NUM_GEARS
  
  simpar$spname     <- c("Outside", Rpath$Group)
  simpar$spnum      <- 0:length(Rpath$BB) 
  
  #Energetics for Living and Dead Groups
  #Reference biomass for calculating YY
  simpar$B_BaseRef <- c(1.0, Rpath$BB) 
  #Mzero proportional to (1-EE)
  # SCALE TO Mzero +/- Mzero * PBpedigree *uniform(-1,1) 
  simpar$MzeroMort <- c(0.0, Rpath$PB * (1.0 - Rpath$EE))  *
        (1 + EBS$pedigree[,2] * runif(EBS$NUM_LIVING+EBS$NUM_DEAD,-1.0,1.0) )
            
  #Unassimilated is the proportion of CONSUMPTION that goes to detritus.  
  simpar$UnassimRespFrac <- c(0.0, Rpath$GS);
  #Active respiration is proportion of CONSUMPTION that goes to "heat"
  #Passive respiration/ VonB adjustment is left out here
  simpar$ActiveRespFrac <-  c(0.0, ifelse(Rpath$QB > 0, 
                                          1.0 - (Rpath$PB / Rpath$QB) - Rpath$GS, 
                                          0.0))
  #Ftime related parameters
  simpar$FtimeAdj   <- rep(0.0, length(simpar$B_BaseRef))
  simpar$FtimeQBOpt <-   c(1.0, Rpath$QB)
  simpar$PBopt      <-   c(1.0, Rpath$PB)           
  
  #Fishing Effort defaults to 0 for non-gear, 1 for gear
  #KYA EFFORT REMOVED FROM PARAMS July 2015
  simpar$fish_Effort <- ifelse(simpar$spnum <= simpar$NUM_LIVING + simpar$NUM_DEAD,
                               0.0,
                               1.0) 
  
  #NoIntegrate
  STEPS_PER_YEAR  <- 12
  STEPS_PER_MONTH <- 1
  simpar$NoIntegrate <- ifelse(c(0, Rpath$PB) / 
                                 (1.0 - simpar$ActiveRespFrac - simpar$UnassimRespFrac) > 
                                 2 * STEPS_PER_YEAR * STEPS_PER_MONTH, 
                               0, 
                               simpar$spnum)  
  
  #Pred/Prey defaults
  simpar$HandleSelf   <- rep(HandleSelfWt,   Rpath$NUM_GROUPS + 1)
  simpar$ScrambleSelf <- rep(ScrambleSelfWt, Rpath$NUM_GROUPS + 1)
  
  #primary production links
  #primTo   <- ifelse(Rpath$PB>0 & Rpath$QB<=0, 1:length(Rpath$PB),0 )
  primTo   <- ifelse(Rpath$type == 1, 
                     1:length(Rpath$PB),
                     0)
  primFrom <- rep(0, length(Rpath$PB))
  primQ    <- Rpath$PB * Rpath$BB 
  
  #Predator/prey links
  preyfrom  <- row(Rpath$DC)
  preyto    <- col(Rpath$DC)	
  predpreyQ <- Rpath$DC * t(matrix(Rpath$QB[1:Rpath$NUM_LIVING] * Rpath$BB[1:Rpath$NUM_LIVING],
                                   Rpath$NUM_LIVING, Rpath$NUM_LIVING + Rpath$NUM_DEAD))
  
  #combined
  simpar$PreyFrom <- c(primFrom[primTo > 0], preyfrom [predpreyQ > 0])
  simpar$PreyTo   <- c(primTo  [primTo > 0], preyto   [predpreyQ > 0])
  simpar$QQ       <- c(primQ   [primTo > 0], predpreyQ[predpreyQ > 0])             	
  
  numpredprey <- length(simpar$QQ)
  
  simpar$DD <- rep(MHANDLE,   numpredprey)
  simpar$VV <- rep(MSCRAMBLE, numpredprey)
  
  #NOTE:  Original in C didn't set handleswitch for primary production groups.  Error?
  #probably not when group 0 biomass doesn't change from 1.
  simpar$HandleSwitch <- rep(PREYSWITCH, numpredprey)
  
  #scramble combined prey pools
  Btmp <- simpar$B_BaseRef
  py   <- simpar$PreyFrom + 1.0
  pd   <- simpar$PreyTo + 1.0
  VV   <- simpar$VV * simpar$QQ / Btmp[py]
  AA   <- (2.0 * simpar$QQ * VV) / (VV * Btmp[pd] * Btmp[py] - simpar$QQ * Btmp[pd])
  simpar$PredPredWeight <- AA * Btmp[pd] 
  simpar$PreyPreyWeight <- AA * Btmp[py] 
  
  simpar$PredTotWeight <- rep(0, length(simpar$B_BaseRef))
  simpar$PreyTotWeight <- rep(0, length(simpar$B_BaseRef))
  
  for(links in 1:numpredprey){
    simpar$PredTotWeight[py[links]] <- simpar$PredTotWeight[py[links]] + simpar$PredPredWeight[links]
    simpar$PreyTotWeight[pd[links]] <- simpar$PreyTotWeight[pd[links]] + simpar$PreyPreyWeight[links]    
  }  
  #simpar$PredTotWeight[]   <- as.numeric(tapply(simpar$PredPredWeight,py,sum))
  #simpar$PreyTotWeight[]   <- as.numeric(tapply(simpar$PreyPreyWeight,pd,sum))
  
  simpar$PredPredWeight <- simpar$PredPredWeight/simpar$PredTotWeight[py]
  simpar$PreyPreyWeight <- simpar$PreyPreyWeight/simpar$PreyTotWeight[pd]
  
  simpar$NumPredPreyLinks <- numpredprey
  simpar$PreyFrom       <- c(0, simpar$PreyFrom)
  simpar$PreyTo         <- c(0, simpar$PreyTo)
  simpar$QQ             <- c(0, simpar$QQ)
  simpar$DD             <- c(0, simpar$DD)
  simpar$VV             <- c(0, simpar$VV) 
  simpar$HandleSwitch   <- c(0, simpar$HandleSwitch) 
  simpar$PredPredWeight <- c(0, simpar$PredPredWeight)
  simpar$PreyPreyWeight <- c(0, simpar$PreyPreyWeight)
  
  #catchlinks
  fishfrom    <- row(as.matrix(Rpath$Catch))
  fishthrough <- col(as.matrix(Rpath$Catch)) + (Rpath$NUM_LIVING + Rpath$NUM_DEAD)
  fishcatch   <- Rpath$Catch
  fishto      <- fishfrom * 0
  
  if(sum(fishcatch) > 0){
    simpar$FishFrom    <- fishfrom   [fishcatch > 0]
    simpar$FishThrough <- fishthrough[fishcatch > 0]
    simpar$FishQ       <- fishcatch  [fishcatch > 0] / simpar$B_BaseRef[simpar$FishFrom + 1]  
    simpar$FishTo      <- fishto     [fishcatch > 0]
  }
  #discard links
  
  for(d in 1:Rpath$NUM_DEAD){
    detfate <- Rpath$DetFate[(Rpath$NUM_LIVING + Rpath$NUM_DEAD + 1):Rpath$NUM_GROUPS, d]
    detmat  <- t(matrix(detfate, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
    
    fishfrom    <-  row(as.matrix(Rpath$Discards))                      
    fishthrough <-  col(as.matrix(Rpath$Discards)) + (Rpath$NUM_LIVING + Rpath$NUM_DEAD)
    fishto      <-  t(matrix(Rpath$NUM_LIVING + d, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
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
  simpar$FishFrom        <- c(0, simpar$FishFrom)
  simpar$FishThrough     <- c(0, simpar$FishThrough)
  simpar$FishQ           <- c(0, simpar$FishQ)  
  simpar$FishTo          <- c(0, simpar$FishTo)   
  
# SET DETRITAL FLOW
  detfrac <- Rpath$DetFate[1:(Rpath$NUM_LIVING + Rpath$NUM_DEAD), ]
  detfrom <- row(as.matrix(detfrac))
  detto   <- col(as.matrix(detfrac)) + Rpath$NUM_LIVING
  
  detout <- 1 - rowSums(as.matrix(Rpath$DetFate[1:(Rpath$NUM_LIVING + Rpath$NUM_DEAD), ]))
  dofrom <- 1:length(detout)
  doto   <- rep(0, length(detout))
  
  simpar$DetFrac <- c(0, detfrac[detfrac > 0], detout[detout > 0])
  simpar$DetFrom <- c(0, detfrom[detfrac > 0], dofrom[detout > 0])
  simpar$DetTo   <- c(0, detto  [detfrac > 0], doto  [detout > 0])
  simpar$NumDetLinks <- length(simpar$DetFrac) - 1
  
# STATE VARIABLE DEFAULT 
  #simpar$state_BB    <- simpar$B_BaseRef
  #simpar$state_Ftime <- rep(1, length(Rpath$BB) + 1)
  simpar$BURN_YEARS <- -1
  simpar$COUPLED    <-  1
  simpar$CRASH_YEAR <- -1
  
  return(simpar)
}

