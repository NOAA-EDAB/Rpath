################################################################################
# R version of ecosim (to work with ecosim.dll c code) by Kerim Aydin
#
# Version 0.04
################################################################################ 

##------------------------------------------------------------------------------ 


load_ecosim   <- function(Years, Pathfile, Dietfile, Pedigreefile, Juvfile){

    # Load and balance Ecopath	
	  testpath  <- ecopath(Pathfile, Dietfile, Pedigreefile, FALSE)
    # Convert Ecopath to ecosim rates
      testrates <- path_to_rates(testpath)
    # Add juveniles
      testjuvs  <- initialize_stanzas(testrates, Juvfile)  
    # Put in arrays appropriate for c
      simout    <- ecosim_pack(testjuvs, Years)
  
   return(simout) 
}  

##------------------------------------------------------------------------------ 
# Load and unload dll (unloading needed before compiling any c-code changes)

# Load
  ll<-function(){
    if (!is.loaded("ecosim_run")){dyn.load("ecosim.dll")}
  } 

# Reload
  rl<-function(){
	if (is.loaded("ecosim_run")){dyn.unload("ecosim.dll")} 
	dyn.load("ecosim.dll")		
  }

# Unload
  ul<-function(){
	if (is.loaded("ecosim_run")){dyn.unload("ecosim.dll")} 	
  }
##------------------------------------------------------------------------------    

path_to_rates<-function(path){

  MSCRAMBLE      <- 2.0;
  MHANDLE        <- 1000.0;
  PREYSWITCH     <- 1.0;
  # For SelfWts 1.0 = no overlap, 0.0 = complete overlap
  ScrambleSelfWt <- 1.0;
  HandleSelfWt   <- 1.0;

    rpar <- NULL
    rpar$NUM_GROUPS <- path$NUM_GROUPS
    rpar$NUM_LIVING <- path$NUM_LIVING
    rpar$NUM_DEAD   <- path$NUM_DEAD
    rpar$NUM_GEARS  <- path$NUM_GEARS
    
    rpar$spname          <- c("Outside",path$Group)
    rpar$spnum           <- 0:length(path$BB) 
# Energetics for Living and Dead Groups
  # Reference biomass for calculating YY
    rpar$B_BaseRef       <- c(1.0,path$BB); 
  # Mzero proportional to (1-EE)
    rpar$MzeroMort       <- c(0.0, path$PB * (1.0 - path$EE)); 
  # Unassimilated is the proportion of CONSUMPTION that goes to detritus.  
    rpar$UnassimRespFrac <- c(0.0, path$GS);
  # Active respiration is proportion of CONSUMPTION that goes to "heat"
  # Passive respiration/ VonB adjustment is left out here
    rpar$ActiveRespFrac <-  c(0.0,
         ifelse(path$QB>0, 1.0 - (path$PB/path$QB) - path$GS, 0.0))
  # Ftime related parameters
    rpar$FtimeAdj   <- rep(0.0, length(rpar$B_BaseRef))
    rpar$FtimeQBOpt <- c(1.0, path$QB)
    rpar$PBopt      <- c(1.0, path$PB)           

  # Fishing Effort defaults to 0 for non-gear, 1 for gear
    rpar$fish_Effort <- ifelse(rpar$spnum<=rpar$NUM_LIVING+rpar$NUM_DEAD,0.0,1.0) 
    
  # NoIntegrate
    STEPS_PER_YEAR <- 12; STEPS_PER_MONTH <- 1;
    rpar$NoIntegrate <- ifelse( c(0,path$PB) /
                 (1.0 - rpar$ActiveRespFrac - rpar$UnassimRespFrac) >
 				 2 * STEPS_PER_YEAR * STEPS_PER_MONTH, 0, rpar$spnum)  

  # Pred/Prey defaults
    rpar$HandleSelf   <- rep(HandleSelfWt,   path$NUM_GROUPS+1)
    rpar$ScrambleSelf <- rep(ScrambleSelfWt, path$NUM_GROUPS+1)

# primary production links
  #primTo   <- ifelse(path$PB>0 & path$QB<=0, 1:length(path$PB),0 )
  primTo   <- ifelse(path$type == 1, 1:length(path$PB),0 )
  primFrom <- rep(0,length(path$PB))
  primQ    <- path$PB * path$BB 
  
# Predator/prey links
  preyfrom  <- row(path$DC)
  preyto    <- col(path$DC)	
  predpreyQ <- path$DC * t(matrix(path$QB[1:path$NUM_LIVING] * path$BB[1:path$NUM_LIVING],
                           path$NUM_LIVING,path$NUM_LIVING+path$NUM_DEAD))

# combined
  rpar$PreyFrom <- c(primFrom[primTo>0] , preyfrom[predpreyQ>0])
  rpar$PreyTo   <- c(primTo[primTo>0]   , preyto[predpreyQ>0])
  rpar$QQ       <- c(primQ[primTo>0]    , predpreyQ[predpreyQ>0])             	
  
  numpredprey <- length(rpar$QQ)

  rpar$DD           <- rep(MHANDLE,    numpredprey)
  rpar$VV           <- rep(MSCRAMBLE,  numpredprey)
# NOTE:  Original in C didn't set handleswitch for primary production groups.  Error?
# probably not when group 0 biomass doesn't change from 1.
  rpar$HandleSwitch <- rep(PREYSWITCH,numpredprey);

# scramble combined prey pools
  Btmp <- rpar$B_BaseRef
  py   <- rpar$PreyFrom + 1.0
  pd   <- rpar$PreyTo + 1.0
  VV <- rpar$VV * rpar$QQ / Btmp[py]
  AA <- (2.0 * rpar$QQ * VV) / (VV * Btmp[pd] * Btmp[py] - rpar$QQ * Btmp[pd])
  rpar$PredPredWeight  <- AA * Btmp[pd]; 
  rpar$PreyPreyWeight  <- AA * Btmp[py]; 

  rpar$PredTotWeight   <- rep(0,length(rpar$B_BaseRef))
  rpar$PreyTotWeight   <- rep(0,length(rpar$B_BaseRef))
  for (links in 1:numpredprey){
    rpar$PredTotWeight[py[links]] = rpar$PredTotWeight[py[links]] + rpar$PredPredWeight[links]
    rpar$PreyTotWeight[pd[links]] = rpar$PreyTotWeight[pd[links]] + rpar$PreyPreyWeight[links]    
  }  
  #rpar$PredTotWeight[]   <- as.numeric(tapply(rpar$PredPredWeight,py,sum))
  #rpar$PreyTotWeight[]   <- as.numeric(tapply(rpar$PreyPreyWeight,pd,sum))

  rpar$PredPredWeight = rpar$PredPredWeight/rpar$PredTotWeight[py]
  rpar$PreyPreyWeight = rpar$PreyPreyWeight/rpar$PreyTotWeight[pd]

  rpar$NumPredPreyLinks <- numpredprey
  rpar$PreyFrom       <- c(0,rpar$PreyFrom)
  rpar$PreyTo         <- c(0,rpar$PreyTo)
  rpar$QQ             <- c(0,rpar$QQ)
  rpar$DD             <- c(0,rpar$DD)
  rpar$VV             <- c(0,rpar$VV) 
  rpar$HandleSwitch   <- c(0,rpar$HandleSwitch) 
  rpar$PredPredWeight <- c(0,rpar$PredPredWeight)
  rpar$PreyPreyWeight <- c(0,rpar$PreyPreyWeight)
  
# catchlinks
  fishfrom    <-  row(path$Catch)                      
  fishthrough <-  col(path$Catch) + (path$NUM_LIVING + path$NUM_DEAD)
  fishcatch   <-  path$Catch
  fishto      <-  fishfrom * 0
  
  if (sum(fishcatch)>0){
      rpar$FishFrom    <- fishfrom[fishcatch>0]
      rpar$FishThrough <- fishthrough[fishcatch>0]
      rpar$FishQ       <- fishcatch[fishcatch>0] / rpar$B_BaseRef[rpar$FishFrom+1]  
      rpar$FishTo      <- fishto[fishcatch>0]
  }
# discard links
  
  for (d in 1:path$NUM_DEAD){
   detfate <- path$DetFate[(path$NUM_LIVING + path$NUM_DEAD+1):path$NUM_GROUPS,d]
   detmat  <- t(matrix(detfate,path$NUM_GEAR,path$NUM_GROUPS))
   
   fishfrom    <-  row(path$Discards)                      
   fishthrough <-  col(path$Discards) + (path$NUM_LIVING + path$NUM_DEAD)
   fishto      <-  t(matrix(path$NUM_LIVING + d,path$NUM_GEAR,path$NUM_GROUPS))
   fishcatch   <-  path$Discards * detmat
   if (sum(fishcatch)>0){
      rpar$FishFrom    <- c(rpar$FishFrom,  fishfrom[fishcatch>0])
      rpar$FishThrough <- c(rpar$FishThrough, fishthrough[fishcatch>0])
      ffrom = fishfrom[fishcatch>0]
      rpar$FishQ       <- c(rpar$FishQ, fishcatch[fishcatch>0] / rpar$B_BaseRef[ffrom+1])  
      rpar$FishTo      <- c(rpar$FishTo, fishto[fishcatch>0])
   }
  } 
  
  rpar$NumFishingLinks <- length(rpar$FishFrom)  
  rpar$FishFrom    <- c(0, rpar$FishFrom)
  rpar$FishThrough <- c(0, rpar$FishThrough)
  rpar$FishQ       <- c(0, rpar$FishQ)  
  rpar$FishTo      <- c(0, rpar$FishTo)   

# Unwound discard fate for living groups
  detfrac <- path$DetFate[1:(path$NUM_LIVING + path$NUM_DEAD),]
  detfrom <- row(detfrac)
  detto   <- col(detfrac) + path$NUM_LIVING
  
  detout   <- 1-rowSums(path$DetFate[1:(path$NUM_LIVING + path$NUM_DEAD),])
  dofrom   <- 1:length(detout)
  doto     <- rep(0,length(detout))
  
  rpar$DetFrac <- c(0,detfrac[detfrac>0],detout[detout>0])
  rpar$DetFrom <- c(0,detfrom[detfrac>0],dofrom[detout>0])
  rpar$DetTo   <- c(0,detto[detfrac>0],doto[detout>0])
  
  rpar$NumDetLinks <- length(rpar$DetFrac) - 1

  rpar$state_BB    <- rpar$B_BaseRef
  rpar$state_Ftime <- rep(1,length(path$BB)+1)
  

return(rpar)

}


###############################################################################

# Initialize juvenile adult or "stanza" age structure in sim
initialize_stanzas<-function(rpar,juv_file){    

  STEPS_PER_YEAR     <- 12
  MAX_MONTHS_STANZA  <- 400
  MIN_REC_FACTOR     <- 6.906754779  #// This is ln(1/0.001 - 1) used to set min. logistic matuirty to 0.001

  jpar       <- rpar
  juv        <- read.csv(juv_file)
  jpar$juv_N <- length (juv$JuvNum) 
    
  # fill monthly vectors for each species, rel weight and consumption at age
  # Loops for weight and numbers
    jpar$WageS      <- matrix(0,MAX_MONTHS_STANZA+1,jpar$juv_N+1)
    jpar$WWa        <- matrix(0,MAX_MONTHS_STANZA+1,jpar$juv_N+1)
    jpar$NageS      <- matrix(0,MAX_MONTHS_STANZA+1,jpar$juv_N+1)
    jpar$SplitAlpha <- matrix(0,MAX_MONTHS_STANZA+1,jpar$juv_N+1)

    jpar$state_NN       <-rep(0.0,rpar$NUM_GROUPS+1)
    jpar$stanzaPred     <-rep(0.0,rpar$NUM_GROUPS+1)
    jpar$stanzaBasePred <-rep(0.0,rpar$NUM_GROUPS+1)

    jpar$SpawnBio       <- rep(0.0,jpar$juv_N+1)
    jpar$EggsStanza     <- rep(0.0,jpar$juv_N+1)
    jpar$stanzaGGJuv    <- rep(0.0,jpar$juv_N+1)
    jpar$stanzaGGAdu    <- rep(0.0,jpar$juv_N+1)
    #jpar$drawn_K[MAX_SPLIT+1]
    #jpar$drawn_AduZ[MAX_SPLIT+1]
    #jpar$drawn_JuvZ[MAX_SPLIT+1]
    jpar$SpawnEnergy <- rep(0.0,jpar$juv_N+1)
    jpar$SpawnX      <- rep(0.0,jpar$juv_N+1)
    jpar$SpawnAllocR <- rep(0.0,jpar$juv_N+1)
    jpar$SpawnAllocG <- rep(0.0,jpar$juv_N+1)

    jpar$recruits       <- rep(0.0,jpar$juv_N+1)
    jpar$RzeroS         <- rep(0.0,jpar$juv_N+1)
    jpar$baseEggsStanza <- rep(0.0,jpar$juv_N+1)

    jpar$baseSpawnBio   <- rep(0.0,jpar$juv_N+1)
    jpar$Rbase          <- rep(0.0,jpar$juv_N+1)
    jpar$RscaleSplit    <- rep(0.0,jpar$juv_N+1)
    
    jpar$JuvNum      <- c(0,juv$JuvNum)
    jpar$AduNum      <- c(0,juv$AduNum)
    jpar$RecMonth    <- c(0,juv$RecMonth)
    jpar$VonBD       <- c(0,juv$VonBD)
	jpar$aduEqAgeZ   <- c(0,juv$JuvZ_BAB)
	jpar$juvEqAgeZ   <- c(0,juv$AduZ_BAB)
    jpar$RecPower    <- c(0,juv$RecPower)
    jpar$Wmat001     <- c(0,juv$Wmat001)
    jpar$Wmat50      <- c(0,juv$Wmat50)
    jpar$Amat001     <- c(0,juv$Amat001)
    jpar$Amat50      <- c(0,juv$Amat50)
    
  # first calculate the number of months in each age group.  
  # For end stanza (plus group), capped in EwE @400, this is 90% of 
  # relative max wt, uses generalized exponent D
    firstMoJuv <- rep(0,jpar$juv_N) 
    lastMoJuv  <- juv$RecAge * STEPS_PER_YEAR - 1;
    firstMoAdu <- juv$RecAge * STEPS_PER_YEAR;
    lastMoAdu  <- floor( log (1.0 - (0.9 ^ (1.0 - juv$VonBD))) /
                              ( -3.0 * juv$VonBK * (1.0 - juv$VonBD) / STEPS_PER_YEAR))
    lastMoAdu <- ifelse(lastMoAdu>MAX_MONTHS_STANZA, MAX_MONTHS_STANZA, lastMoAdu)

  # Energy required to grow a unit in weight (scaled to Winf = 1)
    vBM = 1.0 - 3.0 * juv$VonBK/STEPS_PER_YEAR;

  # Set Q multiplier to 1 (used to rescale consumption rates) 
    Qmult <- rep(1.0,jpar$NUM_GROUPS) 
 
  # Survival rates for juveniles and adults
  # KYA Switched survival rate out of BAB == assume all on adults
    survRate_juv_month = exp(-(juv$JuvZ_BAB)/ STEPS_PER_YEAR);
    survRate_adu_month = exp(-(juv$AduZ_BAB)/ STEPS_PER_YEAR);

    WmatSpread = -(juv$Wmat001 - juv$Wmat50)/MIN_REC_FACTOR
    AmatSpread = -(juv$Amat001 - juv$Amat50)/MIN_REC_FACTOR

  # Save numbers that can be saved at this point
    jpar$WmatSpread  <- c(0,WmatSpread)    
    jpar$AmatSpread  <- c(0,AmatSpread)
    jpar$vBM         <- c(0,vBM)
    jpar$firstMoJuv  <- c(0, firstMoJuv)
    jpar$lastMoJuv   <- c(0, lastMoJuv)
    jpar$firstMoAdu  <- c(0, firstMoAdu)
    jpar$lastMoAdu   <- c(0, lastMoAdu)       
    
	if (jpar$juv_N>0){                        
    for (i in 1:jpar$juv_N){
	  # R rpar lookup number for juvs and adults is shifted by 1
        rJuvNum <- juv$JuvNum[i] + 1
        rAduNum <- juv$AduNum[i] + 1
                             
      # vector of ages in months from month 0
        ageMo <- 0:lastMoAdu[i]
        
      # Weight and consumption at age
        WageS <- (1.0 - exp(-3. * juv$VonBK[i] * (1. - juv$VonBD[i]) * ageMo/STEPS_PER_YEAR)) ^ (1.0/(1.0 - juv$VonBD[i]))
        WWa <- WageS ^ juv$VonBD[i]
        
      # Numbers per recruit
        surv <- ifelse( ageMo<=firstMoAdu[i], survRate_juv_month[i]^ageMo,
                        survRate_juv_month[i]^firstMoAdu[i] * 
                        survRate_adu_month[i]^(ageMo-firstMoAdu[i])
                      )
      # plus group
        surv[lastMoAdu[i]+1] <- surv[lastMoAdu[i]] * survRate_adu_month[i] / (1.0 - survRate_adu_month[i])
            
      # correct if using monthly pulse recruitment
        if (juv$RecMonth[i]>0){
              lastsurv <- surv[lastMoAdu[i]] * (survRate_adu_month[i] ^ 12)/ (1.0 - survRate_adu_month[i] ^ 12)
              surv <- (((ageMo + juv$RecMonth[i]) %% 12) == 0) * surv     
              surv[lastMoAdu[i]+1] <- lastsurv
        }
        
      # Vectors of ages for adults and juvs (added for R)
        AduAges <- (ageMo >= firstMoAdu[i])
        JuvAges <- (ageMo <= lastMoJuv[i])
          
      # Sum up biomass per egg in adult group
        BioPerEgg <- sum( AduAges * surv * WageS )
        
      # Actual Eggs is Biomass of path/Biomass per recruit  
        recruits <- jpar$B_BaseRef[rAduNum]/BioPerEgg
      # KYA Also took BAB out of here (don't know if this is a juv. or adu term)
        #RzeroS[i] <- recruits
        #RzeroS[i] = recruits * exp(juv_aduBAB[i]/ STEPS_PER_YEAR);;  
       
      # Now scale numbers up to total recruits   
        NageS <- surv * recruits 
        
      # calc biomass of juvenile group using adult group input
      # adults should be just a check should be same as input adult bio 
        juvBio <- sum( JuvAges * NageS * WageS )
        aduBio <- sum( AduAges * NageS * WageS )
          
      # calculate adult assimilated consumption   
        sumK <- sum( AduAges * NageS * WWa )

      # SumK is adjusted to QB, ratio between the two is A term in respiration
      # Actually, no, since Winf could scale it up, it could be above or below
      # Q/B
    	sumK  <- jpar$FtimeQBOpt[rAduNum] * jpar$B_BaseRef[rAduNum] / sumK
        aduQB <- jpar$FtimeQBOpt[rAduNum]

      # Calculate initial juvenile assimilated consumption
        juvCons = sum(JuvAges * NageS * WWa)          

      # Scale Juvenile consumption by same A as adults
        juvCons <- juvCons *  sumK 
        juvQB   <- juvCons / juvBio 

     # Following is from SplitSetPred
     # v->stanzaBasePred[juv_JuvNum[i]] = v->stanzaPred[juv_JuvNum[i]];
     # v->stanzaBasePred[juv_AduNum[i]] = v->stanzaPred[juv_AduNum[i]]; 
     # Moved from the "second" juvenile loop - splitalpha calculation
      JuvBasePred <- sum(JuvAges * NageS * WWa) 
      AduBasePred <- sum(AduAges * NageS * WWa) 
      
      # Consumption for growth between successive age classes
        nb       <- length(WageS)
        Diff     <- WageS[2:nb] - vBM[i] * WageS[1:(nb-1)]
        Diff[nb] <- Diff[nb-1]  
      
      # Weighting to splitalpha
        SplitAlpha <- Diff * ( JuvAges * JuvBasePred/( juvQB * juvBio) +
                               AduAges * AduBasePred/( aduQB * aduBio) )
    
      #  Calculate spawning biomass as the amount of biomass over Wmat.
      ## R Comment:  LONGSTANDING Error in C code here? (ageMo > Wmat001[i]))
         SpawnBio <- sum( (( WageS > juv$Wmat001[i]) & (ageMo > juv$Amat001[i]) ) *
                        WageS * NageS /
                        (1. + exp( - ((WageS - juv$Wmat50[i]) / WmatSpread[i])   
                                   - ((ageMo - juv$Amat50[i]) / AmatSpread[i])  )))
		#cat(i,SpawnBio,"\n")
      # Qmult is for scaling to the PreyTo list so uses actual juvnum, not juvnum + 1					
        Qmult[juv$JuvNum[i]] = juvCons / (jpar$FtimeQBOpt[rJuvNum] * jpar$B_BaseRef[rJuvNum])

      # R comment:  FOLLOWING CHANGES RPAR
      # Set NoIntegrate flag (negative in both cases, different than split pools)
        jpar$NoIntegrate[rJuvNum] <- -juv$JuvNum[i]
        jpar$NoIntegrate[rAduNum] <- -juv$AduNum[i]

      # Reset juvenile B and QB in line with age stucture
        jpar$FtimeAdj[rJuvNum]  = 0.0
        jpar$FtimeAdj[rAduNum]  = 0.0
        jpar$FtimeQBOpt[rJuvNum] = juvQB
        jpar$FtimeQBOpt[rAduNum] = aduQB 
        jpar$B_BaseRef[rJuvNum]  = juvBio
    
      # KYA Spawn X is Beverton-Holt.  To turn off Beverton-Holt, set
      # SpawnX to 10000 or something similarly hight.  2.0 is half-saturation.
      # 1.000001 or so is minimum.
        k <- i+1
        jpar$SpawnX[k]            <- 10000.0
        jpar$SpawnEnergy[k]       <- 1.0 
        jpar$SpawnAllocR[k]       <- 1.0
        jpar$SpawnAllocG[k]       <- 1.0
        jpar$SpawnBio[k]       <- SpawnBio
        jpar$baseSpawnBio[k]   <- SpawnBio
        jpar$EggsStanza[k]     <- SpawnBio
        jpar$baseEggsStanza[k] <- SpawnBio
        jpar$Rbase[k] <- NageS[lastMoJuv[i]+1] * WageS[lastMoJuv[i]+1]            
      # leaving out code to rescale recruitment if recruitment is seasonal
      # calculates RscaleSplit[i] to give same annual avg as equal monthly rec
        jpar$RscaleSplit[k]     <- 1.0;

        jpar$recruits[k]       <- recruits
        jpar$RzeroS[k]         <- recruits
        jpar$stanzaGGJuv[k]    <- juvQB * juvBio / JuvBasePred
        jpar$stanzaGGAdu[k]    <- aduQB * aduBio / AduBasePred

        jpar$state_NN[rJuvNum]       <- sum(NageS * JuvAges)
        jpar$stanzaPred[rJuvNum]     <- JuvBasePred
        jpar$stanzaBasePred[rJuvNum] <- JuvBasePred

        jpar$state_NN[rAduNum]       <- sum(NageS * AduAges)
        jpar$stanzaPred[rAduNum]     <- AduBasePred
        jpar$stanzaBasePred[rAduNum] <- AduBasePred

        jpar$WageS[1:(lastMoAdu[i]+1),k]      <- WageS
        jpar$WWa[1:(lastMoAdu[i]+1),k]        <- WWa
        jpar$NageS[1:(lastMoAdu[i]+1),k]      <- NageS
        jpar$SplitAlpha[1:(lastMoAdu[i]+1),k] <- SplitAlpha
                     
    }  # END OF JUVENILE LOOP
	}
  # Reset QQ in line with new comsumtion
    #jpar$QQ <- jpar$QQ * Qmult[jpar$PreyTo]

    return(jpar) 
    
}

################################################################################
ecosim_pack<-function(rpar,YEARS=100){
   
   res <- NULL 
   attach(rpar)
   
   res$YEARS   <- YEARS

   res$spname <- rpar$spname       
   res$spnum  <- rpar$spnum          

   res$numlist <- c(NUM_GROUPS,NUM_LIVING,NUM_DEAD,NUM_GEARS,juv_N,NumPredPreyLinks,NumFishingLinks,NumDetLinks)
   names(res$numlist)<-c("NUM_GROUPS","NUM_LIVING","NUM_DEAD","NUM_GEARS","juv_N","NumPredPreyLinks","NumFishingLinks","NumDetLinks") 

   res$flags        <- c(-1,1)
   names(res$flags) <- c("BURN_YEARS","COUPLED")
   
   res$outflags        <- c(-1)
   names(res$outflags) <- c("CRASH_YEAR")

   res$ratelist <- data.frame(
                B_BaseRef,
                MzeroMort,
                UnassimRespFrac,
                ActiveRespFrac,
                FtimeAdj,
                FtimeQBOpt,
                PBopt, 
                NoIntegrate,
                HandleSelf,
                ScrambleSelf,
                PredTotWeight,
                PreyTotWeight,
                fish_Effort)
                
    res$ppind  <-   data.frame(
        PreyFrom         ,# 175   -none- numeric  
        PreyTo           # 175   -none- numeric  
    )
    
    res$pplist <- data.frame(
        QQ               ,# 175   -none- numeric  
        DD               ,# 175   -none- numeric  
        VV               ,# 175   -none- numeric  
        HandleSwitch     ,# 175   -none- numeric  
        PredPredWeight   ,# 175   -none- numeric  
        PreyPreyWeight   # 175   -none- numeric               
     )

    res$fishind <- data.frame(               
        FishFrom         ,#  89   -none- numeric  
        FishThrough      ,#  89   -none- numeric  
        FishTo           #  89   -none- numeric  
     )

   res$fishlist <- data.frame(  
        FishQ            #  89   -none- numeric  
     )

   res$detind <- data.frame(  
        DetFrom          ,#  53   -none- numeric  
        DetTo            #  53   -none- numeric     
     )
     
   res$detlist <- data.frame(  
        DetFrac          #  53   -none- numeric     
     )
 
    res$juvind <- data.frame(     
        JuvNum           ,#   5   -none- numeric  
        AduNum           ,#   5   -none- numeric 
        firstMoJuv       ,#   5   -none- numeric  
        lastMoJuv        ,#   5   -none- numeric  
        firstMoAdu       ,#   5   -none- numeric  
        lastMoAdu        #   5   -none- numeric  
     )
    res$juvlist <- data.frame(  
        SpawnBio         ,#   5   -none- numeric  
        EggsStanza       ,#   5   -none- numeric  
        stanzaGGJuv      ,#   5   -none- numeric  
        stanzaGGAdu      ,#   5   -none- numeric  
        SpawnEnergy      ,#   5   -none- numeric  
        SpawnX           ,#   5   -none- numeric  
        SpawnAllocR      ,#   5   -none- numeric  
        SpawnAllocG      ,#   5   -none- numeric  
        recruits         ,#   5   -none- numeric  
        RzeroS           ,#   5   -none- numeric  
        baseEggsStanza   ,#   5   -none- numeric  
        baseSpawnBio     ,#   5   -none- numeric  
        Rbase            ,#   5   -none- numeric  
        RecMonth         ,#   5   -none- numeric  
        VonBD            ,#   5   -none- numeric  
        aduEqAgeZ        ,#   5   -none- numeric  
        juvEqAgeZ        ,#   5   -none- numeric  
        RecPower         ,#   5   -none- numeric  
        Wmat001          ,#   5   -none- numeric  
        Wmat50           ,#   5   -none- numeric  
        Amat001          ,#   5   -none- numeric  
        Amat50           ,#   5   -none- numeric  
        WmatSpread       ,#   5   -none- numeric  
        AmatSpread       ,#   5   -none- numeric  
        vBM              ,#   5   -none- numeric  
        RscaleSplit      #   5   -none- numeric 
     )
    res$statelist <- data.frame(  
        state_BB         ,#  35   -none- numeric  
        state_Ftime      ,#  35   -none- numeric  
        state_NN         ,#  35   -none- numeric  
        stanzaPred       ,#  35   -none- numeric  
        stanzaBasePred   #  35   -none- numeric  
     )

  res$NageS      <- rpar$NageS
  res$WageS      <- rpar$WageS
  res$WWa        <- rpar$WWa 
  res$SplitAlpha <- rpar$SplitAlpha
  
  res$derivlist <- data.frame(matrix(0.0,rpar$NUM_GROUPS+1,23))
  names(res$derivlist) <- c("TotGain","TotLoss","LossPropToB","LossPropToQ","DerivT","dyt","biomeq",
                   "FoodGain","DetritalGain","FoodLoss","UnAssimLoss","ActiveRespLoss","MzeroLoss",
                   "DetritalLoss","TotDetOut","FishingLoss","FishingGain","FishingThru","preyYY","predYY",
                   "PredSuite","HandleSuite","TerminalF")  
  
  res$targlist         <-  data.frame(matrix(0.0,rpar$NUM_GROUPS+1,3))
  names (res$targlist) <- c("TARGET_BIO","TARGET_F","ALPHA")
  
  mforcemat     <- data.frame(matrix(1.0,YEARS*12+1,rpar$NUM_GROUPS+1))
  names(mforcemat) <- rpar$spname
    res$force_byprey   <- mforcemat
    res$force_bymort   <- mforcemat
    res$force_byrecs   <- mforcemat
    res$force_bysearch <- mforcemat
  
  yforcemat <- data.frame(matrix(0.0,YEARS+1,rpar$NUM_GROUPS+1))
  names(yforcemat) <- rpar$spname
    res$FORCED_FRATE <- yforcemat
    res$FORCED_CATCH <- yforcemat

  omat <- data.frame(matrix(0.0,YEARS*12+1,rpar$NUM_GROUPS+1))
  names(omat) <- rpar$spname
    res$out_BB  <- omat * 0.0
    res$out_CC  <- omat * 0.0
    res$out_SSB <- omat * 0.0              
    res$out_rec <- omat * 0.0
     
   detach(rpar) 
   return(res)
}
################################################################################
ecosim_run<-function(rpar,BYY=0,EYY=0){
	
   if ((EYY<=0)|(EYY>rpar$YEARS)) {EYY=rpar$YEARS}
   if (BYY<0)                     {BYY=0}
   if (BYY>=rpar$YEARS)           {BYY=rpar$YEARS-1}
   if (EYY<=BYY)                  {EYY=BYY+1}

  out <- rpar
  res<-.C("ecosim_run", 
    # inputs that shouldn't change
      as.integer(rpar$YEARS),
      as.integer(BYY),
      as.integer(EYY),
      as.integer(rpar$numlist),
      as.integer(rpar$flags),
      as.double(as.matrix(rpar$ratelist)),
      as.integer(as.matrix(rpar$ppind)),
      as.double(as.matrix(rpar$pplist)),
      as.integer(as.matrix(rpar$fishind)),
      as.double(as.matrix(rpar$fishlist)),
      as.integer(as.matrix(rpar$detind)),
      as.double(as.matrix(rpar$detlist)),
      as.integer(as.matrix(rpar$juvind)),
      as.double(as.matrix(rpar$juvlist)),
      as.double(as.matrix(rpar$targlist)),
    # outputs that should
      outflags       = as.integer(rpar$outflags),  
      state          = as.double(as.matrix(rpar$statelist)),
      NageS          = as.double(rpar$NageS),
      WageS          = as.double(rpar$WageS),
      WWa            = as.double(rpar$WWa),
      SplitAlpha     = as.double(rpar$SplitAlpha),
      derivlist      = as.double(as.matrix(rpar$derivlist)),
      force_byprey   = as.double(as.matrix(rpar$force_byprey)),
      force_bymort   = as.double(as.matrix(rpar$force_bymort)),
      force_byrecs   = as.double(as.matrix(rpar$force_byrecs)),
      force_bysearch = as.double(as.matrix(rpar$force_bysearch)),
      FORCED_FRATE   = as.double(as.matrix(rpar$FORCED_FRATE)),
      FORCED_CATCH   = as.double(as.matrix(rpar$FORCED_CATCH)),
      out_BB         = as.double(as.matrix(rpar$out_BB)),
      out_CC         = as.double(as.matrix(rpar$out_CC)),
      out_SSB        = as.double(as.matrix(rpar$out_SSB)), 
      out_rec        = as.double(as.matrix(rpar$out_rec)) 
)
  out$outflags[]   <- res$outflags
  out$derivlist[]  <- res$derivlist
  out$statelist[]  <- res$state
  out$NageS[]      <- res$NageS
  out$WageS[]      <- res$WageS
  out$WWa[]        <- res$WWa
  out$SplitAlpha[] <- res$SplitAlpha
  out$out_BB[]     <- res$out_BB
  out$out_CC[]     <- res$out_CC
  out$out_SSB[]     <- res$out_SSB  
  out$out_rec[]     <- res$out_rec
return(out)
}

