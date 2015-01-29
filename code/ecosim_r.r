################################################################################
# R version of ecosim (to work with ecosim.dll c code) by Kerim Aydin
#
# Version 0.04
################################################################################ 

##------------------------------------------------------------------------------ 
#Load and unload dll (unloading needed before compiling any c-code changes)
#Load
ll <- function() if(!is.loaded("ecosim_run")) dyn.load("ecosim.dll")

#Reload
rl<-function(){
  if(is.loaded("ecosim_run")) dyn.unload("ecosim.dll") 
  dyn.load("ecosim.dll")		
}

#Unload
ul<-function() if(is.loaded("ecosim_run")) dyn.unload("ecosim.dll") 	

##------------------------------------------------------------------------------ 
rmat_trans <- function(rvec){
  dim(rvec) <- c(1, nrow(rvec) * ncol(rvec))
  return(rvec)
}

##------------------------------------------------------------------------------ 
#New ecosim function that runs all steps.
ecosim <- function(Rpath, juvfile, YEARS = 100){
  #Old path_to_rates--------------------------------------------------------------------
  MSCRAMBLE      <- 2.0
  MHANDLE        <- 1000.0
  PREYSWITCH     <- 1.0
  # For SelfWts 1.0 = no overlap, 0.0 = complete overlap
  ScrambleSelfWt <- 1.0
  HandleSelfWt   <- 1.0

  rpar <- c()
  
  rpar$NUM_GROUPS <- Rpath$NUM_GROUPS
  rpar$NUM_LIVING <- Rpath$NUM_LIVING
  rpar$NUM_DEAD   <- Rpath$NUM_DEAD
  rpar$NUM_GEARS  <- Rpath$NUM_GEARS
  
  rpar$spname     <- c("Outside", Rpath$Group)
  rpar$spnum      <- 0:length(Rpath$BB) 

  #Energetics for Living and Dead Groups
  #Reference biomass for calculating YY
  rpar$B_BaseRef <- c(1.0, Rpath$BB) 
  #Mzero proportional to (1-EE)
  rpar$MzeroMort <- c(0.0, Rpath$PB * (1.0 - Rpath$EE)) 
  #Unassimilated is the proportion of CONSUMPTION that goes to detritus.  
  rpar$UnassimRespFrac <- c(0.0, Rpath$GS);
  #Active respiration is proportion of CONSUMPTION that goes to "heat"
  #Passive respiration/ VonB adjustment is left out here
  rpar$ActiveRespFrac <-  c(0.0, ifelse(Rpath$QB > 0, 
                                        1.0 - (Rpath$PB / Rpath$QB) - Rpath$GS, 
                                        0.0))
  #Ftime related parameters
  rpar$FtimeAdj   <- rep(0.0, length(rpar$B_BaseRef))
  rpar$FtimeQBOpt <-   c(1.0, Rpath$QB)
  rpar$PBopt      <-   c(1.0, Rpath$PB)           

  #Fishing Effort defaults to 0 for non-gear, 1 for gear
  rpar$fish_Effort <- ifelse(rpar$spnum <= rpar$NUM_LIVING + rpar$NUM_DEAD,
                             0.0,
                             1.0) 
    
  #NoIntegrate
  STEPS_PER_YEAR  <- 12
  STEPS_PER_MONTH <- 1
  rpar$NoIntegrate <- ifelse(c(0, Rpath$PB) / 
                               (1.0 - rpar$ActiveRespFrac - rpar$UnassimRespFrac) > 
                               2 * STEPS_PER_YEAR * STEPS_PER_MONTH, 
                             0, 
                             rpar$spnum)  

  #Pred/Prey defaults
  rpar$HandleSelf   <- rep(HandleSelfWt,   Rpath$NUM_GROUPS + 1)
  rpar$ScrambleSelf <- rep(ScrambleSelfWt, Rpath$NUM_GROUPS + 1)
  
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
  rpar$PreyFrom <- c(primFrom[primTo > 0], preyfrom [predpreyQ > 0])
  rpar$PreyTo   <- c(primTo  [primTo > 0], preyto   [predpreyQ > 0])
  rpar$QQ       <- c(primQ   [primTo > 0], predpreyQ[predpreyQ > 0])             	
  
  numpredprey <- length(rpar$QQ)

  rpar$DD <- rep(MHANDLE,   numpredprey)
  rpar$VV <- rep(MSCRAMBLE, numpredprey)

  #NOTE:  Original in C didn't set handleswitch for primary production groups.  Error?
  #probably not when group 0 biomass doesn't change from 1.
  rpar$HandleSwitch <- rep(PREYSWITCH, numpredprey)

  #scramble combined prey pools
  Btmp <- rpar$B_BaseRef
  py   <- rpar$PreyFrom + 1.0
  pd   <- rpar$PreyTo + 1.0
  VV   <- rpar$VV * rpar$QQ / Btmp[py]
  AA   <- (2.0 * rpar$QQ * VV) / (VV * Btmp[pd] * Btmp[py] - rpar$QQ * Btmp[pd])
  rpar$PredPredWeight <- AA * Btmp[pd] 
  rpar$PreyPreyWeight <- AA * Btmp[py] 

  rpar$PredTotWeight <- rep(0, length(rpar$B_BaseRef))
  rpar$PreyTotWeight <- rep(0, length(rpar$B_BaseRef))
  
  for(links in 1:numpredprey){
    rpar$PredTotWeight[py[links]] <- rpar$PredTotWeight[py[links]] + rpar$PredPredWeight[links]
    rpar$PreyTotWeight[pd[links]] <- rpar$PreyTotWeight[pd[links]] + rpar$PreyPreyWeight[links]    
  }  
  #rpar$PredTotWeight[]   <- as.numeric(tapply(rpar$PredPredWeight,py,sum))
  #rpar$PreyTotWeight[]   <- as.numeric(tapply(rpar$PreyPreyWeight,pd,sum))

  rpar$PredPredWeight <- rpar$PredPredWeight/rpar$PredTotWeight[py]
  rpar$PreyPreyWeight <- rpar$PreyPreyWeight/rpar$PreyTotWeight[pd]

  rpar$NumPredPreyLinks <- numpredprey
  rpar$PreyFrom       <- c(0, rpar$PreyFrom)
  rpar$PreyTo         <- c(0, rpar$PreyTo)
  rpar$QQ             <- c(0, rpar$QQ)
  rpar$DD             <- c(0, rpar$DD)
  rpar$VV             <- c(0, rpar$VV) 
  rpar$HandleSwitch   <- c(0, rpar$HandleSwitch) 
  rpar$PredPredWeight <- c(0, rpar$PredPredWeight)
  rpar$PreyPreyWeight <- c(0, rpar$PreyPreyWeight)
  
  #catchlinks
  fishfrom    <- row(Rpath$Catch)                      
  fishthrough <- col(Rpath$Catch) + (Rpath$NUM_LIVING + Rpath$NUM_DEAD)
  fishcatch   <- Rpath$Catch
  fishto      <- fishfrom * 0
  
  if(sum(fishcatch) > 0){
    rpar$FishFrom    <- fishfrom   [fishcatch > 0]
    rpar$FishThrough <- fishthrough[fishcatch > 0]
    rpar$FishQ       <- fishcatch  [fishcatch > 0] / rpar$B_BaseRef[rpar$FishFrom + 1]  
    rpar$FishTo      <- fishto     [fishcatch > 0]
  }
  #discard links
  
  for(d in 1:Rpath$NUM_DEAD){
    detfate <- Rpath$DetFate[(Rpath$NUM_LIVING + Rpath$NUM_DEAD + 1):Rpath$NUM_GROUPS, d]
    detmat  <- t(matrix(detfate, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
   
    fishfrom    <-  row(Rpath$Discards)                      
    fishthrough <-  col(Rpath$Discards) + (Rpath$NUM_LIVING + Rpath$NUM_DEAD)
    fishto      <-  t(matrix(Rpath$NUM_LIVING + d, Rpath$NUM_GEAR, Rpath$NUM_GROUPS))
    fishcatch   <-  Rpath$Discards * detmat
    if(sum(fishcatch) > 0){
      rpar$FishFrom    <- c(rpar$FishFrom,    fishfrom   [fishcatch > 0])
      rpar$FishThrough <- c(rpar$FishThrough, fishthrough[fishcatch > 0])
      ffrom <- fishfrom[fishcatch > 0]
      rpar$FishQ       <- c(rpar$FishQ,  fishcatch[fishcatch > 0] / rpar$B_BaseRef[ffrom + 1])  
      rpar$FishTo      <- c(rpar$FishTo, fishto   [fishcatch > 0])
   }
  } 
  
  rpar$NumFishingLinks <- length(rpar$FishFrom)  
  rpar$FishFrom        <- c(0, rpar$FishFrom)
  rpar$FishThrough     <- c(0, rpar$FishThrough)
  rpar$FishQ           <- c(0, rpar$FishQ)  
  rpar$FishTo          <- c(0, rpar$FishTo)   

# Unwound discard fate for living groups
  detfrac <- Rpath$DetFate[1:(Rpath$NUM_LIVING + Rpath$NUM_DEAD), ]
  detfrom <- row(detfrac)
  detto   <- col(detfrac) + Rpath$NUM_LIVING
  
  detout <- 1 - rowSums(Rpath$DetFate[1:(Rpath$NUM_LIVING + Rpath$NUM_DEAD), ])
  dofrom <- 1:length(detout)
  doto   <- rep(0, length(detout))
  
  rpar$DetFrac <- c(0, detfrac[detfrac > 0], detout[detout > 0])
  rpar$DetFrom <- c(0, detfrom[detfrac > 0], dofrom[detout > 0])
  rpar$DetTo   <- c(0, detto  [detfrac > 0], doto  [detout > 0])
  
  rpar$NumDetLinks <- length(rpar$DetFrac) - 1

  rpar$state_BB    <- rpar$B_BaseRef
  rpar$state_Ftime <- rep(1, length(Rpath$BB) + 1)
  
  #Old initialize_stanzas------------------------------------------------------------------------
  #Initialize juvenile adult or "stanza" age structure in sim   
  MAX_MONTHS_STANZA <- 400
  MIN_REC_FACTOR    <- 6.906754779  #// This is ln(1/0.001 - 1) used to set min. logistic matuirty to 0.001

  juv        <- read.csv(juvfile)
  rpar$juv_N <- length(juv$JuvNum) 
    
  #fill monthly vectors for each species, rel weight and consumption at age
  #Loops for weight and numbers
  rpar$WageS      <- rep(0, MAX_MONTHS_STANZA + 1 * rpar$juv_N + 1)
  rpar$WWa        <- rep(0, MAX_MONTHS_STANZA + 1 * rpar$juv_N + 1)
  rpar$NageS      <- rep(0, MAX_MONTHS_STANZA + 1 * rpar$juv_N + 1)
  rpar$SplitAlpha <- rep(0, MAX_MONTHS_STANZA + 1 * rpar$juv_N + 1)
  
  rpar$state_NN       <-rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$stanzaPred     <-rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$stanzaBasePred <-rep(0.0, rpar$NUM_GROUPS + 1)

  rpar$SpawnBio       <- rep(0.0, rpar$juv_N + 1)
  rpar$EggsStanza     <- rep(0.0, rpar$juv_N + 1)
  rpar$stanzaGGJuv    <- rep(0.0, rpar$juv_N + 1)
  rpar$stanzaGGAdu    <- rep(0.0, rpar$juv_N + 1)
  #rpar$drawn_K[MAX_SPLIT + 1]
  #rpar$drawn_AduZ[MAX_SPLIT + 1]
  #rpar$drawn_JuvZ[MAX_SPLIT + 1]
  rpar$SpawnEnergy    <- rep(0.0, rpar$juv_N + 1)
  rpar$SpawnX         <- rep(0.0, rpar$juv_N + 1)
  rpar$SpawnAllocR    <- rep(0.0, rpar$juv_N + 1)
  rpar$SpawnAllocG    <- rep(0.0, rpar$juv_N + 1)

  rpar$recruits       <- rep(0.0, rpar$juv_N + 1)
  rpar$RzeroS         <- rep(0.0, rpar$juv_N + 1)
  rpar$baseEggsStanza <- rep(0.0, rpar$juv_N + 1)

  rpar$baseSpawnBio   <- rep(0.0, rpar$juv_N + 1)
  rpar$Rbase          <- rep(0.0, rpar$juv_N + 1)
  rpar$RscaleSplit    <- rep(0.0, rpar$juv_N + 1)
    
  rpar$JuvNum         <- c(0, juv$JuvNum)
  rpar$AduNum         <- c(0, juv$AduNum)
  rpar$RecMonth       <- c(0, juv$RecMonth)
  rpar$VonBD          <- c(0, juv$VonBD)
	rpar$aduEqAgeZ      <- c(0, juv$JuvZ_BAB)
	rpar$juvEqAgeZ      <- c(0, juv$AduZ_BAB)
  rpar$RecPower       <- c(0, juv$RecPower)
  rpar$Wmat001        <- c(0, juv$Wmat001)
  rpar$Wmat50         <- c(0, juv$Wmat50)
  rpar$Amat001        <- c(0, juv$Amat001)
  rpar$Amat50         <- c(0, juv$Amat50)
    
  #first calculate the number of months in each age group.  
  #For end stanza (plus group), capped in EwE @400, this is 90% of 
  #relative max wt, uses generalized exponent D
  firstMoJuv <- rep(0, rpar$juv_N) 
  lastMoJuv  <- juv$RecAge * STEPS_PER_YEAR - 1
  firstMoAdu <- juv$RecAge * STEPS_PER_YEAR
  lastMoAdu  <- floor(log(1.0 - (0.9 ^ (1.0 - juv$VonBD))) /
                        (-3.0 * juv$VonBK * (1.0 - juv$VonBD) / STEPS_PER_YEAR))
  lastMoAdu  <- ifelse(lastMoAdu > MAX_MONTHS_STANZA, 
                       MAX_MONTHS_STANZA, 
                       lastMoAdu)

  #Energy required to grow a unit in weight (scaled to Winf = 1)
  vBM <- 1.0 - 3.0 * juv$VonBK / STEPS_PER_YEAR

  #Set Q multiplier to 1 (used to rescale consumption rates) 
  Qmult <- rep(1.0, rpar$NUM_GROUPS) 
 
  #Survival rates for juveniles and adults
  #KYA Switched survival rate out of BAB == assume all on adults
  survRate_juv_month <- exp(-(juv$JuvZ_BAB) / STEPS_PER_YEAR)
  survRate_adu_month <- exp(-(juv$AduZ_BAB) / STEPS_PER_YEAR)

  WmatSpread <- -(juv$Wmat001 - juv$Wmat50) / MIN_REC_FACTOR
  AmatSpread <- -(juv$Amat001 - juv$Amat50) / MIN_REC_FACTOR

  #Save numbers that can be saved at this point
  rpar$WmatSpread  <- c(0, WmatSpread)    
  rpar$AmatSpread  <- c(0, AmatSpread)
  rpar$vBM         <- c(0, vBM)
  rpar$firstMoJuv  <- c(0, firstMoJuv)
  rpar$lastMoJuv   <- c(0, lastMoJuv)
  rpar$firstMoAdu  <- c(0, firstMoAdu)
  rpar$lastMoAdu   <- c(0, lastMoAdu)       
    
	if(rpar$juv_N>0){                        
    for(i in 1:rpar$juv_N){
	    #R rpar lookup number for juvs and adults is shifted by 1
      rJuvNum <- juv$JuvNum[i] + 1
      rAduNum <- juv$AduNum[i] + 1
                             
      #vector of ages in months from month 0
      ageMo <- 0:lastMoAdu[i]
        
      #Weight and consumption at age
      WageS <- (1.0 - exp(-3 * juv$VonBK[i] * 
                            (1 - juv$VonBD[i]) * 
                            ageMo / STEPS_PER_YEAR)) ^ (1.0 / (1.0 - juv$VonBD[i]))
      WWa <- WageS ^ juv$VonBD[i]
        
      #Numbers per recruit
      surv <- ifelse(ageMo <= firstMoAdu[i], 
                     survRate_juv_month[i] ^ ageMo,
                     survRate_juv_month[i] ^ firstMoAdu[i] * 
                       survRate_adu_month[i] ^ (ageMo - firstMoAdu[i])
                      )
      #plus group
      surv[lastMoAdu[i] + 1] <- surv[lastMoAdu[i]] * survRate_adu_month[i] / 
        (1.0 - survRate_adu_month[i])
            
      #correct if using monthly pulse recruitment
      if(juv$RecMonth[i] > 0){
        lastsurv               <- surv[lastMoAdu[i]] * (survRate_adu_month[i] ^ 12) / 
                                   (1.0 - survRate_adu_month[i] ^ 12)
        surv                   <- (((ageMo + juv$RecMonth[i]) %% 12) == 0) * surv     
        surv[lastMoAdu[i] + 1] <- lastsurv
        }
        
      #Vectors of ages for adults and juvs (added for R)
      AduAges <- (ageMo >= firstMoAdu[i])
      JuvAges <- (ageMo <= lastMoJuv [i])
          
      #Sum up biomass per egg in adult group
      BioPerEgg <- sum(AduAges * surv * WageS)
        
      #Actual Eggs is Biomass of Rpath/Biomass per recruit  
      recruits <- rpar$B_BaseRef[rAduNum] / BioPerEgg
      #KYA Also took BAB out of here (don't know if this is a juv. or adu term)
      #RzeroS[i] <- recruits
      #RzeroS[i] = recruits * exp(juv_aduBAB[i]/ STEPS_PER_YEAR);;  
       
      #Now scale numbers up to total recruits   
      NageS <- surv * recruits 
        
      #calc biomass of juvenile group using adult group input
      #adults should be just a check should be same as input adult bio 
      juvBio <- sum(JuvAges * NageS * WageS)
      aduBio <- sum(AduAges * NageS * WageS)
          
      #calculate adult assimilated consumption   
      sumK <- sum(AduAges * NageS * WWa)

      #SumK is adjusted to QB, ratio between the two is A term in respiration
      #Actually, no, since Winf could scale it up, it could be above or below
      #Q/B
    	sumK  <- rpar$FtimeQBOpt[rAduNum] * rpar$B_BaseRef[rAduNum] / sumK
      aduQB <- rpar$FtimeQBOpt[rAduNum]

      #Calculate initial juvenile assimilated consumption
      juvCons <- sum(JuvAges * NageS * WWa)          

      #Scale Juvenile consumption by same A as adults
      juvCons <- juvCons *  sumK 
      juvQB   <- juvCons / juvBio 

      #Following is from SplitSetPred
      #v->stanzaBasePred[juv_JuvNum[i]] = v->stanzaPred[juv_JuvNum[i]];
      #v->stanzaBasePred[juv_AduNum[i]] = v->stanzaPred[juv_AduNum[i]]; 
      #Moved from the "second" juvenile loop - splitalpha calculation
      JuvBasePred <- sum(JuvAges * NageS * WWa) 
      AduBasePred <- sum(AduAges * NageS * WWa) 
      
      #Consumption for growth between successive age classes
      nb       <- length(WageS)
      Diff     <- WageS[2:nb] - vBM[i] * WageS[1:(nb - 1)]
      Diff[nb] <- Diff[nb - 1]  
      
      #Weighting to splitalpha
      SplitAlpha <- Diff * (JuvAges * JuvBasePred / (juvQB * juvBio) +
                              AduAges * AduBasePred / (aduQB * aduBio))
    
      #Calculate spawning biomass as the amount of biomass over Wmat.
      ##R Comment:  LONGSTANDING Error in C code here? (ageMo > Wmat001[i]))
      SpawnBio <- sum(((WageS > juv$Wmat001[i]) & (ageMo > juv$Amat001[i])) *
                        WageS * NageS /
                        (1 + exp(-((WageS - juv$Wmat50[i]) / WmatSpread[i])   
                                   -((ageMo - juv$Amat50[i]) / AmatSpread[i]))))
		  #cat(i,SpawnBio,"\n")
      #Qmult is for scaling to the PreyTo list so uses actual juvnum, not juvnum + 1					
      Qmult[juv$JuvNum[i]] <- juvCons / (rpar$FtimeQBOpt[rJuvNum] * rpar$B_BaseRef[rJuvNum])

      #R comment:  FOLLOWING CHANGES RPAR
      #Set NoIntegrate flag (negative in both cases, different than split pools)
      rpar$NoIntegrate[rJuvNum] <- -juv$JuvNum[i]
      rpar$NoIntegrate[rAduNum] <- -juv$AduNum[i]

      #Reset juvenile B and QB in line with age stucture
      rpar$FtimeAdj  [rJuvNum] <- 0.0
      rpar$FtimeAdj  [rAduNum] <- 0.0
      rpar$FtimeQBOpt[rJuvNum] <- juvQB
      rpar$FtimeQBOpt[rAduNum] <- aduQB 
      rpar$B_BaseRef [rJuvNum] <- juvBio
    
      #KYA Spawn X is Beverton-Holt.  To turn off Beverton-Holt, set
      #SpawnX to 10000 or something similarly hight.  2.0 is half-saturation.
      #1.000001 or so is minimum.
      k <- i + 1
      rpar$SpawnX[k]         <- 10000.0
      rpar$SpawnEnergy[k]    <- 1.0 
      rpar$SpawnAllocR[k]    <- 1.0
      rpar$SpawnAllocG[k]    <- 1.0
      rpar$SpawnBio[k]       <- SpawnBio
      rpar$baseSpawnBio[k]   <- SpawnBio
      rpar$EggsStanza[k]     <- SpawnBio
      rpar$baseEggsStanza[k] <- SpawnBio
      rpar$Rbase[k]          <- NageS[lastMoJuv[i] + 1] * WageS[lastMoJuv[i] + 1]            
      #leaving out code to rescale recruitment if recruitment is seasonal
      #calculates RscaleSplit[i] to give same annual avg as equal monthly rec
      rpar$RscaleSplit[k]     <- 1.0

      rpar$recruits[k]       <- recruits
      rpar$RzeroS[k]         <- recruits
      rpar$stanzaGGJuv[k]    <- juvQB * juvBio / JuvBasePred
      rpar$stanzaGGAdu[k]    <- aduQB * aduBio / AduBasePred

      rpar$state_NN[rJuvNum]       <- sum(NageS * JuvAges)
      rpar$stanzaPred[rJuvNum]     <- JuvBasePred
      rpar$stanzaBasePred[rJuvNum] <- JuvBasePred

      rpar$state_NN[rAduNum]       <- sum(NageS * AduAges)
      rpar$stanzaPred[rAduNum]     <- AduBasePred
      rpar$stanzaBasePred[rAduNum] <- AduBasePred

      rpar$WageS     [1:(lastMoAdu[i] + 1), k] <- WageS
      rpar$WWa       [1:(lastMoAdu[i] + 1), k] <- WWa
      rpar$NageS     [1:(lastMoAdu[i] + 1), k] <- NageS
      rpar$SplitAlpha[1:(lastMoAdu[i] + 1), k] <- SplitAlpha
                     
    }  #END OF JUVENILE LOOP
	}
  #Reset QQ in line with new comsumtion
  #rpar$QQ <- rpar$QQ * Qmult[rpar$PreyTo]

  #Old ecosim_pack-------------------------------------------------------------------- 
  rpar$YEARS  <- YEARS

  
  rpar$BURN_YEARS <- -1
  rpar$COUPLED    <-  1
 
  rpar$CRASH_YEAR <- -1
  
  #derivlist
  rpar$TotGain        <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$TotLoss        <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$LossPropToB    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$LossPropToQ    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$DerivT         <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$dyt            <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$biomeq         <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$FoodGain       <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$DetritalGain   <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$FoodLoss       <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$UnAssimLoss    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$ActiveRespLoss <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$MzeroLoss      <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$DetritalLoss   <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$TotDetOut      <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$FishingLoss    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$FishingGain    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$FishingThru    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$preyYY         <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$predYY         <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$PredSuite      <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$HandleSuite    <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$TerminalF      <- rep(0.0, rpar$NUM_GROUPS + 1)
 
  
  #targlist     
  rpar$TARGET_BIO <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$TARGET_F   <- rep(0.0, rpar$NUM_GROUPS + 1)
  rpar$ALPHA      <- rep(0.0, rpar$NUM_GROUPS + 1)
  
  mforcemat           <- data.frame(matrix(1.0, YEARS * 12 + 1, rpar$NUM_GROUPS + 1))
  names(mforcemat)    <- rpar$spname
  rpar$force_byprey   <- mforcemat
  rpar$force_bymort   <- mforcemat
  rpar$force_byrecs   <- mforcemat
  rpar$force_bysearch <- mforcemat
  
  yforcemat         <- data.frame(matrix(0.0, YEARS + 1, rpar$NUM_GROUPS + 1))
  names(yforcemat)  <- rpar$spname
  rpar$FORCED_FRATE <- yforcemat
  rpar$FORCED_CATCH <- yforcemat

  omat           <- data.frame(matrix(0.0, YEARS * 12 + 1, rpar$NUM_GROUPS + 1))
  names(omat)    <- rpar$spname
  rpar$out_BB  <- omat 
  rpar$out_CC  <- omat
  rpar$out_SSB <- omat              
  rpar$out_rec <- omat
  
  class(rpar) <- 'Rpath.sim'
  return(rpar)

}
 

##------------------------------------------------------------------------------ 

ecosim_run <- function(rpar, BYY = 0, EYY = 0){
  if((EYY <= 0) | (EYY > rpar$YEARS)) EYY <- rpar$YEARS
  if(BYY <0)                          BYY <- 0
  if(BYY >= rpar$YEARS)               BYY <- rpar$YEARS - 1
  if(EYY <= BYY)                      EYY <- BYY + 1

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
    out_rec        = as.double(as.matrix(rpar$out_rec)))
  
  out$outflags[]   <- res$outflags
  out$derivlist[]  <- res$derivlist
  out$statelist[]  <- res$state
  out$NageS[]      <- res$NageS
  out$WageS[]      <- res$WageS
  out$WWa[]        <- res$WWa
  out$SplitAlpha[] <- res$SplitAlpha
  out$out_BB[]     <- res$out_BB
  out$out_CC[]     <- res$out_CC
  out$out_SSB[]    <- res$out_SSB  
  out$out_rec[]    <- res$out_rec

  return(out)
}

