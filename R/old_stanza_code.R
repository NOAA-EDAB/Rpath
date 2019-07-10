
# 
# simpar$state_BB    <- simpar$B_BaseRef
# simpar$state_Ftime <- rep(1, length(Rpath$BB) + 1)
# return(simpar)
# }
# 
# #'Initialization of stanzas for ecosim
# #'
# #'Sets up the initial conditions for groups with multiple stanzas.
# #'
# #'@family Rpath functions
# #'
# #'@param simpar List of ecosim parameters created by ecosim.rates
# #'@param juvfile Comma deliminated file with multi-stanza parameters. 
# #'@param steps_yr Number of time steps per year.
# #'@param max_months Maximum number of months that a stanza can survive.
# #'@param min_rec_factor Used to set minimum logistic maturity to 0.001.  Default value is ln(1 / 0.001 - 1).
# #'
# #'@return Returns an Rpath.sim object that can be supplied to the ecosim.run function.
# #'@export 
# ecosim.stanzas <- function(simpar, juvfile, steps_yr = 12, max_months = 400, min_rec_factor = 6.906754779){
#   #Juvenile Parameters
#   if(is.character(juvfile)){
#     juv  <- as.data.table(read.csv(juvfile))
#   } else {
#     juv <- as.data.table(juvfile)
#   }
#   simpar$juv_N <- length(juv$JuvNum) 
#   
#   #fill monthly vectors for each species, rel weight and consumption at age
#   #Loops for weight and numbers
#   simpar$WageS      <- matrix(0, max_months + 1, simpar$juv_N + 1)
#   simpar$WWa        <- matrix(0, max_months + 1, simpar$juv_N + 1)
#   simpar$NageS      <- matrix(0, max_months + 1, simpar$juv_N + 1)
#   simpar$SplitAlpha <- matrix(0, max_months + 1, simpar$juv_N + 1)
#   
#   simpar$state_NN       <-rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$stanzaPred     <-rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$stanzaBasePred <-rep(0.0, simpar$NUM_GROUPS + 1)
#   
#   simpar$SpawnBio       <- rep(0.0, simpar$juv_N + 1)
#   simpar$EggsStanza     <- rep(0.0, simpar$juv_N + 1)
#   simpar$stanzaGGJuv    <- rep(0.0, simpar$juv_N + 1)
#   simpar$stanzaGGAdu    <- rep(0.0, simpar$juv_N + 1)
#   #simpar$drawn_K[MAX_SPLIT + 1]
#   #simpar$drawn_AduZ[MAX_SPLIT + 1]
#   #simpar$drawn_JuvZ[MAX_SPLIT + 1]
#   simpar$SpawnEnergy    <- rep(0.0, simpar$juv_N + 1)
#   simpar$SpawnX         <- rep(0.0, simpar$juv_N + 1)
#   simpar$SpawnAllocR    <- rep(0.0, simpar$juv_N + 1)
#   simpar$SpawnAllocG    <- rep(0.0, simpar$juv_N + 1)
#   
#   simpar$recruits       <- rep(0.0, simpar$juv_N + 1)
#   simpar$RzeroS         <- rep(0.0, simpar$juv_N + 1)
#   simpar$baseEggsStanza <- rep(0.0, simpar$juv_N + 1)
#   
#   simpar$baseSpawnBio   <- rep(0.0, simpar$juv_N + 1)
#   simpar$Rbase          <- rep(0.0, simpar$juv_N + 1)
#   simpar$RscaleSplit    <- rep(0.0, simpar$juv_N + 1)
#   
#   simpar$JuvNum         <- c(0, juv$JuvNum)
#   simpar$AduNum         <- c(0, juv$AduNum)
#   simpar$RecMonth       <- c(0, juv$RecMonth)
#   simpar$VonBD          <- c(0, juv$VonBD)
#   simpar$aduEqAgeZ      <- c(0, juv$JuvZ_BAB)
#   simpar$juvEqAgeZ      <- c(0, juv$AduZ_BAB)
#   simpar$RecPower       <- c(0, juv$RecPower)
#   simpar$Wmat001        <- c(0, juv$Wmat001)
#   simpar$Wmat50         <- c(0, juv$Wmat50)
#   simpar$Amat001        <- c(0, juv$Amat001)
#   simpar$Amat50         <- c(0, juv$Amat50)
#   
#   #first calculate the number of months in each age group.  
#   #For end stanza (plus group), capped in EwE @400, this is 90% of 
#   #relative max wt, uses generalized exponent D
#   firstMoJuv <- rep(0, simpar$juv_N) 
#   lastMoJuv  <- juv$RecAge * steps_yr - 1
#   firstMoAdu <- juv$RecAge * steps_yr
#   lastMoAdu  <- floor(log(1.0 - (0.9 ^ (1.0 - juv$VonBD))) / 
#                         (-3.0 * juv$VonBK * (1.0 - juv$VonBD) / 
#                            steps_yr))
#   lastMoAdu  <- ifelse(lastMoAdu > max_months, 
#                        max_months, 
#                        lastMoAdu)
#   
#   #Energy required to grow a unit in weight (scaled to Winf = 1)
#   vBM <- 1.0 - 3.0 * juv$VonBK / steps_yr
#   
#   #Set Q multiplier to 1 (used to rescale consumption rates) 
#   Qmult <- rep(1.0, simpar$NUM_GROUPS) 
#   
#   #Survival rates for juveniles and adults
#   #KYA Switched survival rate out of BAB == assume all on adults
#   survRate_juv_month <- exp(-(juv$JuvZ_BAB) / steps_yr)
#   survRate_adu_month <- exp(-(juv$AduZ_BAB) / steps_yr)
#   
#   WmatSpread <- -(juv$Wmat001 - juv$Wmat50) / min_rec_factor
#   AmatSpread <- -(juv$Amat001 - juv$Amat50) / min_rec_factor
#   
#   #Save numbers that can be saved at this point
#   simpar$WmatSpread  <- c(0, WmatSpread)    
#   simpar$AmatSpread  <- c(0, AmatSpread)
#   simpar$vBM         <- c(0, vBM)
#   simpar$firstMoJuv  <- c(0, firstMoJuv)
#   simpar$lastMoJuv   <- c(0, lastMoJuv)
#   simpar$firstMoAdu  <- c(0, firstMoAdu)
#   simpar$lastMoAdu   <- c(0, lastMoAdu)       
#   
#   if(simpar$juv_N > 0){                        
#     for(i in 1:simpar$juv_N){
#       #R simpar lookup number for juvs and adults is shifted by 1
#       rJuvNum <- juv$JuvNum[i] + 1
#       rAduNum <- juv$AduNum[i] + 1
#       
#       #vector of ages in months from month 0
#       ageMo <- 0:lastMoAdu[i]
#       
#       #Weight and consumption at age
#       WageS <- (1.0 - exp(-3 * juv$VonBK[i] * 
#                             (1 - juv$VonBD[i]) * 
#                             ageMo / steps_yr)) ^ (1.0 / (1.0 - juv$VonBD[i]))
#       WWa <- WageS ^ juv$VonBD[i]
#       
#       #Numbers per recruit
#       surv <- ifelse(ageMo <= firstMoAdu[i], 
#                      survRate_juv_month[i] ^ ageMo,
#                      survRate_juv_month[i] ^ firstMoAdu[i] * 
#                        survRate_adu_month[i] ^ (ageMo - firstMoAdu[i])
#       )
#       #plus group
#       surv[lastMoAdu[i] + 1] <- surv[lastMoAdu[i]] * survRate_adu_month[i] / 
#         (1.0 - survRate_adu_month[i])
#       
#       #correct if using monthly pulse recruitment
#       if(juv$RecMonth[i] > 0){
#         lastsurv               <- surv[lastMoAdu[i]] * (survRate_adu_month[i] ^ 12) / 
#           (1.0 - survRate_adu_month[i] ^ 12)
#         surv                   <- (((ageMo + juv$RecMonth[i]) %% 12) == 0) * surv     
#         surv[lastMoAdu[i] + 1] <- lastsurv
#       }
#       
#       #Vectors of ages for adults and juvs (added for R)
#       AduAges <- (ageMo >= firstMoAdu[i])
#       JuvAges <- (ageMo <= lastMoJuv [i])
#       
#       #Sum up biomass per egg in adult group
#       BioPerEgg <- sum(AduAges * surv * WageS)
#       
#       #Actual Eggs is Biomass of Rpath/Biomass per recruit  
#       recruits <- simpar$B_BaseRef[rAduNum] / BioPerEgg
#       #KYA Also took BAB out of here (don't know if this is a juv. or adu term)
#       #RzeroS[i] <- recruits
#       #RzeroS[i] = recruits * exp(juv_aduBAB[i]/ STEPS_PER_YEAR);;  
#       
#       #Now scale numbers up to total recruits   
#       NageS <- surv * recruits 
#       
#       #calc biomass of juvenile group using adult group input
#       #adults should be just a check should be same as input adult bio 
#       juvBio <- sum(JuvAges * NageS * WageS)
#       aduBio <- sum(AduAges * NageS * WageS)
#       
#       #calculate adult assimilated consumption   
#       sumK <- sum(AduAges * NageS * WWa)
#       
#       #SumK is adjusted to QB, ratio between the two is A term in respiration
#       #Actually, no, since Winf could scale it up, it could be above or below
#       #Q/B
#       sumK  <- simpar$FtimeQBOpt[rAduNum] * simpar$B_BaseRef[rAduNum] / sumK
#       aduQB <- simpar$FtimeQBOpt[rAduNum]
#       
#       #Calculate initial juvenile assimilated consumption
#       juvCons <- sum(JuvAges * NageS * WWa)          
#       
#       #Scale Juvenile consumption by same A as adults
#       juvCons <- juvCons *  sumK 
#       juvQB   <- juvCons / juvBio 
#       
#       #Following is from SplitSetPred
#       #v->stanzaBasePred[juv_JuvNum[i]] = v->stanzaPred[juv_JuvNum[i]];
#       #v->stanzaBasePred[juv_AduNum[i]] = v->stanzaPred[juv_AduNum[i]]; 
#       #Moved from the "second" juvenile loop - splitalpha calculation
#       JuvBasePred <- sum(JuvAges * NageS * WWa) 
#       AduBasePred <- sum(AduAges * NageS * WWa) 
#       
#       #Consumption for growth between successive age classes
#       nb       <- length(WageS)
#       Diff     <- WageS[2:nb] - vBM[i] * WageS[1:(nb - 1)]
#       Diff[nb] <- Diff[nb - 1]  
#       
#       #Weighting to splitalpha
#       SplitAlpha <- Diff * (JuvAges * JuvBasePred / (juvQB * juvBio) +
#                               AduAges * AduBasePred / (aduQB * aduBio))
#       
#       #Calculate spawning biomass as the amount of biomass over Wmat.
#       ##R Comment:  LONGSTANDING Error in C code here? (ageMo > Wmat001[i]))
#       SpawnBio <- sum(((WageS > juv$Wmat001[i]) & (ageMo > juv$Amat001[i])) *
#                         WageS * NageS /
#                         (1 + exp(-((WageS - juv$Wmat50[i]) / WmatSpread[i])   
#                                  -((ageMo - juv$Amat50[i]) / AmatSpread[i]))))
#       #cat(i,SpawnBio,"\n")
#       #Qmult is for scaling to the PreyTo list so uses actual juvnum, not juvnum + 1					
#       Qmult[juv$JuvNum[i]] <- juvCons / (simpar$FtimeQBOpt[rJuvNum] * simpar$B_BaseRef[rJuvNum])
#       
#       #R comment:  FOLLOWING CHANGES simpar
#       #Set NoIntegrate flag (negative in both cases, different than split pools)
#       simpar$NoIntegrate[rJuvNum] <- -juv$JuvNum[i]
#       simpar$NoIntegrate[rAduNum] <- -juv$AduNum[i]
#       
#       #Reset juvenile B and QB in line with age stucture
#       simpar$FtimeAdj  [rJuvNum] <- 0.0
#       simpar$FtimeAdj  [rAduNum] <- 0.0
#       simpar$FtimeQBOpt[rJuvNum] <- juvQB
#       simpar$FtimeQBOpt[rAduNum] <- aduQB 
#       simpar$B_BaseRef [rJuvNum] <- juvBio
#       
#       #KYA Spawn X is Beverton-Holt.  To turn off Beverton-Holt, set
#       #SpawnX to 10000 or something similarly hight.  2.0 is half-saturation.
#       #1.000001 or so is minimum.
#       k <- i + 1
#       simpar$SpawnX[k]         <- 10000.0
#       simpar$SpawnEnergy[k]    <- 1.0 
#       simpar$SpawnAllocR[k]    <- 1.0
#       simpar$SpawnAllocG[k]    <- 1.0
#       simpar$SpawnBio[k]       <- SpawnBio
#       simpar$baseSpawnBio[k]   <- SpawnBio
#       simpar$EggsStanza[k]     <- SpawnBio
#       simpar$baseEggsStanza[k] <- SpawnBio
#       simpar$Rbase[k]          <- NageS[lastMoJuv[i] + 1] * WageS[lastMoJuv[i] + 1]            
#       #leaving out code to rescale recruitment if recruitment is seasonal
#       #calculates RscaleSplit[i] to give same annual avg as equal monthly rec
#       simpar$RscaleSplit[k]     <- 1.0
#       
#       simpar$recruits[k]       <- recruits
#       simpar$RzeroS[k]         <- recruits
#       simpar$stanzaGGJuv[k]    <- juvQB * juvBio / JuvBasePred
#       simpar$stanzaGGAdu[k]    <- aduQB * aduBio / AduBasePred
#       
#       simpar$state_NN[rJuvNum]       <- sum(NageS * JuvAges)
#       simpar$stanzaPred[rJuvNum]     <- JuvBasePred
#       simpar$stanzaBasePred[rJuvNum] <- JuvBasePred
#       
#       simpar$state_NN[rAduNum]       <- sum(NageS * AduAges)
#       simpar$stanzaPred[rAduNum]     <- AduBasePred
#       simpar$stanzaBasePred[rAduNum] <- AduBasePred
#       
#       simpar$WageS     [1:(lastMoAdu[i] + 1), k] <- WageS
#       simpar$WWa       [1:(lastMoAdu[i] + 1), k] <- WWa
#       simpar$NageS     [1:(lastMoAdu[i] + 1), k] <- NageS
#       simpar$SplitAlpha[1:(lastMoAdu[i] + 1), k] <- SplitAlpha
#       
#     }  #END OF JUVENILE LOOP
#   }
#   #Reset QQ in line with new comsumtion
#   #simpar$QQ <- simpar$QQ * Qmult[simpar$PreyTo]
#   return(simpar)
# }
# 
# #'Adds final variables to be used in ecosim.run
# #'
# #'Final step of ecosim.init that adds output matrices to the Rpath.sim object
# #'that will be used during ecosim.run.
# #'
# #'@family Rpath functions
# #'
# #'@param simpar List of ecosim parameters created by ecosim.rates and modified by ecosim.stanzas
# #'@param years Integer value to set maximum number of years for the simulation.
# #'@param burn_years
# #'@param coupled
# #'@param crash_yr
# #'@param max_months Maximum number of months that a stanza can survive.
# #'@param min_rec_factor Used to set minimum logistic maturity to 0.001.  Default value is ln(1 / 0.001 - 1).
# #'
# #'@return Returns an Rpath.sim object that can be supplied to the ecosim.run function.
# #'@export
# ecosim.pack <- function(simpar, years, burn_years = -1, coupled = 1, crash_yr = -1){
#   simpar$YEARS      <- years  
#   simpar$BURN_YEARS <- burn_years
#   simpar$COUPLED    <- coupled
#   simpar$CRASH_YEAR <- crash_yr
#   
#   #derivlist
#   simpar$TotGain        <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$TotLoss        <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$LossPropToB    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$LossPropToQ    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$DerivT         <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$dyt            <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$biomeq         <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$FoodGain       <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$DetritalGain   <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$FoodLoss       <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$UnAssimLoss    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$ActiveRespLoss <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$MzeroLoss      <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$DetritalLoss   <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$TotDetOut      <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$FishingLoss    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$FishingGain    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$FishingThru    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$preyYY         <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$predYY         <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$PredSuite      <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$HandleSuite    <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$TerminalF      <- rep(0.0, simpar$NUM_GROUPS + 1)
#   
#   #targlist     
#   simpar$TARGET_BIO <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$TARGET_F   <- rep(0.0, simpar$NUM_GROUPS + 1)
#   simpar$ALPHA      <- rep(0.0, simpar$NUM_GROUPS + 1)
#   
#   mforcemat             <- as.data.frame(matrix(1.0, years * 12 + 1, simpar$NUM_GROUPS + 1))
#   names(mforcemat)      <- simpar$spname
#   simpar$force_byprey   <- mforcemat
#   simpar$force_bymort   <- mforcemat
#   simpar$force_byrecs   <- mforcemat
#   simpar$force_bysearch <- mforcemat
#   
#   yforcemat           <- as.data.frame(matrix(0.0, years + 1, simpar$NUM_GROUPS + 1))
#   names(yforcemat)    <- simpar$spname
#   simpar$FORCED_FRATE <- yforcemat
#   simpar$FORCED_CATCH <- yforcemat
#   
#   omat           <- as.data.frame(matrix(0.0, years * 12 + 1, simpar$NUM_GROUPS + 1))
#   names(omat)    <- simpar$spname
#   simpar$out_BB  <- omat 
#   simpar$out_CC  <- omat
#   simpar$out_SSB <- omat              
#   simpar$out_rec <- omat
#   
#   return(simpar)
# }
# 
# 
# ##------------------------------------------------------------------------------ 
# #'Ecosim modual of Rpath
# #'
# #'Runs the ecosim modual.
# #'
# #'@family Rpath functions
# #'
# #'@param simpar Rpath.sim object containing ecosim parameters (Generated by ecosim.init).
# #'@param BYY, EYY Integer values for the beginning/end year.
# #'@param init_run Optional flag to declare if this is an inital run of the ecosim model.
# #'
# #'@return Returns an Rpath.sim object.
# #'@useDynLib Rpath
# #'@importFrom Rcpp sourceCpp
# #'@export
# ecosim.run <- function(simpar, BYY = 0, EYY = 0, init_run = 0){
#   
#   if((EYY <= 0) | (EYY > simpar$YEARS)) EYY <- simpar$YEARS
#   if(BYY < 0)                           BYY <- 0
#   if(BYY >= simpar$YEARS)               BYY <- simpar$YEARS - 1
#   if(EYY <= BYY)                        EYY <- BYY + 1
#   
#   #Assign initial run flag
#   simpar$init_run <- init_run
#   
#   #Rcpp doesn't handle data frames well so need to convert to matrices
#   simpar$force_byprey   <- as.matrix(simpar$force_byprey)
#   simpar$force_bymort   <- as.matrix(simpar$force_bymort)
#   simpar$force_byrecs   <- as.matrix(simpar$force_byrecs)
#   simpar$force_bysearch <- as.matrix(simpar$force_bysearch)
#   simpar$FORCED_FRATE   <- as.matrix(simpar$FORCED_FRATE)
#   simpar$FORCED_CATCH   <- as.matrix(simpar$FORCED_CATCH)
#   simpar$out_BB         <- as.matrix(simpar$out_BB)
#   simpar$out_CC         <- as.matrix(simpar$out_CC)
#   simpar$out_SSB        <- as.matrix(simpar$out_SSB)       
#   simpar$out_rec        <- as.matrix(simpar$out_rec)
#   
#   #Run C code
#   Adams_Basforth(simpar, BYY, EYY)
#   
#   #Convert outputs to data frames
#   simpar$force_byprey   <- as.data.frame(simpar$force_byprey)
#   simpar$force_bymort   <- as.data.frame(simpar$force_bymort)
#   simpar$force_byrecs   <- as.data.frame(simpar$force_byrecs)
#   simpar$force_bysearch <- as.data.frame(simpar$force_bysearch)
#   simpar$FORCED_FRATE   <- as.data.frame(simpar$FORCED_FRATE)
#   simpar$FORCED_CATCH   <- as.data.frame(simpar$FORCED_CATCH)
#   simpar$out_BB         <- as.data.frame(simpar$out_BB)
#   simpar$out_CC         <- as.data.frame(simpar$out_CC)
#   simpar$out_SSB        <- as.data.frame(simpar$out_SSB)       
#   simpar$out_rec        <- as.data.frame(simpar$out_rec)
#   
#   names(simpar$force_byprey)   <- simpar$spname
#   names(simpar$force_bymort)   <- simpar$spname
#   names(simpar$force_byrecs)   <- simpar$spname
#   names(simpar$force_bysearch) <- simpar$spname
#   names(simpar$FORCED_FRATE)   <- simpar$spname
#   names(simpar$FORCED_CATCH)   <- simpar$spname
#   names(simpar$out_BB)         <- simpar$spname
#   names(simpar$out_CC)         <- simpar$spname
#   names(simpar$out_SSB)        <- simpar$spname
#   names(simpar$out_rec)        <- simpar$spname