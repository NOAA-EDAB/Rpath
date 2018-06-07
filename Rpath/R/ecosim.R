 ################################################################################
# R version of ecosim 
# originally developed by Kerim Aydin
# modified by Sean Lucey
#
#####################################################################################


#'Rsim modual of Rpath
#'
#'Prepares a balanced Rpath model and creates a scenario consisting of 5 list objects:
#'params, start_state, forcing, fishing, and stanzas.
#'
#'@family Rpath functions
#'
#'@param Rpath R object containing a balanced Rpath model.
#'@param Rpath.params R object containing the Rpath parameters.  This is generated
#'  either by the create.rpath.params or read.rpath.params functions.
#'@param years A vector of each year of the simulation.
#'
#'@return Returns an Rsim.scenario object that can be supplied to the rsim.run function.
#'@import data.table
#'@useDynLib Rpath
#'@importFrom Rcpp sourceCpp
#'@export
rsim.scenario <- function(Rpath, Rpath.params, years = 1:100){
  # KYA 11/1/17 modifying so years can take a vector of actual years
  # for row labels (for fitting and general printing out)
  if (length(years)<2){stop("Years needs to be a vector of numeric year labels.")}
  
  params      <- rsim.params(Rpath)
  start_state <- rsim.state(params)
  forcing     <- rsim.forcing(params, years)  
  fishing     <- rsim.fishing(params, years)
  stanzas     <- rsim.stanzas(Rpath.params, start_state, params)
  
  #Set NoIntegrate Flags
  ieco <- as.vector(stanzas$EcopathCode[which(!is.na(stanzas$EcopathCode))])
  params$NoIntegrate[ieco + 1] <- -1 * ieco 
  
  rsim = list(params      = params, 
              start_state = start_state,
              forcing     = forcing,
              fishing     = fishing,
              stanzas     = stanzas)
  
  class(rsim) <- 'Rsim.scenario'
  attr(rsim, 'eco.name') <- attr(Rpath, 'eco.name')
  attr(rsim, 'Start year') <- years[1]
  
  return(rsim)   
}
#####################################################################################
# Runs Ecosim
#'@export
rsim.run <- function(Rpath.scenario, method = 'RK4', years = 1:100){
  scene <- copy(Rpath.scenario) 

  # Figure out starting and ending years for run
  # KYA 4/23/18 single year (e.g. 1971:1971) reduces to a scalar so take out length trap
    #if (length(years)<2){stop("Years should be a vector of year labels")}
    scene.years <- row.names(Rpath.scenario$fishing$FRATE)
    syear <- which(as.character(head(years,1))==scene.years)
    eyear <- which(as.character(tail(years,1))==scene.years)
    if (eyear<syear)     {stop("End year cannot be less than start year.")}
    if (length(syear)!=1){stop("Starting year not found in scenario (or more than once).")}
    if (length(eyear)!=1){stop("Ending year not found in scenario (or more than once).")}
  
  # KYA adds run date and some random salt to ensure uniquieness  
    scene$rundate <- paste(Sys.time(),":salt:",runif(1))    
    
  if(method == 'RK4'){
    rout <- rk4_run(scene$params,  scene$start_state, 
                    scene$forcing, scene$fishing,
                    scene$stanzas, syear, eyear)
  }
  if(method == 'AB'){
    #Run initial derivative
    derv <- deriv_vector(scene$params, scene$start_state, 
                         scene$forcing, scene$fishing, 
                         scene$stanzas, syear, 0, 0)
    #KYA added for first step bump correction 9/20/17
     
    #Run Adams Bashforth Alogrithm
    rout <- Adams_run(scene$params,  scene$start_state, 
                      scene$forcing, scene$fishing,
                      scene$stanzas, syear, eyear, derv)
  }
  # Nicely Name output vectors
  sps <- scene$params$spname[1:(1+scene$params$NUM_BIO)]
  colnames(rout$out_BB)    <- sps
  colnames(rout$out_CC)    <- sps
  colnames(rout$annual_CC) <- sps
  colnames(rout$annual_BB) <- sps
  colnames(rout$annual_QB) <- sps
  colnames(rout$annual_Qlink)<-1:(length(rout$annual_Qlink[1,]))

  # drop the last row (should be always 0; negative index is entry to drop)
    #lastone <- length(rout$annual_CC[,1])
    #rout$annual_CC <-    rout$annual_CC[-lastone,]
    #rout$annual_BB <-    rout$annual_BB[-lastone,]
    #rout$annual_QB <-    rout$annual_QB[-lastone,]
    #rout$annual_Qlink <- rout$annual_Qlink[-lastone,]  
  
  #if(!is.null(Rpath.scenario$fitting$years)){
  #  ylist <- seq(scene$fitting$years[1],length.out=length(rout$annual_CC[,1]))
    # put years in row names
     ys <- min(as.numeric(rownames(scene$fishing$CATCH)))
     ylist <- seq(ys,length.out=length(rout$annual_CC[,1]))
      rownames(rout$annual_CC) <-    ylist #Rpath.scenario$fitting$years
      rownames(rout$annual_BB) <-    ylist #Rpath.scenario$fitting$years
      rownames(rout$annual_QB) <-    ylist #Rpath.scenario$fitting$years
      rownames(rout$annual_Qlink) <- ylist #Rpath.scenario$fitting$years   
  #}
  
  rout$pred <- scene$params$spname[scene$params$PreyTo+1] 
  rout$prey <- scene$params$spname[scene$params$PreyFrom+1] 
  rout$Gear_CC_sp   <- scene$params$spname[scene$params$FishFrom+1] 
  rout$Gear_CC_gear <- scene$params$spname[scene$params$FishThrough+1] 
  rout$Gear_CC_disp <- ifelse(scene$params$FishTo == 0, 'Landing', 'Discard')
  
  rout$start_state       <- scene$start_state
  rout$params$NUM_LIVING <- scene$params$NUM_LIVING
  rout$params$NUM_DEAD   <- scene$params$NUM_DEAD
  rout$params$NUM_GEARS  <- scene$params$NUM_GEARS
  rout$params$spname     <- scene$params$spname
  
  class(rout) <- 'Rsim.output'
  attr(rout, 'eco.name') <- attr(scene, 'eco.name')
  
  return(rout)
}
 
################################################################################
#'@export
rsim.fishing <- function(params, years){
  # Yearly index defaulting to to 0.0, for fishing forcing list
  nyrs <- length(years)
  #if (length(years)>1){nyrs <- length(years)} else {nyrs <- years}
  
  MF <- (matrix(1.0, nyrs * 12, params$NUM_GEARS + 1))
  YF <- (matrix(0.0, nyrs, params$NUM_BIO   + 1))
  #GF <- (matrix(1.0, nyrs, params$NUM_GEARS + 1))
  #if (length(years)>1){
  rownames(YF) <- years
  year.m <- c()
  for(iyear in 1:nyrs){
    iyear.m <- paste0(years[iyear], '.', 1:12)
    year.m  <- c(year.m, iyear.m)
  }
  rownames(MF) <- year.m
  #}
  colnames(YF) <- params$spname[1:(params$NUM_BIO+1)]
  colnames(MF) <- c("Outside",params$spname[(params$NUM_BIO+2):(params$NUM_GROUPS+1)])
  fishing <- list(EFFORT = MF,
                  FRATE  = YF,
                  CATCH  = YF)   
  
  class(fishing) <- "Rsim.fishing"
  return (fishing)
}

#####################################################################################
#'@export
rsim.forcing <- function(params, years){
# Monthly index defaulting to to 1.0, for environmental forcing list
  
  nyrs <- length(years)
  MF <- (matrix(1.0, nyrs * 12, params$NUM_GROUPS + 1))
  year.m <- c()
  for(iyear in 1:nyrs){
    iyear.m <- paste0(years[iyear], '.', 1:12)
    year.m  <- c(year.m, iyear.m)
  }
  rownames(MF) <- year.m
  colnames(MF) <- params$spname
  forcing <- list(byprey    = MF, 
                  bymort    = MF, 
                  byrecs    = MF, 
                  bysearch  = MF,
                  bymigrate = MF * 0,
                  bybio     = MF * -1)
  
  class(forcing) <- "Rsim.forcing"
  return (forcing)
}

#####################################################################################
#'@export
rsim.state <- function(params){
  state  <- list(BB    = params$B_BaseRef, 
                 NN    = rep(0, params$NUM_GROUPS + 1),
                 Ftime = rep(1, length(params$B_BaseRef)))
  class(state) <- "Rsim.state"
  return(state)
}

#####################################################################################
#'@export
rsim.diet <- function(rout,predator){
   pmat <-  rout$annual_Qlink[,rout$pred==predator]
   colnames(pmat) <- rout$prey[rout$pred==predator]
   return(pmat)
}  

#####################################################################################
#'@export
rsim.mort <- function(rout,prey){
  pmat <-  rout$annual_Qlink[,rout$prey==prey]
  colnames(pmat) <- rout$pred[rout$prey==prey]
  return(pmat)
}  
#####################################################################################
#'@export
rsim.deriv <- function(Rpath.scenario, year=0, month=0, tstep=0){
  scene <- copy(Rpath.scenario)
  rout <- deriv_vector(scene$params,  scene$start_state, 
                       scene$forcing, scene$fishing,
                       scene$stanzas, year, month, tstep)
  
  rtab <- data.frame(scene$params$spname,rout$DerivT,rout$TotGain,rout$TotLoss, 
                    rout$FoodGain, rout$DetritalGain, rout$FishingGain,      
                    rout$UnAssimLoss,rout$ActiveRespLoss,
                    rout$FoodLoss,rout$MzeroLoss,rout$FishingLoss,rout$DetritalLoss)
  colnames(rtab)<-c("Species","DerivT","TotGain","TotLoss","FoodGain","DetritalGain","FishingGain",
                    "UnAssimLoss","ActiveRespLoss","FoodLoss","MzeroLoss","FishingLoss","DetritalLoss")
  return(rtab)
}

#####################################################################################
#'Initial set up for Ecosim modual of Rpath
#'
#'Converts the outputs from ecopath into rates for use in ecosim.
#'
#'@family Rpath functions
#'
#'@param Rpath Rpath object containing a balanced model.
#'@param mscramble WILL REMOVE
#'@param mhandle WILL REMOVE
#'@param preyswitch WILL REMOVE - Adjust with adjust.scenario
#'@param scrambleselfwt Value of 1 indicates no overlap while 0 indicates complete overlap.
#'@param handleselfwt Value of 1 indicates no overlap while 0 indicates complete overlap.
#'@param steps_yr Number of time steps per year.
#'@param steps_m Number of time steps per month.
#'
#'@return Returns an Rpath.sim object that can be supplied to the ecosim.run function.
#'@export
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
  simpar$spnum      <- 0:length(Rpath$BB) 
  
  #Energetics for Living and Dead Groups
  #Reference biomass for calculating YY
  simpar$B_BaseRef <- c(1.0, Rpath$BB) 
  #Mzero proportional to (1-EE)
  simpar$MzeroMort <- c(0.0, Rpath$PB * (1.0 - Rpath$EE)) 
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
  simpar$fish_Effort <- ifelse(simpar$spnum <= nliving + ndead,
                               0.0,
                               1.0) 
  
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
  primQ    <- Rpath$PB * Rpath$BB
  #Change production to consumption for mixotrophs
  mixotrophs <- which(Rpath$type > 0 & Rpath$type < 1)
  primQ[mixotrophs] <- primQ[mixotrophs] / Rpath$GE[mixotrophs] * 
    Rpath$type[mixotrophs] 
  
  #Import links
  # importTo <- ifelse(Rpath$DC[nrow(Rpath$DC), ] > 0,
  #                    1:ncol(Rpath$DC),
  #                    0)
  # importFrom <- rep(0, ncol(Rpath$DC))
  # importQ <- Rpath$DC[nrow(Rpath$DC), ] * Rpath$QB[1:nliving] * Rpath$BB[1:nliving]
  
  #Predator/prey links
  preyfrom  <- row(Rpath$DC)
  preyto    <- col(Rpath$DC)	
  predpreyQ <- Rpath$DC[1:(nliving + ndead + 1), ] * 
    t(matrix(Rpath$QB[1:nliving] * Rpath$BB[1:nliving],
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
  fishthrough <- col(as.matrix(Rpath$Catch)) + (nliving + ndead)
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
  simpar$FishFrom        <- c(0, simpar$FishFrom)
  simpar$FishThrough     <- c(0, simpar$FishThrough)
  simpar$FishQ           <- c(0, simpar$FishQ)  
  simpar$FishTo          <- c(0, simpar$FishTo)   
  
# SET DETRITAL FLOW
  detfrac <- Rpath$DetFate[1:(nliving + ndead), ]
  detfrom <- row(as.matrix(detfrac))
  detto   <- col(as.matrix(detfrac)) + nliving
  
  detout <- 1 - rowSums(as.matrix(Rpath$DetFate[1:(nliving + ndead), ]))
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
  simpar$RK4_STEPS  <- 4.0 
  
  class(simpar) <- "Rsim.params"
  return(simpar)
}
 
 #####################################################################################
 #'@export
 rsim.stanzas <- function(Rpath.params, state, params){
   juvfile <- Rpath.params$stanzas
   if(Rpath.params$stanzas$NStanzaGroups > 0){
     #Set up multistanza parameters to pass to C
     rstan <- list()
     rstan$Nsplit      <- juvfile$NStanzaGroups
     #Need leading zeros (+1 col/row) to make indexing in C++ easier
     rstan$Nstanzas    <- c(0, juvfile$stgroups$nstanzas)
     rstan$EcopathCode <- matrix(NA, rstan$Nsplit + 1, max(rstan$Nstanzas) + 1)
     rstan$Age1        <- matrix(NA, rstan$Nsplit + 1, max(rstan$Nstanzas) + 1)
     rstan$Age2        <- matrix(NA, rstan$Nsplit + 1, max(rstan$Nstanzas) + 1)
     rstan$WageS       <- matrix(NA, max(juvfile$stindiv$Last) + 1, rstan$Nsplit + 1)
     rstan$NageS       <- matrix(NA, max(juvfile$stindiv$Last) + 1, rstan$Nsplit + 1)
     rstan$WWa         <- matrix(NA, max(juvfile$stindiv$Last) + 1, rstan$Nsplit + 1)
     
     sPred <- rep(0, params$NUM_GROUPS + 1) #rstan$stanzaPred  <- rep(0, params$NUM_GROUPS + 1)
     
     for(isp in 1:rstan$Nsplit){
       for(ist in 1:rstan$Nstanzas[isp + 1]){
         rstan$EcopathCode[isp + 1, ist + 1] <- juvfile$stindiv[StGroupNum == isp &
                                                          StanzaNum == ist, GroupNum]
         rstan$Age1[isp + 1, ist + 1] <- juvfile$stindiv[StGroupNum == isp & 
                                                   StanzaNum == ist, First]
         rstan$Age2[isp + 1, ist + 1] <- juvfile$stindiv[StGroupNum == isp & 
                                                   StanzaNum == ist, Last]
       }
       rstan$WageS[1:nrow(juvfile$StGroup[[isp]]), isp + 1] <- juvfile$StGroup[[isp]]$WageS
       rstan$NageS[1:nrow(juvfile$StGroup[[isp]]), isp + 1] <- juvfile$StGroup[[isp]]$NageS
       rstan$WWa[  1:nrow(juvfile$StGroup[[isp]]), isp + 1] <- juvfile$StGroup[[isp]]$WWa
     }
     
     #Maturity
     MIN_REC_FACTOR     <- 6.906754779  #// This is ln(1/0.001 - 1) used to set min. logistic matuirty to 0.001
     rstan$Wmat     <- c(0, juvfile$stgroup$Wmat)
     #rstan$Wmat001  <- c(0, juvfile$stgroup$Wmat001)
     #rstan$Wmat50   <- c(0, juvfile$stgroup$Wmat50)
     #rstan$Amat001  <- c(0, juvfile$stgroup$Amat001)
     #rstan$Amat50   <- c(0, juvfile$stgroup$Amat50)
     #rstan$WmatSpread <- -(rstan$Wmat001 - rstan$Wmat50)/MIN_REC_FACTOR
     #rstan$AmatSpread <- -(rstan$Amat001 - rstan$Amat50)/MIN_REC_FACTOR
     rstan$RecPower <- c(0, juvfile$stgroup$RecPower)
     rstan$recruits <- c(0, juvfile$stgroup$r)
     rstan$vBGFd    <- c(0, juvfile$stgroup$VBGF_d)
     rstan$RzeroS   <- rstan$recruits
     
     #Energy required to grow a unit in weight(scaled to Winf = 1)
     rstan$vBM <- c(0, (1 - 3 * juvfile$stgroups$VBGF_Ksp / 12))
     
     # Spawning Biomass and Eggs
     
     eggs <- c(0) #rstan$baseEggsStanza <- c(0)
     for(isp in 1:rstan$Nsplit){
       #id which weight at age is higher than Wmat
       #rstan$baseEggsStanza[isp + 1]
       # KYA 5/3/18 switchout
       #if ((is.na(rstan$Wmat[isp+1]))|(rstan$Wmat[isp+1]<0)){
       #     eggs[isp + 1] <- 
       #    juvfile$StGroup[[isp]][(WageS>rstan$Wmat001[isp+1])&(age>rstan$Amat001[isp+1]),
       #    sum(WageS*NageS / (1. + exp( -((WageS - rstan$Wmat50[isp+1]) / rstan$WmatSpread[isp+1])   
       #                                  -((age   - rstan$Amat50[isp+1]) / rstan$AmatSpread[isp+1]) ))) ]
       #     rstan$Wmat[isp+1] <- -1.0
       #}
       #else{
          eggs[isp + 1] <- juvfile$StGroup[[isp]][WageS > rstan$Wmat[isp+1], 
                                                   sum(NageS * (WageS - rstan$Wmat[isp+1]))]
       #}
     }
     #rstan$EggsStanza <- rstan$baseEggsStanza
     
     #initialize splitalpha growth coefficients using pred information and
     rstan$SplitAlpha <- matrix(NA, max(juvfile$stindiv$Last) + 1, rstan$Nsplit + 1)
     for(isp in 1:rstan$Nsplit){
       for(ist in 1:rstan$Nstanzas[isp + 1]){
         ieco  <- rstan$EcopathCode[isp + 1, ist + 1]
         first <- rstan$Age1[isp + 1, ist + 1]
         last  <- rstan$Age2[isp + 1, ist + 1]
         pred  <- sum(juvfile$StGroup[[isp]][age %in% first:last, NageS * WWa])
         StartEatenBy <- juvfile$stindiv[StGroupNum == isp & StanzaNum == ist, Cons]
  
         SplitAlpha <- (juvfile$StGroup[[isp]][, shift(WageS, type = 'lead')] - 
           rstan$vBM[isp + 1] * juvfile$StGroup[[isp]][, WageS]) * pred / StartEatenBy
         rstan$SplitAlpha[(first + 1):(last + 1), isp + 1] <- SplitAlpha[(first + 1):
                                                                       (last + 1)]
         sPred[ieco + 1] <- pred
       }
       #Carry over final split alpha to plus group
       rstan$SplitAlpha[rstan$Age2[isp + 1, rstan$Nstanzas[isp + 1] + 1] + 1, isp + 1] <- 
         rstan$SplitAlpha[rstan$Age2[isp + 1, rstan$Nstanza[isp + 1]], isp + 1]
     }
     
     #Misc parameters for C
     #KYA Spawn X is Beverton-Holt.  To turn off set to 10000. 2 is half saturation.
     #1.00001 or so is minimum
     rstan$SpawnX         <- c(0, rep(10000, rstan$Nsplit))
     rstan$SpawnEnergy    <- c(0, rep(1, rstan$Nsplit))
     eggs1<-eggs+1; rstan$baseEggsStanza <- eggs1-1
     eggs1<-eggs+1; rstan$EggsStanza     <- eggs1-1
     eggs1<-eggs+1; rstan$SpawnBio       <- eggs1-1
     eggs1<-eggs+1; rstan$baseSpawnBio   <- eggs1-1
     rstan$RscaleSplit    <- c(0, rep(1, rstan$Nsplit))
     sPred1<-sPred+1; rstan$stanzaPred     <- sPred1-1
     sPred1<-sPred+1; rstan$stanzaBasePred <- sPred1-1
     
    # SplitSetPred(rstan, state)
   }
   
   if(juvfile$NStanzaGroups == 0){
     rstan <- list()
     rstan$Nsplit      <- juvfile$NStanzaGroups
     #Need leading zeros (+1 col/row) to make indexing in C++ easier
     rstan$Nstanzas       <- c(0, 0)
     rstan$EcopathCode    <- matrix(rep(0, 4), 2, 2)
     rstan$Age1           <- matrix(rep(0, 4), 2, 2)
     rstan$Age2           <- matrix(rep(0, 4), 2, 2)
     rstan$WageS          <- matrix(rep(0, 4), 2, 2)
     rstan$NageS          <- matrix(rep(0, 4), 2, 2)
     rstan$WWa            <- matrix(rep(0, 4), 2, 2)
     rstan$stanzaPred     <- c(0, 0)
     rstan$Wmat           <- c(0, 0)
     #rstan$Wmat001        <- c(0, 0)
     #rstan$Wmat50         <- c(0, 0)
     #rstan$Amat001        <- c(0, 0)
     #rstan$Amat50         <- c(0, 0)
     #rstan$WmatSpread <- -(rstan$Wmat001 - rstan$Wmat50)/MIN_REC_FACTOR
     #rstan$AmatSpread <- -(rstan$Amat001 - rstan$Amat50)/MIN_REC_FACTOR
     rstan$RecPower       <- c(0, 0)
     rstan$recruits       <- c(0, 0)
     rstan$VBGFd          <- c(0, 0)
     rstan$RzeroS         <- c(0, 0)
     rstan$vBM            <- c(0, 0)
     rstan$baseEggsStanza <- c(0, 0)
     rstan$EggsStanza     <- c(0, 0)
     rstan$SplitAlpha     <- matrix(rep(0, 4), 2, 2)
     rstan$SpawnX         <- c(0, 0)
     rstan$SpawnEnergy    <- c(0, 0)
     rstan$SpawnBio       <- c(0, 0)
     rstan$baseSpawnBio   <- c(0, 0)
     rstan$RscaleSplit    <- c(0, 0)
     rstan$stanzaBasePred <- c(0, 0)
     
   }
     
   
   return(rstan)
 }

 
