 ###############################################################################
# Ecosense Routines (Aydin et al. 2005) as
# implemented for the R version of ecosim (Rsim for Rpath) 
#
# Current version based on Whitehouse et al (submitted in 2019)
#
################################################################################

################################################################################
#' Ecosense function for rpath (rsim.sense)
#'
#' Kerim Aydin 31 December 2019
#'@family Rpath functions
#'
#' This function generates random parameters around a given rsim scenario
#' object.  
#'
#' It has been tested to give the same results as Whitehouse and Aydin 
#' (Submitted) when supplied with a straight-from-ecopath rsim secnario.
#' and the same random seed - though some minor bug corrections have been
#' found in the previous version (previous version is rsim.sense.path)
#'
#'@param Rsim.scenario Rsim scenario object used as center of distributions
#'   (base model) during random parameter generation.
#'@param Rpath.params Rpath parameter object (unbalanced model) used for
#'    data pedigree input.
#'@param Vvary length-2 vector with (lower,upper) bounds of vulnerability
#'   generation in log-space - 1, Vvary = log(X-1) so Walters et al. 1997
#'   range of (1..inf) centered on 2 becomes (-inf,+inf) centered on 0.
#'@param Dvary length-2 vector with (lower,upper) bounds of handling time
#'   generation in log-space - 1, scaled as Vvary, above.
#'
#'@return Returns an Rsim.params object that can be substituted for the params
#'    in an Rsim.scenario object.
#'@useDynLib Rpath
#' @export
rsim.sense <- function(Rsim.scenario, Rpath.params, 
                         Vvary=c(0,0), Dvary=c(0,0)){
 
  # A "read only" version of the input params, for reference  
    orig.params  <- Rsim.scenario$params

  # The output params, initialized as copy of input params
    sense.params <- Rsim.scenario$params

  # handy numbers
    nliving <- orig.params$NUM_LIVING
    ndead   <- orig.params$NUM_DEAD
    ngears  <- orig.params$NUM_GEARS
    ngroups <- orig.params$NUM_GROUPS

  # HERE SHOULD BE THE ONLY USE OF THE Rpath.params object - 
  # do all extractions from Rpath.params in this section 
  # Pedigree vectors, including zeroes for gear groups.
    BBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,2])),rep(0,ngears))
    PBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,3])),rep(0,ngears))
    QBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,4])),rep(0,ngears))
    DCVAR   <-   as.numeric(unlist(Rpath.params$pedigree[,5]))  
    #DCVAR  <- c(as.numeric(unlist(Rpath.params$pedigree[,5])),rep(0,Rpath$NUM_GEARS))
    TYPE    <- Rpath.params$model$Type 

  # KYA 12/30/19 when drawing from a scenario, the Outside is already added.
  # But we have to strip it off to align with pedigree using IND of 2..(ngroups+1) 
  # (generating an extra number would break the random seed setting alignment)
    IND <- 2:(ngroups+1)

  # Biomass (the most straightforward).  VARiance is broken out so runif can
  # easily be replaced by other standard distributions.
    ranBB <- orig.params$B_BaseRef[IND]  * (1 + BBVAR*runif(ngroups,-1.0,1.0))

  # PB, QB, Mzero, and the dependent respiration fractions
    ranPB <- orig.params$PBopt[IND]      * (1 + PBVAR*runif(ngroups,-1.0,1.0))                        
    ranQB <- orig.params$FtimeQBOpt[IND] * (1 + QBVAR*runif(ngroups,-1.0,1.0))

    # To keep Mzero scaled to EE, we back-calculate EE then use that
    # EE = 1 - Mzero/PB - Note that this introduces tiny numerical
    # differences (1e-16 order) compared to doing it from Rpath EE directly.
    ranM0 <- ranPB * (orig.params$MzeroMort[IND] / orig.params$PBopt[IND]) *
                                           (1 + PBVAR*runif(ngroups,-1.0,1.0))

    # Active respiration is proportion of CONSUMPTION that goes to "heat"
    # Passive respiration/ VonB adjustment is left out here
    # Switched ranQB>0 test to TYPE test as QB for prim prods isn't 0 in rsim
# TODO Fix negative respiration rates?
    ranActive <- 1.0 - (ranPB/ranQB) - orig.params$UnassimRespFrac[IND]

  # Now copy these the output scenario, prepending the "Outside" values
  # and fixing values for different group types.
    sense.params$B_BaseRef      <- c(1.0, ranBB)
    sense.params$PBopt          <- c(1.0, ranPB)      
    sense.params$FtimeQBOpt     <- c(1.0, ifelse(TYPE==1, ranPB, ranQB))     
    sense.params$MzeroMort      <- c(0.0, ifelse(TYPE==3, 0, ranM0))
    sense.params$ActiveRespFrac <- c(0.0, ifelse(TYPE<1 , ranActive, 0))

  # Version of new QB saved without "Outside" group, used later 
    QBOpt <- sense.params$FtimeQBOpt[IND]  

# TODO IMPORTANT - confirm that NoIntegrate gets switched appropriately for stanzas?
  # No Integrate for A-B, fixing 2*steps_yr*steps_m as 24
    sense.params$NoIntegrate <- 
      ifelse(sense.params$MzeroMort*sense.params$B_BaseRef > 24, 0, sense.params$spnum)
  
  # KYA 12/31/19 - We don't need to recaculated predprey links, just Q(links)
  # So the section from 'primTo <- ifelse' to  'sense.params$PreyTo   <- c(primTo'
  # has been deleted.
  
  ##### This is where we add uncertainty to diet  #####
  # Original version created this diet comp vector from Rpath$DC and
  # recalculated all the links (i.e. PreyFrom, PreyTo).  
  # KYA 12/31/19 new version calculates DC from QQ vector instead.

    # Old Version got DC vector from Rpath like this:
    #    DCvector <- c(rep(0.0, sum(Rpath$type==1)), Rpath$DC[Rpath$DC>0])
    #
    # For new version, first strip the Outside link off the predprey link lists
    PPIND <- 2:(orig.params$NumPredPreyLinks+1) 
    sp.PreyTo <- orig.params$PreyTo[PPIND]
    Qvector   <-     orig.params$QQ[PPIND]
    # Then convert to DC proportions based on Q summed by PreyTo groups
    # We don't actually need to replace PP groups with 0 (using type) but 
    # setting those to 0's means rgamma's aren't generated so the set.seed 
    # alignment is off for comparisons with the old method unless we do this.
    # IMPORTANT:  The tapply method words for PreyTo (e.g. predators) because
    # there is no predator group 0 - 0 breaks array lookups.  It doesn't work 
    # for Prey - for prey, convert to character as in the preyprey loop below.
    QDCvector <- ifelse(TYPE[sp.PreyTo]==1,0, 
                 Qvector/(tapply(Qvector, sp.PreyTo, "sum")[sp.PreyTo]) )
    DCvector  <- QDCvector 
    # KYA 12/31/19 the remainer of this section (drawing from a gamma that
    # based on DC that is then normalized) is unchanged from previous version.
    # Diet comp pedigree
    DCpedigree <- DCVAR[sp.PreyTo]
    ## Random diet comp
# TODO: make EPSILON lower??  (1e-16 might work?)
    EPSILON <- 1*10^-8
    betascale <- 1.0
    DCbeta <- betascale * DCpedigree * DCpedigree
    alpha <- ifelse(DCbeta>0,DCvector/DCbeta, DCvector)
    DClinks <- rgamma(length(DCvector), shape=alpha, rate=DCbeta)
    # We do need to check cases were Beta=0 (implying no variance)   
    DClinks1 <- ifelse(DCbeta<=0, DCvector, DClinks)
    DClinks2 <- ifelse(DClinks1 < EPSILON, 2 * EPSILON, DClinks1)
    # DClinks2 prevents random diet comps from becoming too low, effectively
    # equal to zero. Zeros in DClinks will produce NaN's in sense.params$QQ, and
    # others, ultimately preventing ecosim.
    DCtot <- tapply(DClinks2, sp.PreyTo, "sum")    
    # Normalized diet comp
    DCnorm <- ifelse(TYPE[sp.PreyTo]==1, 1.0, DClinks2/DCtot[sp.PreyTo])
    # The "if" part of DCnorm is so the DC of phytoplankton (type==1) won't equal zero
    DCQB  <- QBOpt[sp.PreyTo]
    DCBB  <- ranBB[sp.PreyTo]  
    ranQQ <- DCnorm * DCQB * DCBB             	
   
  # Sarah used the following formula to vary vulnerability in Gaichas et al. (2012)
  # That paper states that "vulnerability" (also known as X*predprey) has an
  # effective range from 1.01 to 91 in EwE.
  # KYA 12/30/19: changed to allow specification of V range and D range
  # by user (in log space deviations from original) 
  # Note the old version of this set DD to 1001, not 1000 (minor bug)  
    NPP    <-  length(PPIND)
    log.v  <-  log(orig.params$VV[PPIND] - 1) 
    log.d  <-  log(orig.params$DD[PPIND] - 1) 
    ranVV  <-  1 + exp(runif(NPP, log.v+Vvary[1], log.v+Vvary[2]))
    ranDD  <-  1 + exp(runif(NPP, log.d+Dvary[1], log.d+Dvary[2]))

  # Scramble combined prey pools
  # KYA 12/30/19 PredTotWeight and PreyTotWeight are no longer exported in the
  # creation of an Rsim scenario - so made local only.  Then changed to use
  # tapply instead of a sum looper.
  # IMPORTANT:  since 0 appears as an index for tapply (for prey=0), we need
  # to convert the pred and prey indices to characters and use a character
  # lookup to account for the 0.
    namedBB <- c(1.0,ranBB); names(namedBB)<-0:ngroups
    py   <- as.character(orig.params$PreyFrom[PPIND]) # + 1.0
    pd   <- as.character(orig.params$PreyTo[PPIND])   # + 1.0
    prey.BB <- namedBB[py]
    pred.BB <- namedBB[pd]
    VV   <- as.numeric(ranVV*ranQQ / prey.BB)
    AA   <- as.numeric((2.0 * ranQQ * VV) / (VV*pred.BB*prey.BB - ranQQ*pred.BB))
    ranPredPred <- as.numeric(AA * pred.BB)
    ranPreyPrey <- as.numeric(AA * prey.BB)
    ranPredPredWeight <- as.numeric(ranPredPred/(tapply(ranPredPred, py, "sum")[py]))
    ranPreyPreyWeight <- as.numeric(ranPreyPrey/(tapply(ranPreyPrey, pd, "sum")[pd]))

    #sense.params$PreyFrom       <- c(0, sense.params$PreyFrom)
    #sense.params$PreyTo         <- c(0, sense.params$PreyTo)
    sense.params$QQ             <- c(0, ranQQ)
    sense.params$DD             <- c(1000, ranDD) # using default of 1000 for first group 
    sense.params$VV             <- c(2, ranVV)    # using default of 2 for first group
    sense.params$PredPredWeight <- c(0, ranPredPredWeight)
    sense.params$PreyPreyWeight <- c(0, ranPreyPreyWeight)

  # Fisheries Catch
  # The Whitehouse et al. version didn't actualy vary catch, the only thing it 
  # does it renomalize FishQ.  Therefore only Q calculation lines are needed.  
  # Since Q is now calculated slightly differently (not from rpath) there are 
  # slight (10^-18) differences from the previous version.
    CIND = orig.params$FishFrom + 1
    sense.params$FishQ <- orig.params$FishQ * 
                          orig.params$B_BaseRef[CIND]/sense.params$B_BaseRef[CIND] 
  
#TODO add catch uncertainty

  class(sense.params) <- 'Rsim.params'
  return(sense.params)   
  
}


################################################################################
#' Ecosense function for rpath (rsim.sense.path)
#'
#' 5 November 2019
#'@family Rpath functions
#'
#' rsim.sense.path is depreciated and included for testing only, the preferred
#' function is rsim.sense().
#' Whitehouse and Aydin (Submitted) Assessing the sensitivity of three Alaska marine
#' food webs to perturbations: an example of Ecosim simulations using Rpath
#'
#' This function generates random parameters around the Ecopath baseline,
#' not the scenario.  
#'
#'@return Returns an Rsim.scenario object that can be supplied to the rsim.run function.
#'@useDynLib Rpath
#' @export
rsim.sense.path <- function(Rsim.scenario, Rpath, Rpath.params,
                            steps_yr = 12, steps_m = 1){
   
  sense.params <- Rsim.scenario$params

  nliving <- Rpath$NUM_LIVING
  ndead   <- Rpath$NUM_DEAD
  
  # Set-up pedigree vectors, including zeroes for gear groups.
  BBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,2])),rep(0,Rpath$NUM_GEARS))
  PBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,3])),rep(0,Rpath$NUM_GEARS))
  QBVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,4])),rep(0,Rpath$NUM_GEARS))
  DCVAR	<- c(as.numeric(unlist(Rpath.params$pedigree[,5])),rep(0,Rpath$NUM_GEARS))
  

  # Biomass
  ranBB <- Rpath$BB * (1 + BBVAR * runif(Rpath$NUM_GROUPS,-1.0,1.0))
  sense.params$B_BaseRef <- c(1.0, ranBB)
  # PB
  ranPB <- Rpath$PB * (1 + PBVAR * runif(Rpath$NUM_GROUPS,-1.0,1.0))
  # QB
  ranQB <- Rpath$QB * (1 + QBVAR * runif(Rpath$NUM_GROUPS,-1.0,1.0))

  # Mzero
  sense.params$MzeroMort <- c(0.0, ranPB * (1.0 - Rpath$EE)  *
                          (1 + PBVAR * runif(Rpath$NUM_GROUPS,-1.0,1.0)))
    # The zero at the start of the M0 vector is for "outside"
    # same for the other vectors

  #Active respiration is proportion of CONSUMPTION that goes to "heat"
  #Passive respiration/ VonB adjustment is left out here
  sense.params$ActiveRespFrac <-  c(0.0, ifelse(ranQB > 0, 
                                          1.0 - (ranPB / ranQB) - Rpath$GS, 
                                          0.0))

  #####	rQBPrimProd is a added here because we need a non-zero QB for
  #####	primary production groups. QB's of zero produce NaN's in ecosim
  #####	and just bring the whole thing to a screaching hault. rQBPrimProd
  #####	sets the QB of PP groups equal to their PB.
  rQBPrimProd	<- ifelse(Rpath$type==1,ranPB,ranQB)
  sense.params$FtimeQBOpt <-   c(1.0, rQBPrimProd)
  #sense.params$FtimeQBOpt <-   c(1.0, ranQB)
  sense.params$PBopt      <-   c(1.0, ranPB)           

  #No Integrate
  sense.params$NoIntegrate <- ifelse(sense.params$MzeroMort * sense.params$B_BaseRef > 
                                 2 * steps_yr * steps_m, 0, sense.params$spnum)
  
  #primary production links
  primTo   <- ifelse(Rpath$type > 0 & Rpath$type <= 1, 
                     1:length(ranPB),
                     0)
  primFrom <- rep(0, length(Rpath$PB))
  primQ    <- ranPB * ranBB 

  # Change production to consusmption for mixotrophs
  mixotrophs <- which(Rpath$type > 0 & Rpath$type < 1)
  primQ[mixotrophs] <- primQ[mixotrophs] / Rpath$GE[mixotrophs] * 
          Rpath$type[mixotrophs]

  #Predator/prey links
  preyfrom  <- row(Rpath$DC)
  preyto    <- col(Rpath$DC)	
  predpreyQ <- Rpath$DC[1:(nliving + ndead + 1), ] * 
    t(matrix(rQBPrimProd[1:Rpath$NUM_LIVING] * ranBB[1:Rpath$NUM_LIVING],
            nliving, nliving + ndead + 1))
  
  #combined
  sense.params$PreyFrom <- c(primFrom[primTo > 0], preyfrom [predpreyQ > 0])
  # Changed import prey number to 0
  sense.params$PreyFrom[which(sense.params$PreyFrom == nrow(Rpath$DC))] <- 0
  sense.params$PreyTo   <- c(primTo  [primTo > 0], preyto   [predpreyQ > 0])
  
  ##### This is where we add uncertainty to diet  #####
  # Diet comp vector
  DCvector <- c(rep(0.0, sum(Rpath$type==1)), Rpath$DC[Rpath$DC>0])
  # Diet comp pedigree
  DCped <- as.numeric(unlist(Rpath.params$pedigree[,5]))
  DCpedigree <- DCped[sense.params$PreyTo]
  ## Random diet comp
  EPSILON <- 1*10^-8
  betascale <- 1.0
  DCbeta <- betascale * DCpedigree * DCpedigree
  alpha <- DCvector/DCbeta
  DClinks <- rgamma(length(DCvector), shape=alpha, rate=DCbeta)
  DClinks2 <- ifelse(DClinks < EPSILON, 2 * EPSILON, DClinks)
  # DClinks2 prevents random diet comps from becoming too low, effectively
  # equal to zero. Zeros in DClinks will produce NaN's in sense.params$QQ, and
  # others, ultimately preventing ecosim.
  DCtot <- tapply(DClinks2, sense.params$PreyTo, "sum")    
  # Normalized diet comp
  DCnorm <- ifelse(Rpath$type[sense.params$PreyTo]==1, 1.0, DClinks2/DCtot[sense.params$PreyTo])
  # The "if" part of DCnorm is so the DC of phytoplankton (type==1) won't equal zero
  DCQB <- rQBPrimProd[sense.params$PreyTo]
  DCBB <- ranBB[sense.params$PreyTo]  
  sense.params$QQ <- DCnorm * DCQB * DCBB             	
   
  numpredprey <- length(sense.params$QQ)

  # Sarah used the following formula to vary vulnerability in Gaichas et al. (2012)
  # That paper states that "vulnerability" (also known as X*predprey) has an
  # effective range from 1.01 to 91 in EwE.
  sense.params$VV	<-	1 + exp(9 * (runif(length(sense.params$QQ))-0.5))
  sense.params$DD	<-	1000 + exp(0 * (runif(length(sense.params$QQ))-0.5))

  # Scramble combined prey pools
  Btmp <- sense.params$B_BaseRef
  py   <- sense.params$PreyFrom + 1.0
  pd   <- sense.params$PreyTo + 1.0
  VV   <- sense.params$VV * sense.params$QQ / Btmp[py]
  AA   <- (2.0 * sense.params$QQ * VV) / (VV * Btmp[pd] * Btmp[py] - sense.params$QQ * Btmp[pd])
  sense.params$PredPredWeight <- AA * Btmp[pd] 
  sense.params$PreyPreyWeight <- AA * Btmp[py] 
  
  sense.params$PredTotWeight <- rep(0, length(sense.params$B_BaseRef))
  sense.params$PreyTotWeight <- rep(0, length(sense.params$B_BaseRef))
  
  for(links in 1:numpredprey){
    sense.params$PredTotWeight[py[links]] <- sense.params$PredTotWeight[py[links]] + sense.params$PredPredWeight[links]
    sense.params$PreyTotWeight[pd[links]] <- sense.params$PreyTotWeight[pd[links]] + sense.params$PreyPreyWeight[links]    
  }  

  sense.params$PredPredWeight <- sense.params$PredPredWeight/sense.params$PredTotWeight[py]
  sense.params$PreyPreyWeight <- sense.params$PreyPreyWeight/sense.params$PreyTotWeight[pd]

  sense.params$PreyFrom       <- c(0, sense.params$PreyFrom)
  sense.params$PreyTo         <- c(0, sense.params$PreyTo)
  sense.params$QQ             <- c(0, sense.params$QQ)
  sense.params$DD             <- c(0, sense.params$DD)
  sense.params$VV             <- c(0, sense.params$VV) 
  sense.params$PredPredWeight <- c(0, sense.params$PredPredWeight)
  sense.params$PreyPreyWeight <- c(0, sense.params$PreyPreyWeight)


  #catchlinks
  fishfrom    <- row(as.matrix(Rpath$Catch))
  fishthrough <- col(as.matrix(Rpath$Catch)) + (nliving + ndead)
  fishcatch   <- Rpath$Catch
  fishto      <- fishfrom * 0
  
  if(sum(fishcatch) > 0){
    sense.params$FishFrom    <- fishfrom   [fishcatch > 0]
    sense.params$FishThrough <- fishthrough[fishcatch > 0]
    sense.params$FishQ       <- fishcatch  [fishcatch > 0] / sense.params$B_BaseRef[sense.params$FishFrom + 1]  
    sense.params$FishTo      <- fishto     [fishcatch > 0]
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
      sense.params$FishFrom    <- c(sense.params$FishFrom,    fishfrom   [fishcatch > 0])
      sense.params$FishThrough <- c(sense.params$FishThrough, fishthrough[fishcatch > 0])
      ffrom <- fishfrom[fishcatch > 0]
      sense.params$FishQ       <- c(sense.params$FishQ,  fishcatch[fishcatch > 0] / sense.params$B_BaseRef[ffrom + 1])  
      sense.params$FishTo      <- c(sense.params$FishTo, fishto   [fishcatch > 0])
    }
  } 

  sense.params$FishFrom        <- c(0, sense.params$FishFrom)
  sense.params$FishThrough     <- c(0, sense.params$FishThrough)
  sense.params$FishQ           <- c(0, sense.params$FishQ)  
  sense.params$FishTo          <- c(0, sense.params$FishTo)   

 
  class(sense.params) <- 'Rsim.params'
  return(sense.params)   
  
}

################################################################################ 
