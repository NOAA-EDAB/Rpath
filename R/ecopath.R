#'Ecopath module of Rpath
#'
#'Performs initial mass balance using a model parameter file and diet
#'matrix file.
#'
#'@family Rpath functions
#'
#'@param Rpath.params R object containing the Rpath parameters.  This is generated
#'  either by the create.rpath.params or read.rpath.params functions.
#'@param eco.name Optional name of the ecosystem which becomes an attribute of
#'    rpath object.
#'@param eco.area Optional area of the ecosystem which becomes an attribute of the
#'    rpath object.
#'
#'@return Returns an Rpath object that can be supplied to the rsim.scenario function.
#'@import data.table
#'@export
rpath <- function(Rpath.params, eco.name = NA, eco.area = 1){
  #Need to define variables to eliminate check() note about no visible binding
  Type <- Group <- DetInput <- ProdCons <- PB <- QB <- noB <- noEE <- alive <- NULL
  BEE <- Biomass <- Q <- BioAcc <- BioQB <- diag.a <- EEa <- B <- M0 <- NULL
  QBloss <- Unassim <- Ex <- NULL
  
  # Model Parameters - Basic parameters, detritus fate, catch, discards in that order
  model <- copy(Rpath.params$model)
  
  #Diet Parameters - diet matrix, predators as columns, prey as rows - include
  #producers as predators even though they do not consume any groups
  diet <- copy(Rpath.params$diet)
  
  #Check that all columns of model are numeric and not logical
  if(length(which(sapply(model, class) == 'logical')) > 0){
    logic.col <- which(sapply(model, class) == 'logical')
    for(i in 1:length(logic.col)){
      set(model, j = logic.col[i], value = as.numeric(model[[logic.col[i]]]))
    }
  }
    
  #Remove first column if names (factor or character)
  if(sapply(diet, class)[1] == 'factor')    diet[, 1 := NULL]
  if(sapply(diet, class)[1] == 'character') diet[, 1 := NULL]

  #Adjust diet comp of mixotrophs
  mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  mix.Q <- 1 - model[mixotrophs, Type]
  for(i in seq_along(mixotrophs)){
    new.dc <- diet[, mixotrophs[i], with = F] * mix.Q[i]
    diet[, mixotrophs[i] := new.dc]
  }
  
  #Convert NAs to zero in diet matrix
  diet[is.na(diet)] <- 0
  
  # Get number of groups, living, dead, and gear
  ngroups <- nrow(model)
  nliving <- nrow(model[Type <  2, ])
  ndead   <- nrow(model[Type == 2, ])
  ngear   <- nrow(model[Type == 3, ])

  nodetrdiet <- diet[1:nliving, ]
  model[is.na(DetInput), DetInput := 0]

  # fill in GE(PQ), QB, or PB from other inputs
  GE   <- ifelse(is.na(model[, ProdCons]), model[, PB / QB],       model[, ProdCons])
  QB.1 <- ifelse(is.na(model[, QB]),       model[, PB / GE],       model[, QB])
  PB.1 <- ifelse(is.na(model[, PB]),       model[, ProdCons * QB], model[, PB])
  model[, QB := QB.1]
  model[, PB := PB.1]
  
  # define landings, discards, necessary sums
  landmat     <- model[, (10 + ndead + 1):(10 + ndead + ngear), with = F]
  discardmat  <- model[, (10 + ndead + 1 + ngear):(10 + ndead + (2 * ngear)), with = F]
  totcatchmat <- landmat + discardmat
    
  # KYA 1/16/14 Need if statement here because rowSums fail if only one 
  # fishery (catch is vector instead of matrix)     ##FIX PROPAGATION HERE
  if (is.data.frame(totcatchmat)){
    totcatch <- rowSums(totcatchmat)
    landings <- rowSums(landmat)    
    discards <- rowSums(discardmat)  
    gearland <- colSums(landmat,   na.rm = T)
    geardisc <- colSums(discardmat, na.rm = T)
  }else{
    totcatch <- totcatchmat
    landings <- landmat    
    discards <- discardmat 
    gearland <- sum(landmat,    na.rm = T)
    geardisc <- sum(discardmat, na.rm = T)                     
  }   
  
  geartot <- gearland + geardisc
  model[, landings := landings]
  model[, discards := discards]
  model[, totcatch := totcatch]

  # flag missing pars and subset for estimation
  model[, noB   := 0]
  model[, noEE  := 0]
  model[, alive := 0]
  model[, BEE   := 0]
  model[is.na(Biomass), noB   := 1]
  model[is.na(EE),      noEE  := 1]
  model[Type < 2,       alive := 1]
  model[noB == 0 & noEE == 0, BEE := 1]
  
  # define detritus fate matrix
  detfate <- model[, (10 + 1):(10 + ndead), with = F]

  # set up and solve the system of equations for living group B or EE
  living  <- model[alive == 1, ]
  
  #Set up right hand side b
  living[, Ex := totcatch + BioAcc]
  living[, BioQB := Biomass * QB]
  cons  <- as.matrix(nodetrdiet) * living$BioQB[col(as.matrix(nodetrdiet))]
  living[, b := Ex + rowSums(cons, na.rm = T)] 
  
  #Set up A matrix
  living[noEE == 1, diag.a := Biomass * PB]
  living[noEE == 0, diag.a := PB * EE]
  
  #Special case where B and EE are known then need to solve for BA
  #living[BEE == 1, b := b - (Biomass * PB * EE)]
  #living[BEE  == 1, diag.a := 0] #Need to work on this solution
  
  A       <- matrix(0, nliving, nliving)
  diag(A) <- living[, diag.a]
  QBDC    <- as.matrix(nodetrdiet) * living$QB[col(as.matrix(nodetrdiet))]
  dimnames(QBDC) <- list(NULL, NULL)
  QBDC[is.na(QBDC)] <- 0
  #Flip noB flag for known B and EE
  #living[BEE == 1, noB := 1]
  QBDCa <- as.matrix(QBDC) * living$noB[col(as.matrix(QBDC))]
  A     <- A - QBDCa 
  #Switch flag back
  #living[BEE == 1, noB := 0]
   
  # Generalized inverse does the actual solving
  #Invert A and multiple by b to get x (unknowns)
  x <- MASS::ginv(A, tol = .Machine$double.eps) %*% living[, b]
  
  #Assign unknown values
  living[, EEa := x * noEE]
  living[is.na(EE), EE := EEa]
  
  living[, B := x * noB]
  living[is.na(Biomass), Biomass := B]

  # detritus EE calcs
  living[, M0 := PB * (1 - EE)]
  living[, QBloss := QB]
  living[is.na(QBloss), QBloss := 0]
  loss <- c((living[, M0] * living[, Biomass]) + 
              (living[, Biomass] * living[, QBloss] * living[, Unassim]),
            model[Type ==2, DetInput], 
            geardisc)
  detinputs  <- colSums(loss * detfate)
  detdiet    <- diet[(nliving + 1):(nliving + ndead), ]
  BQB        <- living[, Biomass * QB]
  detcons    <- as.matrix(detdiet) * BQB[col(as.matrix(detdiet))]
  detoutputs <- rowSums(detcons, na.rm = T)
  EE         <- c(living[, EE], as.vector(detoutputs / detinputs))

  # added by kya
  # if a detritus biomass is put into the spreadsheet, use that and 
  # calculate PB.  If no biomass, but a PB, use that pb with inflow to 
  # calculate biomass.  If neither, use default PB=0.5, Bio = inflow/PB  
  # This is done because Ecosim requires a detrital biomass.
  Default_Detrital_PB <- 0.5 
  inDetPB <- model[(nliving + 1):(nliving + ndead), PB] 
  inDetB  <- model[(nliving + 1):(nliving + ndead), Biomass]
  DetPB   <- ifelse(is.na(inDetPB), Default_Detrital_PB, inDetPB)
  DetB    <- ifelse(is.na(inDetB), detinputs / DetPB, inDetB)
  DetPB   <- detinputs / DetB

  # Trophic Level calcs
  b             <- rep(1, ngroups)
  TLcoeff       <- matrix(0, ngroups, ngroups)
  diag(TLcoeff) <- rep(1, ngroups)
  gearcons      <- as.matrix(totcatchmat) / geartot[col(as.matrix(totcatchmat))]
  dimnames(gearcons) <- list(NULL, NULL)
  gearcons[is.na(gearcons)] <- 0
  dietplus <- as.matrix(diet)
  dimnames(dietplus) <- list(NULL, NULL)

  #Adjust for mixotrophs (partial primary producers) - #Moved this code up so that
  #it also impacted the EE calculation
  # mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  # mix.Q <- 1 - model[mixotrophs, Type]
  # for(i in seq_along(mixotrophs)){
  #   dietplus[, mixotrophs[i]] <- dietplus[, mixotrophs[i]] * mix.Q[i]
  # }
  #Adjust for diet import (Consumption outside model)
  import <- which(dietplus[nrow(diet), ] > 0)
  for(i in seq_along(import)){
    import.denom <- 1 - dietplus[nrow(diet), import[i]]
    dietplus[, import[i]] <- dietplus[, import[i]] / import.denom
  }
  dietplus <- dietplus[1:(nliving + ndead), ]
  dietplus <- rbind(dietplus, matrix(0, ngear, nliving))
  dietplus <- cbind(dietplus, matrix(0, ngroups, ndead), gearcons)
  TLcoeffA <- TLcoeff - dietplus
  TL       <- solve(t(TLcoeffA), b)     

  #kya changed these following four lines for detritus, and removing NAs
  #to match header file format (replacing NAs with 0.0s)
  Bplus  <- c(living[, Biomass], DetB, rep(0.0, ngear))
  
  PBplus <- model[, PB] 
  PBplus[(nliving + 1):(nliving + ndead)] <- DetPB
  PBplus[is.na(PBplus)] <- 0.0
  
  EEplus <- c(EE, rep(0.0, ngear))
  
  QBplus <- model[, QB]
  QBplus[is.na(QBplus)] <- 0.0
  
  GE[is.na(GE)] <- 0.0
  
  RemPlus <- model[, totcatch]
  RemPlus[is.na(RemPlus)] <- 0.0
  
  balanced <- list(Group    = model[, Group], 
                   TL       = TL, 
                   Biomass  = Bplus, 
                   PB       = PBplus, 
                   QB       = QBplus, 
                   EE       = EEplus, 
                   GE       = GE, 
                   Removals = RemPlus)

  M0plus  <- c(living[, M0], as.vector(detoutputs / detinputs))
  gearF   <- as.matrix(totcatchmat) / living[, Biomass][row(as.matrix(totcatchmat))]
  #newcons <- as.matrix(nodetrdiet)  * living[, BQB][col(as.matrix(nodetrdiet))]
  newcons <- as.matrix(nodetrdiet)  * BQB[col(as.matrix(nodetrdiet))]
  predM   <- as.matrix(newcons) / living[, Biomass][row(as.matrix(newcons))]
  predM   <- rbind(predM, detcons)
  morts   <- list(Group = model[Type < 3, Group], 
                  PB    = model[Type < 3, PB], 
                  M0    = M0plus, 
                  F     = gearF[1:(nliving + ndead), ], 
                  M2    = predM)
  
  # convert from levels to characters
  gnames <- as.character(balanced$Group)
  
  # cleanup before sending to sim -- C code wants 0 as missing value, not NA
  balanced$Biomass[is.na(balanced$Biomass)] <- 0
  balanced$PB[is.na(balanced$PB)]     <- 0
  balanced$QB[is.na(balanced$QB)]     <- 0
  balanced$EE[is.na(balanced$EE)]     <- 0
  balanced$GE[is.na(balanced$GE)]     <- 0
  model$BioAcc[is.na(model$BioAcc)]   <- 0
  model$Unassim[is.na(model$Unassim)] <- 0
  dietm                               <- as.matrix(diet)
  dimnames(dietm)                     <- list(c(gnames[1:(nliving+ndead)],"Import"), gnames[1:nliving])
  dietm[is.na(dietm)]                 <- 0
  landmatm                            <- as.matrix(landmat)
  dimnames(landmatm)                  <- list(gnames, gnames[(ngroups-ngear+1):ngroups])
  landmatm[is.na(landmatm)]           <- 0
  discardmatm                         <- as.matrix(discardmat)
  dimnames(discardmatm)               <- list(gnames, gnames[(ngroups-ngear+1):ngroups])
  discardmatm[is.na(discardmatm)]     <- 0
  detfatem                            <- as.matrix(detfate)
  dimnames(detfatem)                  <- list(gnames, gnames[(nliving+1):(nliving+ndead)])
  detfatem[is.na(detfatem)]           <- 0

  # KYA April 2020 - added names for output list
    out.Group   <- gnames;           names(out.Group) <- gnames
    out.type    <- model[, Type];    names(out.type) <- gnames
    out.TL      <- TL;               names(out.TL) <- gnames
    out.Biomass <- balanced$Biomass; names(out.Biomass) <- gnames
    out.PB      <- balanced$PB;      names(out.PB) <- gnames
    out.QB      <- balanced$QB;      names(out.QB) <- gnames
    out.EE      <- balanced$EE;      names(out.EE) <- gnames
    out.BA      <- model[, BioAcc];  names(out.BA) <- gnames
    out.Unassim <- model[, Unassim]; names(out.Unassim) <- gnames
    out.GE      <- balanced$GE;      names(out.GE) <- gnames    
    
  # list structure for sim inputs
  path.model <- list(NUM_GROUPS = ngroups,
                     NUM_LIVING = nliving,
                     NUM_DEAD   = ndead,
                     NUM_GEARS  = ngear,
                     Group      = out.Group,
                     type       = out.type,
                     TL         = out.TL,
                     Biomass    = out.Biomass,
                     PB         = out.PB,
                     QB         = out.QB,
                     EE         = out.EE,
                     BA         = out.BA,
                     Unassim    = out.Unassim,
                     GE         = out.GE,
                     DC         = dietm,
                     DetFate    = detfatem,
                     Landings   = landmatm,
                     Discards   = discardmatm)      

#Define class of output
class(path.model) <- 'Rpath'
attr(path.model, 'eco.name') <- eco.name
attr(path.model, 'eco.area') <- eco.area

return(path.model)
}


#'Calculate biomass and consumption for multistanza groups
#'
#'Uses the leading stanza to calculate the biomass and consumption of other stanzas
#'necessary to support the leading stanza.
#'
#'@family Rpath functions
#'
#'@inheritParams rpath
#'
#'@return Calculates and adds biomass and consumption for trailing stanza groups.  
#'  Also adds weight at age and number at age for multi-staza groups.
#'  
#'@import data.table
#'@export 
rpath.stanzas <- function(Rpath.params){
  #Need to define variables to eliminate check() note about no visible binding
  StGroupNum <- First <- StanzaNum <- VBGF_d <- VBGF_Ksp <- Last <- GroupNum <- NULL
  WageS <- age <- QageS <- Survive <- Z <- survive_L <- bs.num <- qs.num <- Leading <- Oldest <- NULL
  Group <- Biomass <- R <- NageS <- bs.denom <- bs <- qs.denom <- qs <- Cons <- NULL
  QB <- BAB <- Ex <- NULL

  #Determine the total number of groups with multistanzas
  Nsplit     <- Rpath.params$stanza$NStanzaGroups
  groupfile  <- Rpath.params$stanza$stgroups
  stanzafile <- Rpath.params$stanza$stindiv
  
  # Add Oldest column so age groups will be properly calculated
  stanzafile[Last == max(Last), Oldest := T, by = StGroupNum]
  stanzafile[is.na(Oldest), Oldest := F]
  
  #Need to add vector of stanza number
  lastmonth <- rep(NA,Nsplit)
  for(isp in 1:Nsplit){
    #Put the stanzas in order for each split species
    stnum <- order(stanzafile[StGroupNum == isp, First])
    
    stanzafile[StGroupNum == isp, StanzaNum := stnum]

    
    #Calculate the last month for the final ("leading") stanza
    #KYA Aug 2021:
    # Formerly used fraction of Winf, but that didn't work for species
    # with rapid growth but low mortality (e.g. marine mammals).
    # So instead, calculate biomass out for a very long time and
    # taking 0.99999 of cumulative biomass as a cutoff.
    
    #this selects all of the stanza lines, then picks the last one
    #(maybe data table has a better way...)
    stmax <- max(stanzafile[StGroupNum == isp, StanzaNum])
    st <- stanzafile[StGroupNum == isp & StanzaNum==stmax,]

    gp <- groupfile[isp,]
    #Max age class in months should be one less than a multiple of 12
    #(trying 5999 - probably overkill but for safety)
    AGE <- st$First:5999
    mz <- (st$Z + gp$BAB)/12
    k  <- gp$VBGF_Ksp
    d  <- gp$VBGF_d
    NN <- shift(cumprod(rep(exp(-1*mz),length(AGE))),1,1.0)
    BB <- NN * (1 - exp(-k * (1 - d) * (AGE))) ^ (1 / (1 - d))
    BBcum <- cumsum(BB)/sum(BB)
    #Age at which 0.99999 of cumulative leading stanza biomass is represented,
    #rounded to nearest higher multiple of 12 (-1 since index starts at 0)
    lastmonth[isp] <- ceiling(AGE[min(which(BBcum>0.99999))]/12) * 12 - 1
  }
  
  #Save the maximum month vector in the table
  groupfile[, last := lastmonth]
  
  for(isp in 1:Nsplit){
    nstanzas <- groupfile[
      StGroupNum == isp, nstanzas]
    stanzafile[StGroupNum == isp & Oldest, Last := lastmonth[isp]]
    
    #Grab ecopath group codes
    group.codes <- stanzafile[StGroupNum == isp, GroupNum]
    
    #Grab index for first and last months for stanzas
    first   <- stanzafile[StGroupNum == isp, First]
    second  <- stanzafile[StGroupNum == isp, Last]

    #Calculate weight and consumption at age    
    StGroup <- data.table(age = stanzafile[StGroupNum == isp & StanzaNum == 1, First]:
                            lastmonth[isp])
    #Calculate monthly generalized k: (Ksp * 3) / 12
    k <- (groupfile[StGroupNum == isp, VBGF_Ksp] * 3) / 12
    d <-  groupfile[StGroupNum == isp, VBGF_d]
    StGroup[, WageS := (1 - exp(-k * (1 - d) * (age))) ^ (1 / (1 - d))]
    StGroup[, QageS   := WageS ^ d]
    
    #Calculate the relative number of animals at age a
    #Vector of survival rates from 1 stanza to the next
    
    #Unwind the by-stanza mortality rates into by-month survival rates
    #by looping through the stanzas
    survive_L <- rep(NA,length(StGroup$age))
    for(ist in 1:nstanzas){
      #Convert Z to a monthly Z
      month.z <- (stanzafile[StGroupNum == isp & StanzaNum == ist, Z] + 
                   groupfile[StGroupNum == isp, BAB]) / 12
      survive_L[which(StGroup$age %in% first[ist]:second[ist])] <- exp(-1*month.z)
    } 
    
    #Shift survival rates forward one - first survival is 1.0
    survive_L <- shift(survive_L,1,1.0)
    #Use cumulative product to get overall survival to each month
    StGroup[, Survive := cumprod(survive_L)]
      
    StGroup[, B := Survive * WageS]
    StGroup[, Q := Survive * QageS]
    for(ist in 1:nstanzas){      
      #Numerator for the relative biomass/consumption calculations
      b.num <- StGroup[age %in% first[ist]:second[ist], sum(B)]
      q.num <- StGroup[age %in% first[ist]:second[ist], sum(Q)]
      
      stanzafile[StGroupNum == isp & StanzaNum == ist, bs.num := b.num]
      stanzafile[StGroupNum == isp & StanzaNum == ist, qs.num := q.num]
    }
      
    #Scale numbers up to total recruits
    BaseStanza <- stanzafile[StGroupNum == isp & Leading == T, ]
    BioPerEgg <- StGroup[age %in% BaseStanza[, First]:BaseStanza[, Last], sum(B)]
    recruits <- Rpath.params$model[Group == BaseStanza[, Group], Biomass] / BioPerEgg
    #Save recruits
    groupfile[StGroupNum == isp, R := recruits]

    #Numbers at age S
    StGroup[, NageS := Survive * recruits]
    
    #Calculate relative biomass
    stanzafile[StGroupNum == isp, bs.denom := sum(bs.num)]
    stanzafile[StGroupNum == isp, bs := bs.num / bs.denom]
    
    #Calculate relative consumption
    stanzafile[StGroupNum == isp, qs.denom := sum(qs.num)]
    stanzafile[StGroupNum == isp, qs := qs.num / qs.denom]
    
    #Use leading group to calculate other biomasses
    stanzafile[StGroupNum == isp & Leading == T, 
               Biomass := Rpath.params$model[Group == BaseStanza[, Group], Biomass]]
    B <- stanzafile[StGroupNum == isp & Leading == T, Biomass / bs]
    stanzafile[StGroupNum == isp, Biomass := bs * B]
    
    #Use leading group to calculate other consumption
    stanzafile[StGroupNum == isp & Leading == T, 
               Cons := Rpath.params$model[Group == BaseStanza[, Group], QB] *
                 Rpath.params$model[Group == BaseStanza[, Group], Biomass]]
    Q <- stanzafile[StGroupNum == isp & Leading == T, Cons / qs]
    stanzafile[StGroupNum == isp, Cons := qs * Q]
    stanzafile[, QB := Cons / Biomass]
    
  Rpath.params$stanzas$StGroup[[isp]] <- StGroup
  
  }
  
  #Drop extra columns
  stanzafile[, c('bs.num', 'bs.denom', 'bs', 'qs.num', 'qs.denom', 'qs') := NULL]
  
  #Push biomass to modfile
  for(i in 1:nrow(stanzafile)){
    Rpath.params$model[Group == stanzafile[i, Group], Biomass := stanzafile[i, Biomass]]
  }
  
  #Push consumption to modfile
  for(i in 1:nrow(stanzafile)){
    Rpath.params$model[Group == stanzafile[i, Group], QB := stanzafile[i, QB]]
  }
  
  return(Rpath.params)
} 
