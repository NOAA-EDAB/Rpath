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
  #mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  #mix.Q <- 1 - model[mixotrophs, Type]
  # for(i in seq_along(mixotrophs)){
  #   new.dc <- diet[, mixotrophs[i], with = F] * mix.Q[i]
  #   diet[, mixotrophs[i] := new.dc, with = F]
  # }
  
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
  
  # define catch, discards, necessary sums
  catchmat    <- model[, (10 + ndead + 1):(10 + ndead + ngear), with = F]
  discardmat  <- model[, (10 + ndead + 1 + ngear):(10 + ndead + (2 * ngear)), with = F]
  totcatchmat <- catchmat + discardmat
    
  # KYA 1/16/14 Need if statement here because rowSums fail if only one 
  # fishery (catch is vector instead of matrix)     ##FIX PROPAGATION HERE
  if (is.data.frame(totcatchmat)){
    totcatch  <- rowSums(totcatchmat)
    catch     <- rowSums(catchmat)    
    discards  <- rowSums(discardmat)  
    gearcatch <- colSums(catchmat,   na.rm = T)
    geardisc  <- colSums(discardmat, na.rm = T)
  }else{
    totcatch  <- totcatchmat
    catch     <- catchmat    
    discards  <- discardmat 
    gearcatch <- sum(catchmat,   na.rm = T)
    geardisc  <- sum(discardmat, na.rm = T)                     
  }   
  
  geartot <- gearcatch + geardisc
  model[, catch    := catch]
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
  
  #Set up right hand side Q
  living[, Q := totcatch + BioAcc]
  living[, BioQB := Biomass * QB]
  cons  <- as.matrix(nodetrdiet) * living$BioQB[col(as.matrix(nodetrdiet))]
  living[, Q := Q + rowSums(cons, na.rm = T)] 
  living[BEE == 1, Q := Q - (Biomass * PB * EE)]
  
  #Set up A matrix
  living[noEE == 1, diag.a := Biomass * PB]
  living[noEE == 0, diag.a := PB * EE]
  living[noEE == 0 & noB == 0, diag.a := 0]
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
  pars <- MASS::ginv(A, tol = .Machine$double.eps) %*% living[, Q]
  
  living[, EEa := pars * noEE]
  living[is.na(EE), EE := EEa]
  living[, EEa := NULL]
  living[, B := pars * noB]
  living[!is.na(Biomass), B := Biomass]

  # detritus EE calcs
  living[, M0 := PB * (1 - EE)]
  living[, QBloss := QB]
  living[is.na(QBloss), QBloss := 0]
  loss <- c((living[, M0] * living[, B]) + (living[, B] * living[, QBloss] * 
                                              living[, Unassim]),
            model[Type ==2, DetInput], 
            geardisc)
  detinputs  <- colSums(loss * detfate)
  detdiet    <- diet[(nliving + 1):(nliving + ndead), ]
  BQB        <- living[, B * QB]
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
  #Adjust for mixotrophs (partial primary producers)
  mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  mix.Q <- 1 - model[mixotrophs, Type]
  for(i in seq_along(mixotrophs)){
    dietplus[, mixotrophs[i]] <- dietplus[, mixotrophs[i]] * mix.Q[i]
  }
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
  Bplus  <- c(living[, B], DetB, rep(0.0, ngear))
  
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
  gearF   <- as.matrix(totcatchmat) / living[, B][row(as.matrix(totcatchmat))]
  #newcons <- as.matrix(nodetrdiet)  * living[, BQB][col(as.matrix(nodetrdiet))]
  newcons <- as.matrix(nodetrdiet)  * BQB[col(as.matrix(nodetrdiet))]
  predM   <- as.matrix(newcons) / living[, B][row(as.matrix(newcons))]
  predM   <- rbind(predM, detcons)
  morts   <- list(Group = model[Type < 3, Group], 
                  PB    = model[Type < 3, PB], 
                  M0    = M0plus, 
                  F     = gearF[1:(nliving + ndead), ], 
                  M2    = predM)
     
  # cleanup before sending to sim -- C code wants 0 as missing value, not NA
  balanced$Biomass[is.na(balanced$Biomass)] <- 0
  balanced$PB[is.na(balanced$PB)]     <- 0
  balanced$QB[is.na(balanced$QB)]     <- 0
  balanced$EE[is.na(balanced$EE)]     <- 0
  balanced$GE[is.na(balanced$GE)]     <- 0
  model$BioAcc[is.na(model$BioAcc)]   <- 0
  model$Unassim[is.na(model$Unassim)] <- 0
  dietm                               <- as.matrix(diet)
  dimnames(dietm)                     <- list(NULL, NULL)
  dietm[is.na(dietm)]                 <- 0
  catchmatm                           <- as.matrix(catchmat)
  dimnames(catchmatm)                 <- list(NULL, NULL)
  catchmatm[is.na(catchmatm)]         <- 0
  discardmatm                         <- as.matrix(discardmat)
  dimnames(discardmatm)               <- list(NULL, NULL)
  discardmatm[is.na(discardmatm)]     <- 0
  detfatem                            <- as.matrix(detfate)
  dimnames(detfatem)                  <- list(NULL, NULL)
  detfatem[is.na(detfatem)]           <- 0

  # list structure for sim inputs
  path.model <- list(NUM_GROUPS = ngroups,     #define NUM_GROUPS 80  INCLUDES GEAR
                NUM_LIVING = nliving,          #define NUM_LIVING 60
                NUM_DEAD   = ndead,            #define NUM_DEAD 3
                NUM_GEARS  = ngear,            #define NUM_GEARS 17
                Group      = as.character(balanced$Group),
                type       = model[, Type],
                TL         = TL,
                BB         = balanced$Biomass, #float path_BB[1..NUM_GROUPS] vector
                PB         = balanced$PB,      #float path_PB[1..NUM_GROUPS] vector
                QB         = balanced$QB,      #float path_QB[1..NUM_GROUPS] vector
                EE         = balanced$EE,      #float path_EE[1..NUM_GROUPS] vector
                BA         = model[, BioAcc],  #float path_BA[1..NUM_GROUPS] vector
                GS         = model[, Unassim], #float path_GS[1..NUM_GROUPS] vector
                GE         = balanced$GE,      #float path_GS[1..NUM_GROUPS] vector
                DC         = dietm,            #float path_DC[1..NUM_GROUPS][1..NUM_GROUPS]  matrix in [prey][pred] order     NUM_LIVING?
                DetFate    = detfatem,         #float path_DetFate[1..NUM_DEAD][1..NUM_GROUPS]  matrix in [det][groups] order
                Catch      = catchmatm,        #float path_Catch[1..NUM_GEARS][1..NUM_GROUPS]  matrix
                Discards   = discardmatm)      #float path_Discards[1..NUM_GEARS][1..NUM_GROUPS] matrix


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
#'@param Rpath.params Object containing the Rpath parameters generated either by 
#'  create.rpath.params or read.rpath.params
#'
#'@return Calculates and adds biomass and consumption for trailing stanza groups.  
#'  Also adds weight at age and number at age for multi-staza groups.
#'  
#'@import data.table
#'@export 
rpath.stanzas <- function(Rpath.params){
  
  #Determine the total number of groups with multistanzas
  Nsplit     <- Rpath.params$stanza$NStanzaGroups
  groupfile  <- Rpath.params$stanza$stgroups
  stanzafile <- Rpath.params$stanza$stindiv
  
  #Need to add vector of stanza number
  for(isp in 1:Nsplit){
    stnum <- order(stanzafile[StGroupNum == isp, First])
    stanzafile[StGroupNum == isp, StanzaNum := stnum]
  }
  
  #Calculate the last month for the final stanza
  #Months to get to 99% Winf (We don't use an accumulator function like EwE)
  groupfile[, last := floor(log(1 - 0.9999^(1 - VBGF_d)) / 
                                  (-1 * (VBGF_Ksp * 3 / 12) * (1 - VBGF_d)))]
  
  for(isp in 1:Nsplit){
    nstanzas <- groupfile[StGroupNum == isp, nstanzas]
    t99 <- groupfile[StGroupNum == isp, last]
    stanzafile[StGroupNum == isp & StanzaNum == nstanzas, Last := t99]
    
    #Grab ecopath group codes
    group.codes <- stanzafile[StGroupNum == isp, GroupNum]
    
    #Grab index for first and last months for stanzas
    first   <- stanzafile[StGroupNum == isp, First]
    second  <- stanzafile[StGroupNum == isp, Last]
        
    #Calculate weight and consumption at age    
    StGroup <- data.table(age = stanzafile[StGroupNum == isp & StanzaNum == 1, First]:
                            t99)
    #Calculate monthly generalized k: (Ksp * 3) / 12
    k <- (groupfile[StGroupNum == isp, VBGF_Ksp] * 3) / 12
    d <-  groupfile[StGroupNum == isp, VBGF_d]
    StGroup[, WageS := (1 - exp(-k * (1 - d) * (age))) ^ (1 / (1 - d))]
    StGroup[, WWa   := WageS ^ d]
    
    #Calculate the relative number of animals at age a
    #Vector of survival rates from 1 stanza to the next
    StGroup[age == 0, Survive := 1]
    prev.surv <- 1
    for(ist in 1:nstanzas){
      #Convert Z to a monthly Z
      month.z <- stanzafile[StGroupNum == isp & StanzaNum == ist, Z] / 12
      StGroup[age %in% first[ist]:second[ist], surv := exp(-1*month.z)]
      
      if(first[ist] > 0){
        StGroup[age == first[ist], Survive := StGroup[age == first[ist] - 1, 
                                                    Survive] * prev.surv]
      }
      
      for(a in (first[ist] + 1):second[ist]){
        StGroup[age == a, Survive := StGroup[age == a - 1, Survive] * surv]
      }
      
      StGroup[, B := Survive * WageS]
      StGroup[, Q := Survive * WWa]
      
      #Numerator for the relative biomass/consumption calculations
      b.num <- StGroup[age %in% first[ist]:second[ist], sum(B)]
      q.num <- StGroup[age %in% first[ist]:second[ist], sum(Q)]
      
      stanzafile[StGroupNum == isp & StanzaNum == ist, bs.num := b.num]
      stanzafile[StGroupNum == isp & StanzaNum == ist, qs.num := q.num]
      
      prev.surv <- exp(-1 * month.z)
    }
      
    #Scale numbers up to total recruits
    BaseStanza <- stanzafile[StGroupNum == isp & Leading == T, ]
    BioPerEgg <- StGroup[age %in% BaseStanza[, First]:BaseStanza[, Last], sum(B)]
    recruits <- Rpath.params$model[Group == BaseStanza[, Group], Biomass] / BioPerEgg
    #Save recruits
    groupfile[StGroupNum == isp, r := recruits]

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