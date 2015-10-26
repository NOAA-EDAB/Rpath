## R version of Ecopath balance by Sarah Gaichas and Kerim Aydin
## Modified by Sean Lucey
## Version controled by git
## Function ecopathR takes as input 3 csv files and optional
## ecosystem name  


#'Ecopath modual of Rpath
#'
#'Performs initial mass balance using a model parameter file and diet
#'matrix file.
#'
#'@family Rpath functions
#'
#'@param modfile Comma deliminated model parameter file.
#'@param dietfile Comma deliminated diet matrix file.
#'@param eco.name Optional name of the ecosystem which becomes an attribute of
#'    rpath object.
#'
#'@return Returns an Rpath object that can be supplied to the ecosim.init function.
#'@import data.table
#'@export
ecopath <- function(modfile, dietfile, eco.name = NA){
  
  #Read in parameter files - either as file path or data.table object
  # Model Parameters - Basic parameters, detritus fate, catch, discards in that order
  if(is.character(modfile)){
    model <- as.data.table(read.csv(modfile))
  } else {
    model <- as.data.table(modfile)
  }
  #Diet Parameters - diet matrix, predators as columns, prey as rows - include
  #producers as predators even though they do not consume any groups
  if(is.character(dietfile)){
    diet  <- as.data.table(read.csv(dietfile))
  } else {
    diet <- as.data.table(dietfile)
  }
  
  #Check that all columns of model are numeric and not logical
  if(length(which(sapply(model, class) == 'logical')) > 0){
    logic.col <- which(sapply(model, class) == 'logical')
    for(i in 1:length(logic.col)){
      set(model, j = logic.col[i], value = as.numeric(model[[logic.col[i]]]))
    }
  }
    
  #Remove first column if names (factor or character)
  if(sapply(diet, class)[1] == 'factor')    diet <- diet[, 1 := NULL, with = F]
  if(sapply(diet, class)[1] == 'character') diet <- diet[, 1 := NULL, with = F]
  
  #Convert NAs to zero in diet matrix
  diet[is.na(diet)] <- 0
  
  # Get number of groups, living, dead, and gear
  ngroups <- nrow(model)
  nliving <- nrow(model[Type <  2, ])
  ndead   <- nrow(model[Type == 2, ])
  ngear   <- nrow(model[Type == 3, ])

  nodetrdiet <- diet[1:nliving, ]
  model[is.na(DetInput), DetInput := 0]

  # fill in GE and QB from inputs
  GE <- ifelse(is.na(model[, ProdCons]), model[, PB / QB],       model[, ProdCons])
  QB <- ifelse(is.na(model[, QB]),       model[, QB := PB / GE], model[, QB])

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
  model[is.na(Biomass), noB   := 1]
  model[is.na(EE),      noEE  := 1]
  model[Type < 2,       alive := 1]
  
  # define detritus fate matrix
  detfate <- model[, (10 + 1):(10 + ndead), with = F]

  # set up and solve the system of equations for living group B or EE
  living  <- model[alive == 1, ]
  living[, Q := totcatch + BioAcc]
  living[noEE == 1, diag.a := Biomass * PB]
  living[noEE == 0, diag.a := PB * EE]
  A       <- matrix(0, nliving, nliving)
  diag(A) <- living[, diag.a]
  QBDC    <- as.matrix(nodetrdiet) * living$QB[col(as.matrix(nodetrdiet))]
  dimnames(QBDC) <- list(NULL, NULL)
  QBDC[is.na(QBDC)] <- 0
  QBDCa <- as.matrix(QBDC) * living$noB[col(as.matrix(QBDC))]
  A     <- A - QBDCa 
  living[, BioQB := Biomass * QB]
  cons  <- as.matrix(nodetrdiet) * living$BioQB[col(as.matrix(nodetrdiet))]
  living[, Q := Q + rowSums(cons, na.rm = T)]  

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
  TL            <- rep(1, ngroups)
  TLcoeff       <- matrix(0, ngroups, ngroups)
  diag(TLcoeff) <- rep(1, ngroups)
  gearcons      <- as.matrix(totcatchmat) / geartot[col(as.matrix(totcatchmat))]
  dimnames(gearcons) <- list(NULL, NULL)
  gearcons[is.na(gearcons)] <- 0
  dietplus <- as.matrix(diet)
  dimnames(dietplus) <- list(NULL, NULL)
  dietplus <- rbind(dietplus, matrix(0, ngear, nliving))
  dietplus <- cbind(dietplus, matrix(0, ngroups, ndead), gearcons)
  TLcoeffA <- TLcoeff - dietplus
  TL       <- solve(t(TLcoeffA), TL)     

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
  newcons <- as.matrix(nodetrdiet)  * living[, BQB][col(as.matrix(nodetrdiet))]
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

return(path.model)
}


#'Calculate biomass and consumption for multistanza groups
#'
#'Uses the leading stanza to calculate the biomass and consumption of other stanzas
#'necessary to support the leading stanza.
#'
#'@family Rpath functions
#'
#'@param modfile Object containing the model parameters.
#'@param groupfile Object containing the group specific characteristics.
#'@param stanzafile Object containing the stanza specific characteristics.
#'@param output Logical value indicating whether relative abundance, weight at age,
#'and consumption at age should be output from function.  This is for use with the
#'stanzaplot function.
#'
#'@return Adds the biomass and consumption to the relative groups in the model parameter
#'object
#'@import data.table
#'@export 
rpath.stanzas <- function(modfile, groupfile, stanzafile, output = F){
  #Determine the total number of groups with multistanzas
  ngroups <- max(groupfile[, StGroupNum])
  
  #Set up output for plotting
  stanza.out <- list()
  
  for(i in 1:ngroups){  
    #Calculate last month adult
    last <- groupfile[StGroupNum == i, nstanzas]
    #Convert to generalized k from Ksp and make monthly
    k <- (groupfile[StGroupNum == i, VBGF_Ksp] * 3) / 12
    d <-  groupfile[StGroupNum == i, VBGF_d]
    #Months to get to 90% Winf
    t90 <- floor(log(1 - 0.9^(1 - d)) / (-1 * k * (1 - d)))
    stanzafile[StGroupNum == i & Stanza == last, Last := t90]
    
    #Vector of survival rates from 1 stanza to the next
    prev.surv <- 1
    
    StGroup <- data.table(age = stanzafile[StGroupNum == i & Stanza == 1, First]:
                            t90)
    StGroup[, wa := (1 - exp(-k * (1 - d) * (age))) ^ (1 / (1 - d))]
    StGroup[, wwa := wa ^ d]
    
    #Calculate the relative number of animals at age a
    for(j in 1:last){
      #Grab the first and last month within the stanza
      a <- stanzafile[StGroupNum == i & Stanza == j, First]:
           stanzafile[StGroupNum == i & Stanza == j, Last]
      
      #Convert Z to a monthly Z
      month.z <- stanzafile[StGroupNum == i & Stanza == j, Z] / 12
      
      StGroup[age %in% a, z := month.z]
      for(age in a[1]:(a[length(a)] + 1)){
        StGroup[age == age, la := prev.surv[j] * (exp(-z) ^ ((age + 1) - a[1]))]
        
      }
      
      StGroup[age %in% a, lawa  := la * wa]
      StGroup[age %in% a, lawwa := la * wwa]
      
      #Numerator for the relative biomass/consumption calculations
      wa.num  <- StGroup[age %in% a, sum(lawa)]
      wwa.num <- StGroup[age %in% a, sum(lawwa)]
      
      stanzafile[StGroupNum == i & Stanza == j, bs.num := wa.num]
      stanzafile[StGroupNum == i & Stanza == j, qs.num := wwa.num]
      
      prev.surv[j + 1] <- StGroup[age == a[length(a)], la]
    }
    #Save relative calculations for plot routine
    stanza.out[[i]] <- StGroup
    
    #Calculate relative biomass
    stanzafile[StGroupNum == i, bs.denom := sum(bs.num)]
    stanzafile[StGroupNum == i, bs := bs.num / bs.denom]
    
    #Calculate relative consumption
    stanzafile[StGroupNum == i, qs.denom := sum(qs.num)]
    stanzafile[StGroupNum == i, qs := qs.num / qs.denom]
    
    #Use leading group to calculate other biomasses
    stanzafile[StGroupNum == i & Leading == T, 
               Biomass := modfile[Group == stanzafile[StGroupNum == i & 
                                                        Leading == T, Group], Biomass]]
    B <- stanzafile[StGroupNum == i & Leading == T, Biomass / bs]
    stanzafile[StGroupNum == i, Biomass := bs * B]
    
    #Use leading group to calculate other consumption
    stanzafile[StGroupNum == i & Leading == T, 
               Cons := modfile[Group == stanzafile[StGroupNum == i & Leading == T, 
                                                   Group], QB] *
                 modfile[Group == stanzafile[StGroupNum == i & Leading == T, Group],
                         Biomass]]
    Q <- stanzafile[StGroupNum == i & Leading == T, Cons / qs]
    stanzafile[StGroupNum == i, Cons := qs * Q]
    stanzafile[, QB := Cons / Biomass]
  }
  
  #Drop extra columns
  stanzafile[, c('bs.num', 'bs.denom', 'bs', 'qs.num', 'qs.denom', 'qs') := NULL]
  
  #Push biomass to modfile
  for(i in 1:nrow(stanzafile)){
    modfile[Group == stanzafile[i, Group], Biomass := stanzafile[i, Biomass]]
  }
  
  #Push consumption to modfile
  for(i in 1:nrow(stanzafile)){
    modfile[Group == stanzafile[i, Group], QB := stanzafile[i, QB]]
  }
  
  if(output == T) return(stanza.out)
} 