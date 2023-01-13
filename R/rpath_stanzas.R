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
  WageS <- age <- QageS <- Survive <- Z <- survive_L <- bs.num <- qs.num <- Leading <- NULL
  Group <- Biomass <- R <- NageS <- bs.denom <- bs <- qs.denom <- qs <- Cons <- NULL
  QB <- BAB <- Ex <- NULL
  
  #Determine the total number of groups with multistanzas
  Nsplit     <- Rpath.params$stanza$NStanzaGroups
  groupfile  <- Rpath.params$stanza$stgroups
  stanzafile <- Rpath.params$stanza$stindiv
  
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
    nstanzas <- groupfile[StGroupNum == isp, nstanzas]
    stanzafile[StGroupNum == isp & Leading, Last := lastmonth[isp]]
    
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