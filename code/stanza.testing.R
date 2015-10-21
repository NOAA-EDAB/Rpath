#Test multi-stanza inputs/calculations using simple 2 stage from REco
library(data.table)

#First two stanzas from REco
groups <- data.table(GroupNum = c(1, 2),
                     Group    = c('Roundfish1', 'Roundfish2'),
                     nstanzas = c(2, 2),
                     VBGF_Ksp = c(0.145, 0.295),
                     VBGF_d   = c(0.66667, 0.66667),
                     Wmat     = c(0.0577, .421))

stanzas <- data.table(GroupNum = c(1, 1, 2, 2),
                      Stanza   = c(1, 2, 1, 2),
                      First    = c(0, 24, 0, 24),
                      Last     = c(23, 400, 23, 400),
                      Z        = c(2.026, 0.42, 2.1, 0.425),
                      Leading  = c(F, T, F, T),
                      Biomass  = c(0, 1.39, 0, 5.553))

ngroup <- max(groups[, GroupNum])

#Vector of survival rates from 1 stanza to the next
prev.surv <- 1


for(i in 1:ngroup){
  
  #Calculate last month adult
  stan.last <- groups[GroupNum == i, nstanzas]
  #Convert to generalized k from Ksp and make monthly
  k <- (groups[GroupNum == i, VBGF_Ksp] * 3) / 12
  d <- groups[GroupNum == i, VBGF_d]
  #Months to get to 90% Winf
  t90 <- floor(log(1 - 0.9^(1 - d)) / (-1 * k * (1 - d)))
  stanzas[GroupNum == i & Stanza == stan.last, Last := t90]
  
  #Calculate the relative number of animals at age a
  for(j in 1:stan.last){
    #Grab the first and last month within the stanza
    a <- stanzas[GroupNum == i & Stanza == j, First]:stanzas[GroupNum == i & 
                                                               Stanza == j, Last]
    #Convert Z to a monthly Z
    z <- stanzas[GroupNum == i & Stanza == j, Z] / 12
    
    la <- data.table(age  = a,
                     z    = z,
                     surv = prev.surv[j] * (exp(-z) ^ ((a + 1) - a[1])))

    la[, wa := (1 - exp(-k * (1 - d) * (age))) ^ (1 / (1 - d))]
    la[, lawa := surv * wa]
    stanzas[GroupNum == i & Stanza == j, bs.num := la[, sum(lawa)]]
    prev.surv[j + 1] <- la[nrow(la), surv]
  }
  
  stanzas[GroupNum == i, bs.denom := sum(bs.num)]
  stanzas[GroupNum == i, bs := bs.num / bs.denom]
  
  #Use leading group to calculate other biomasses
  B <- stanzas[GroupNum == i & Leading == T, Biomass / bs]
  stanzas[GroupNum == i, Biomass := bs * B]
}

stanzas[, c('bs.num', 'bs.denom', 'bs') := NULL]
r.stanzas <- copy(stanzas)

#EwE version
#Note EwE code indexes to 0 but R can not...so all indexes for Survive are one
#higher in this code
for(i in 1:ngroup){
  BaseStanza <- stanzas[GroupNum == i & Leading == T, Stanza]
  Bio <- stanzas[GroupNum == i & Leading == T, Biomass]
  k <- (groups[GroupNum == i, VBGF_Ksp] * 3) / 12
  d <- groups[GroupNum == i, VBGF_d]
  
  SplitWage <- c()
  for(Age in 0:stanzas[GroupNum == i, max(Last)]){
    SplitWage[Age] <- (1 - exp(-k * (1 - d) * (Age))) ^ (1 / (1 - d))  
  }
  
  Survive <- 1
  PrevSurv <- 1
  
  for(Grp in 1:groups[GroupNum == i, nstanzas]){
    Surv <- exp(-1 * stanzas[GroupNum == i & Stanza == Grp, Z] / 12)
    first  <- stanzas[GroupNum == i & Stanza == Grp, First]
    second <- stanzas[GroupNum == i & Stanza == Grp, Last]
    if(Surv > 0){
      if(first > 0){
        Survive[first + 1] <- Survive[first] * PrevSurv
      }
      for(Age in (first + 1):second){
        Survive[Age + 1] <- Survive[Age] * Surv
      }
      PrevSurv <- Surv
    }
  }
  #Not sure why this is in the code
  if(Surv < 1) Survive[Age] <- Survive[Age] / (1 - Surv)
  
  SumB <- 0
  first <- stanzas[GroupNum == i & Stanza == BaseStanza, First]
  second <- stanzas[GroupNum == i & Stanza == BaseStanza, Last]
  for(Age in first:second){
    SumB <- SumB + Survive[Age + 1] * SplitWage[Age]
  }
  
  Recruits <- Bio / SumB
  
  SplitNo <- c()
  for(Age in 0:stanzas[GroupNum == i, max(Last)]){
    SplitNo[Age] <- Recruits * Survive[Age + 1]
  }
  
  for(Grp in 1:groups[GroupNum == i, nstanzas]){
    first  <- stanzas[GroupNum == i & Stanza == Grp, First]
    second <- stanzas[GroupNum == i & Stanza == Grp, Last]
    biot <- 0
    for(Age in first:second){
      if(Age > 0) biot <- biot + SplitNo[Age] * SplitWage[Age]
    }
    stanzas[GroupNum == i & Stanza == Grp, Biomass := biot]
  } 
}


r.stanzas
stanzas






