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
                      First    = c(0, 24, 0, 25),
                      Last     = c(23, 400, 24, 400),
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
                     surv = exp(-z) ^ (a - a[1]))

    la[, wa := (1 - exp(-k * (1 - d) * (age))) ^ (1 / 1 - d)]
    la[, lawa := surv * wa]
    stanzas[GroupNum == i & Stanza == j, bs.num := la[, sum(lawa)]]
    prev.surv[j + 1] <- exp(-z) ^ (max(a) - min(a))
  }
  
  stanzas[GroupNum == i, bs.denom := sum(bs.num)]
  stanzas[GroupNum == i, bs := bs.num / bs.denom]
  
  #Use leading group to calculate other biomasses
  B <- stanzas[GroupNum == i & Leading == T, Biomass / bs]
  stanzas[GroupNum == i, Biomass := bs * B]
}

#Calculate biomass in other stanzas



stanza.out <- data.table(isp = unique(stanzas[, isp]),
                         s1  = as.numeric(NA),
                         s2  = as.numeric(NA))
for(i in unique(stanzas[, isp])){
  ind.stanza <- stanzas[isp == i, ]
  ind.stanza[, month.z := z/12]
  n.stanzas <- nrow(ind.stanza)
  for(j in 1:n.stanzas){
    l <- data.table(a = ind.stanza[j, first]:ind.stanza[j, second])
    l[, z := ind.stanza[j, month.z]]
    l[, bab := ind.stanza[j, bab]]
    l[, la := exp(-1 * z * a)]
    l[, wa := (1-exp(-1 * vbgf[isp == ind.stanza[1, isp], k / 12] * a)) ^ 3]
    l[, bs.num := sum(la*wa)]
    stanza.out[i, j + 1] <- l[1, bs.num]
  }
}
    
stanza.out[, bs.denom := rowSums(.SD), .SDcols = c('s1', 's2')]
stanza.out[, bs1 := s1 / bs.denom]
stanza.out[, bs2 := s2 / bs.denom]

Bio <- data.table(isp = c(1, 2),
                  bio = c(1.39, 5.553))
stanza.bio <- merge(stanza.out, Bio, by = 'isp')

stanza.bio[, B := bio / bs2]
stanza.bio[, B1 := bs1 * B]
stanza.bio[, B2 := bio]


