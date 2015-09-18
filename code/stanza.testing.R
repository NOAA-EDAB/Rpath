#Test multi-stanza inputs/calculations using simple 2 stage from REco
library(data.table)

#First two stanzas from REco
stanzas <- data.table(isp    = c(1, 1, 2, 2),
                      first  = c(0, 25, 0, 25),
                      second = c(24, 400, 24, 400),
                      z      = c(2.026, 0.42, 2.1, 0.425),
                      bab    = c(0, 0, 0, 0))

vbgf <- data.table(isp  = c(1, 2),
                   k    = c(0.145, 0.295),
                   d    = c(0.66667, 0.66667),
                   Wmat = c(0.0577, 0.421))

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
    
