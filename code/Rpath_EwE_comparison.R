#Compare EwE and Rpath 

#User parameters
windows <- F
if(windows == T){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}
if(windows == F){
  data.dir <- "/home/slucey/slucey/Rpath/outputs/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#------------------------------------------------------------------
#User created functions
ewe.2.rpath <- function(x){
  ewe <- as.data.table(read.csv(x, skip = 9))
  ewe[, Outside := 0]
  setcolorder(ewe, c('Outside', names(ewe)[which(names(ewe) != 'Outside')]))
  Bmat <- as.matrix(ewe)
  out <- list(spname = names(ewe), out_BB = Bmat)
  return(out)
}

#------------------------------------------------------------------
#Required Packages
library(data.table)

#------------------------------------------------------------------
#Rpath
rpath.REco <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Parameters.csv', sep = '')))
rpath.REco.mort <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Mortalities.csv', sep = '')))
rpath.s1 <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Ecosim_s1.csv', sep = '')))
rpath.s2 <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Ecosim_s2.csv', sep = '')))
rpath.s3 <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Ecosim_s3.csv', sep = '')))

#EwE
ewe.REco <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Parameters.csv', sep = '')))
ewe.REco.mort <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_1.csv', sep = '')))
ewe.REco.mort.pred  <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_2.csv', sep = '')))
ewe.REco.mort.fleet <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_3.csv', sep = '')))
ewe.s1 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s1_Biomass.csv', sep = ''), skip = 9))
ewe.s2 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s2_Biomass.csv', sep = ''), skip = 9))
ewe.s3 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s3_Biomass.csv', sep = ''), skip = 9))
ewe.s1.catch <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s1_Yield.csv', sep = ''), skip = 9))
ewe.s2.catch <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s2_Yield.csv', sep = ''), skip = 9))
ewe.s3.catch <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_s3_Yield.csv', sep = ''), skip = 9))

#To plot EwE scenarios
ewe.s1.list <- ewe.2.rpath(paste(data.dir, 'EwE_R_Ecosystem_s1_Biomass.csv', sep = ''))
ewe.s2.list <- ewe.2.rpath(paste(data.dir, 'EwE_R_Ecosystem_s2_Biomass.csv', sep = ''))
ewe.s3.list <- ewe.2.rpath(paste(data.dir, 'EwE_R_Ecosystem_s3_Biomass.csv', sep = ''))

#Clean up EwE outputs
#Parameters
ewe.REco[, c(1, 4, 5, 12) := NULL, with = F]
setnames(ewe.REco, names(ewe.REco), c('Group', 'TL', 'Biomass', 'Z', 'PB', 'QB', 'EE', 'GE'))
ewe.REco <- ewe.REco[!Group %in% c('Roundfish1', 'Roundfish2', 'Flatfish1', 'Flatfish2'), ]
ewe.REco[is.na(PB), PB := Z]
ewe.REco[, Z := NULL]

#Mortalities
ewe.REco.mort <- ewe.REco.mort[, list(Group.name, Prod.biom.or.Z, X..Other.mort..rate...year.)]
setnames(ewe.REco.mort, names(ewe.REco.mort), c('Group', 'PB', 'M0'))
ewe.REco.mort <- ewe.REco.mort[!Group %in% c('Roundfish1', 'Roundfish2', 'Flatfish1', 'Flatfish2'), ]
setnames(ewe.REco.mort.fleet, names(ewe.REco.mort.fleet), c('X', 'Group', 'F.Trawlers', 'F.Midwater', 'F.Dredgers', 'X2'))
ewe.REco.mort.fleet[, c('X', 'X2') := NULL]
ewe.REco.mort <- merge(ewe.REco.mort, ewe.REco.mort.fleet, by = 'Group', all = T)
groups <- ewe.REco.mort.pred[, Prey...predator]
setnames(ewe.REco.mort.pred, names(ewe.REco.mort.pred), c('X', 'Group', paste('M2.', groups, sep = '')))
ewe.REco.mort.pred[, c('X', 'M2.Phytoplankton') := NULL]
ewe.REco.mort <- merge(ewe.REco.mort, ewe.REco.mort.pred, by = 'Group', all = T)

#Scenarios
ewe.s1.summary <- data.table(Group = rpath.s1[, Group])
for(i in 2:nrow(ewe.s1.summary)){
  group <- as.character(rpath.s1[i, Group])
  ewe.s1.summary[i, StartBio := ewe.s1[1, get(group)]]
  ewe.s1.summary[i, EndBio   := ewe.s1[nrow(ewe.s1) - 1, get(group)]]
  ewe.s1.summary[i, StartCatch := ewe.s1.catch[2, get(group)]]
  ewe.s1.summary[i, EndCatch   := ewe.s1.catch[nrow(ewe.s1.catch) - 1, get(group)]]
}
ewe.s2.summary <- data.table(Group = rpath.s2[, Group])
for(i in 2:nrow(ewe.s2.summary)){
  group <- as.character(rpath.s1[i, Group])
  ewe.s2.summary[i, StartBio := ewe.s2[1, get(group)]]
  ewe.s2.summary[i, EndBio   := ewe.s2[nrow(ewe.s2) - 1, get(group)]]
  ewe.s2.summary[i, StartCatch := ewe.s2.catch[2, get(group)]]
  ewe.s2.summary[i, EndCatch   := ewe.s2.catch[nrow(ewe.s2.catch) - 1, get(group)]]
}
ewe.s3.summary <- data.table(Group = rpath.s3[, Group])
for(i in 2:nrow(ewe.s3.summary)){
  group <- as.character(rpath.s1[i, Group])
  ewe.s3.summary[i, StartBio := ewe.s3[1, get(group)]]
  ewe.s3.summary[i, EndBio   := ewe.s3[nrow(ewe.s3) - 1, get(group)]]
  ewe.s3.summary[i, StartCatch := ewe.s3.catch[2, get(group)]]
  ewe.s3.summary[i, EndCatch   := ewe.s3.catch[nrow(ewe.s3.catch) - 1, get(group)]]
}

#Calculate differences
#Ecopath
ecopath.diff <- merge(rpath.REco, ewe.REco, by = 'Group')
ecopath.diff[, Biomass := ((Biomass.x - Biomass.y) / Biomass.y) * 100]
ecopath.diff[, TL      := ((TL.x - TL.y) / TL.y) * 100]
ecopath.diff[, PB      := ((PB.x - PB.y) / PB.y) * 100]
ecopath.diff[, QB      := ((QB.x - QB.y) / QB.y) * 100]
ecopath.diff[, EE      := ((EE.x - EE.y) / EE.y) * 100]
ecopath.diff[, GE      := ((GE.x - GE.y) / GE.y) * 100]
ecopath.diff[, c(paste(c('Biomass', 'TL', 'PB', 'QB', 'EE', 'GE'), '.x', sep = ''), 
                 paste(c('Biomass', 'TL', 'PB', 'QB', 'EE', 'GE'), '.y', sep = ''),
                 c('X', 'type', 'Removals')) := NULL]

png(file = paste(out.dir, 'Ecopath_differences.png', sep = ''), height = 800, width = 1063, res = 300)
opar <- par(mar = c(3, 3, 1, 1))
boxplot(ecopath.diff[, 2:7, with = F], axes = F)
axis(1, at = axTicks(1), labels = c('Biomass', 'TL', 'PB', 'QB', 'EE', 'GE'), cex.axis = 0.7, padj = -1.5)
axis(2, las = T, cex.axis = 0.8, hadj = .7)
box(lwd = 2)
mtext(1, text = 'Ecopath parameters',  line = 1.3)
mtext(2, text = 'Percent difference', line = 2.1)
dev.off()

#Ecosim scenario 1
#Note - EwE catch is an annual step while Rpath is monthly
ecosim.1.diff <- merge(rpath.s1, ewe.s1.summary, by = 'Group')
ecosim.1.diff <- ecosim.1.diff[!Group %in% c('Outside', 'Detritus')]
ecosim.1.diff[, StartBio := ((StartBio.x - StartBio.y) / StartBio.y) * 100]
ecosim.1.diff[, EndBio   := ((EndBio.x   - EndBio.y)   / EndBio.y)   * 100]
ecosim.1.diff[, StartCatch := (((StartCatch.x * 12) - StartCatch.y) / StartCatch.y) * 100]
ecosim.1.diff[, EndCatch   := (((EndCatch.x   * 12) - EndCatch.y)   / EndCatch.y)   * 100]
ecosim.1.diff[, c(paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.x', sep = ''),
                  paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.y', sep = ''),
                  'X') := NULL]

png(file = paste(out.dir, 'Ecosim_differences_s1.png', sep = ''), height = 800, width = 1063, res = 300)
opar <- par(mar = c(3, 3, 1, 1))
boxplot(ecosim.1.diff[, 2:5, with = F], axes = F)
axis(1, at = axTicks(1), labels = c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), cex.axis = 0.7, padj = -1.5)
axis(2, las = T, cex.axis = 0.8, hadj = .7)
box(lwd = 2)
mtext(1, text = 'Ecosim outputs',  line = 1.3)
mtext(2, text = 'Percent difference', line = 2.1)
dev.off()

#Ecosim scenario 2
#Note - EwE catch is an annual step while Rpath is monthly
ecosim.2.diff <- merge(rpath.s2, ewe.s2.summary, by = 'Group')
ecosim.2.diff <- ecosim.2.diff[!Group %in% c('Outside', 'Detritus')]
ecosim.2.diff[, StartBio := ((StartBio.x - StartBio.y) / StartBio.y) * 100]
ecosim.2.diff[, EndBio   := ((EndBio.x   - EndBio.y)   / EndBio.y)   * 100]
ecosim.2.diff[, StartCatch := (((StartCatch.x * 12) - StartCatch.y) / StartCatch.y) * 100]
ecosim.2.diff[, EndCatch   := (((EndCatch.x   * 12) - EndCatch.y)   / EndCatch.y)   * 100]
ecosim.2.diff[, c(paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.x', sep = ''),
                  paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.y', sep = ''),
                  'X') := NULL]

png(file = paste(out.dir, 'Ecosim_differences_s2.png', sep = ''), height = 800, width = 1063, res = 300)
opar <- par(mar = c(3, 3, 1, 1))
boxplot(ecosim.2.diff[, 2:5, with = F], axes = F)
axis(1, at = axTicks(1), labels = c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), cex.axis = 0.7, padj = -1.5)
axis(2, las = T, cex.axis = 0.8, hadj = .7)
box(lwd = 2)
mtext(1, text = 'Ecosim outputs',  line = 1.3)
mtext(2, text = 'Percent difference', line = 2.1)
dev.off()

#Ecosim scenario 3
#Note - EwE catch is an annual step while Rpath is monthly
ecosim.3.diff <- merge(rpath.s3, ewe.s3.summary, by = 'Group')
ecosim.3.diff <- ecosim.3.diff[!Group %in% c('Outside', 'Detritus')]
ecosim.3.diff[, StartBio := ((StartBio.x - StartBio.y) / StartBio.y) * 100]
ecosim.3.diff[, EndBio   := ((EndBio.x   - EndBio.y)   / EndBio.y)   * 100]
ecosim.3.diff[, StartCatch := (((StartCatch.x * 12) - StartCatch.y) / StartCatch.y) * 100]
ecosim.3.diff[, EndCatch   := (((EndCatch.x   * 12) - EndCatch.y)   / EndCatch.y)   * 100]
ecosim.3.diff[, c(paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.x', sep = ''),
                  paste(c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), '.y', sep = ''),
                  'X') := NULL]

png(file = paste(out.dir, 'Ecosim_differences_s3.png', sep = ''), height = 800, width = 1063, res = 300)
opar <- par(mar = c(3, 3, 1, 1))
boxplot(ecosim.3.diff[, 2:5, with = F], axes = F)
axis(1, at = axTicks(1), labels = c('StartBio', 'EndBio', 'StartCatch', 'EndCatch'), cex.axis = 0.7, padj = -1.5)
axis(2, las = T, cex.axis = 0.8, hadj = .7)
box(lwd = 2)
mtext(1, text = 'Ecosim outputs',  line = 1.3)
mtext(2, text = 'Percent difference', line = 2.1)
dev.off()
