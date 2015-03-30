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

#------------------------------------------------------------------
#Required Packages
library(data.table)

#------------------------------------------------------------------
#Rpath
rpath.REco <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Parameters.csv', sep = '')))
rpath.REco.mort <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Mortalities.csv', sep = '')))
rpath.s1 <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Ecosim_s1.csv', sep = '')))
rpath.s2 <- as.data.table(read.csv(paste(data.dir, 'R_Ecosystem_Ecosim_s2.csv', sep = '')))

#EwE
ewe.REco <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Parameters.csv', sep = '')))
ewe.REco.mort <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_1.csv', sep = '')))
ewe.REco.mort.pred  <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_2.csv', sep = '')))
ewe.REco.mort.fleet <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Mortalities_3.csv', sep = '')))
ewe.s1 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s1.csv', sep = ''), skip = 1))
ewe.s2 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s2.csv', sep = ''), skip = 1))
ewe.s1.catch <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s1_catch.csv', sep = ''), skip = 1))
ewe.s2.catch <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s2_catch.csv', sep = ''), skip = 1))

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
ewe.s1.fix <- data.table(Group = rpath.s1[, Group])
for(i in 1:nrow(ewe.s1.fix)){
  ewe.s1.fix[i, StartBio := ewe.s1[1, i, with = F]]
  ewe.s1.fix[i, EndBio   := ewe.s1[nrow(ewe.s1) - 1, i, with = F]]
  ewe.s1.fix[i, StartCatch := ewe.s1.catch[2, i, with = F]]
  ewe.s1.fix[i, EndCatch   := ewe.s1.catch[nrow(ewe.s1.catch) - 1, i, with = F]]
}
ewe.s2.fix <- data.table(Group = rpath.s2[, Group])
for(i in 1:nrow(ewe.s2.fix)){
  ewe.s2.fix[i, StartBio := ewe.s2[1, i, with = F]]
  ewe.s2.fix[i, EndBio   := ewe.s2[nrow(ewe.s2) - 1, i, with = F]]
  ewe.s2.fix[i, StartCatch := ewe.s2.catch[2, i, with = F]]
  ewe.s2.fix[i, EndCatch   := ewe.s2.catch[nrow(ewe.s2.catch) - 1, i, with = F]]
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
ecosim.1.diff <- merge(rpath.s1, ewe.s1.fix, by = 'Group')
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
ecosim.2.diff <- merge(rpath.s2, ewe.s2.fix, by = 'Group')
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
