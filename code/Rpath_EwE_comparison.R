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
ewe.s1 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s1.csv', sep = '')))
ewe.s2 <- as.data.table(read.csv(paste(data.dir, 'EwE_R_Ecosystem_Ecosim_s2.csv', sep = '')))

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
  ewe.s1.fix[i, StartBio := ewe.s1[2, i, with = F]]
}
