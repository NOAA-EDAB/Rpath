#Georges_Bank_Rpath.r
#Compile data from surveys/EMAX for Georges Bank
#8/14
#SML

#User parameters
data.dir <- "L:\\PhD\\Rpath\\data\\"
out.dir  <- "L:\\PhD\\Rpath\\data\\"
r.dir    <- "L:\\Rworkspace\\Survey\\"
gis.dir  <- "L:\\Rworkspace\\GIS_files"

#-------------------------------------------------------------------------------
#Required packages
library(data.table); library(RODBC); library(rgdal)

#-------------------------------------------------------------------------------
#User created functions
count <- function(x){
    num <- rep(1, length(x))
    out <- sum(num)
    return(out)
    }
    
#-------------------------------------------------------------------------------
memory.size(4000)

load(paste(data.dir, "Rpath_survey_data.RData", sep = ''))

#Using 5 year average from 1996 - 2000 similar to EMAX
GB <- rpath.bio[EPU == 'GB' & YEAR %in% 1996:2000, ]
GB.mod <- GB[, mean(biomass.area), by = RPATH]
GB.mod[, Type := 0]
GB.mod[, Biomass := V1 * .001] #convert kg/km^2 to g/m^2
GB.mod[, V1 := NULL]
setnames(GB.mod, "RPATH", "Group")

#Add PB
GB.mod[Group %in% c('AtlHer', 'AtlMack', 'RivHer', 'Shad', 'OthPels'), PB := 0.420]
GB.mod[Group %in% c('AtlScal', 'Clams'), PB := 0.800]
GB.mod[Group %in% c('But', 'SmPels', 'Scup', 'MesoPels', 'Loligo', 'Illex', 'OthSquid'), PB := 0.950]
GB.mod[Group %in% c('AmLob', 'Megaben', 'RedCrab', 'OthShrimp'), PB := 1.500]
GB.mod[Group %in% c('Blu', 'Cod', 'LgSkates', 'Red', 'SmDog', 'SpDog', 'SumFl', 'AtlHal',
                    'Goose', 'Had', 'OffHake', 'OthDem', 'Pol', 'Pout', 'WhHake',
                    'WinFl', 'YTFl', 'StrBass', 'Rays', 'OthFlats', 'BSB'), PB := 0.450]
GB.mod[Group %in% c('SmFlats', 'RedHake', 'SilHake', 'SmSkates', 'OtherSci', 'Weak'), PB := 0.550]
GB.mod[Group == 'Macroben', PB := 2.500]
GB.mod[Group %in% grep('_Juv', GB.mod[, Group], value = T), PB := 15.00]

#Add QB
GB.mod[Group %in% c('AtlHer', 'AtlMack', 'RivHer', 'Shad', 'OthPels'), QB := 2]
GB.mod[Group %in% c('AtlScal', 'Clams'), QB := 10]
GB.mod[Group %in% c('But', 'SmPels', 'Scup'), QB := 2]
GB.mod[Group %in% c('Loligo', 'Illex', 'OthSquid'), QB := 2.75]
GB.mod[Group == 'MesoPels', QB := 1.83]
GB.mod[Group %in% c('AmLob', 'Megaben', 'RedCrab', 'OthShrimp'), QB := 14]
GB.mod[Group %in% c('Blu', 'Cod', 'LgSkates', 'Red', 'SmDog', 'SpDog', 'SumFl', 'AtlHal',
                    'Goose', 'Had', 'OffHake', 'OthDem', 'Pol', 'Pout', 'WhHake',
                    'WinFl', 'YTFl', 'StrBass', 'Rays', 'OthFlats', 'BSB'), QB := 0.83]
GB.mod[Group %in% c('SmFlats', 'RedHake', 'SilHake', 'SmSkates', 'OtherSci', 'Weak'), QB := 0.92]
GB.mod[Group == 'Macroben', QB := 14]
GB.mod[Group %in% grep('_Juv', GB.mod[, Group], value = T), QB := 45]
GB.mod <- GB.mod[Group != 'Sharks', ]

#Add values from EMAX
emax <- data.table(Group   = c('Sharks', 'LgPels', 'Seabirds', 'Seals',
                               'BalWhale', 'ToothWhale', 'Dolphins',
                               'Micronekton', 'GelZoo', 'Zoop', 'Microzoop',
                               'Bact', 'PP'),
                   Type    = c(rep(0, 12), 1),
                   Biomass = c(0.024, 0.035, 0.014, 0.060,
                               0.310, 0.075, 0.045,
                               7.400, 5.200, 1.450, 0.261,
                               0.345, 20.11),
                   PB      = c(0.417, 0.686, 0.286, 0.065,
                               0.040, 0.040, 0.040,
                               15.00, 40.00, 0.190, 0.200,
                               35.00, 163.178),
                   QB      = c(0.550, 2.050, 4.380, 0.025,
                               4.500, 14.40, 14.40,
                               36.50, 146.0, 115.0, 242.4,
                               380.0, 0.000))
                               
#Merge survey data with emax data
GB.mod <- rbindlist(list(GB.mod, emax))

#Add detritus groups
detritus <- data.table(Group = c('Discards', 'Detritus'), Type = 2, B = NaN, PB = NaN, QB = NaN)
GB.mod <- rbindlist(list(GB.mod, detritus))

#--------------------------------------------------------------------------------------
#Catch data
load(paste(data.dir, "Rpath_catch_data.RData", sep = ''))

#Using 5 year average from 1996 - 2000 similar to EMAX
GB.catch <- rpath.catch[EPU == 'GB' & YEAR %in% 1996:2000, ]
setkey(GB.catch, GEAR, RPATH)
GB.mod.land <- GB.catch[, mean(land.area, na.rm = T), by = key(GB.catch)]
GB.mod.disc <- GB.catch[, mean(disc.area, na.rm = T), by = key(GB.catch)]

GB.mod.land[, Landings := V1 * .001] #convert kg/km^2 to g/m^2
GB.mod.disc[, Discards := V1 * .001] #convert kg/km^2 to g/m^2

GB.mod.land[, V1 := NULL]
GB.mod.disc[, V1 := NULL]

GB.mod.catch <- merge(GB.mod.land, GB.mod.disc, by = key(GB.catch))
setnames(GB.mod.catch, "RPATH", "Group")

#Add fleets to bottom of GB.mod
fleets <- data.table(Group = unique(GB.mod.catch[, GEAR]), Type = 3, B = NaN, PB = NaN, QB = NaN)
GB.mod <- rbindlist(list(GB.mod, fleets))

#Add in additional columns
GB.mod[, EE       := NaN]
GB.mod[, ProdCons := NaN]
GB.mod[Type %in% 0:2, BioAcc   := 0]
GB.mod[Type == 0,     Unassim  := 0.2]
GB.mod[Type %in% 1:2, Unassim  := 0]
GB.mod[, DetInput := NaN]
GB.mod[Type == 3, Discards := 1]
GB.mod[Type %in% 0:2, Discards := 0]
GB.mod[Type %in% 0:1, Detritus := 1]
GB.mod[Type %in% 2:3, Detritus := 0]

#Merge landings
for(i in 1:nrow(fleets)){
  fleet.land <- GB.mod.catch[GEAR == fleets[i, Group], list(Group, Landings)]
  setnames(fleet.land, "Landings", as.character(fleets[i, Group]))
  GB.mod <- merge(GB.mod, fleet.land, by = 'Group', all.x = T)
  }

#Merge discards
for(i in 1:nrow(fleets)){
  fleet.disc <- GB.mod.catch[GEAR == fleets[i, Group], list(Group, Discards)]
  setnames(fleet.disc, "Discards", paste(as.character(fleets[i, Group]), "_disc", sep = ''))
  GB.mod <- merge(GB.mod, fleet.disc, by = 'Group', all.x = T)
  }

setkey(GB.mod, Type, Group)

save(GB.mod, file = paste(out.dir, "GB.RData", sep = ''))
write.csv(GB.mod, file = paste(out.dir, "GB.csv", sep = ''), row.names = F)

#------------------------------------------------------------------------------------------------
#Diet File
load(paste(data.dir, "Rpath_diet_data.RData", sep = ''))

GB.diet <- rpath.diet[EPU == 'GB' & YEAR %in% 1996:2000, ]

#Number of stomachs
setkey(GB.diet, CRUISE6, STRATUM, STATION, RPATH)
GB.stom <- unique(GB.diet)
GB.stom[, n     := count(STATION), by = 'RPATH']
GB.stom <- GB.stom[, list(RPATH, n)]
setkey(GB.stom, RPATH)
GB.stom <- unique(GB.stom)
GB.diet <- merge(GB.diet, GB.stom, by = 'RPATH')

#Calculate DC by haul
setkey(GB.diet, CRUISE6, STRATUM, STATION, RPATH)
GB.diet[, TotPW := sum(PYAMTW), by = key(GB.diet)]
GB.diet[, DCh   := PYAMTW / TotPW]

#Calculate overall DC
setkey(GB.diet, RPATH, PYPATH)
GB.diet[, sumDC := sum(DCh), by = key(GB.diet)]
GB.diet[, DC    := sumDC / n]
GB.diet <- unique(GB.diet)
GB.diet <- GB.diet[, list(RPATH, PYPATH, DC)]

#Add data from EMAX for missing diets
emax.diet <- data.table(RPATH  = c(rep('Bact', 2),                                
                                   rep('Microzoop', 4),
                                   rep('Zoop', 6),
                                   rep('GelZoo', 12),
                                   rep('Micronekton', 4),
                                   rep('Macroben', 9),
                                   rep('Megaben', 7),
                                   rep('AtlScal', 3),
                                   rep('Clams', 3),
                                   rep('AmLob', 5),
                                   rep('RedCrab', 5),
                                   rep('OthShrimp', 7),
                                   rep('Illex', 10),
                                   rep('Loligo', 10),
                                   rep('OthSquid', 10),
                                   rep('MesoPels', 6),
                                   rep('SmPels', 6),
                                   rep('RivHer', 5),
                                   rep('OthSci', 15),
                                   rep('SmFlats', 15),                           
                                   rep('Weak', 16),           
                                   rep('AtlHal', 20),
                                   rep('Rays', 15),
                                   rep('Sharks', 43),
                                   rep('LgPels', 7),
                                   rep('Seals', 31),
                                   rep('BalWhale', 13),
                                   rep('ToothWhale', 22),
                                   rep('Dolphins', 21),
                                   rep('Seabirds', 17)),
                        PYPATH = c('PP', 'Detritus',
                                   'PP', 'Bact', 'Microzoop', 'Detritus',
                                   'PP', 'Microzoop', 'Zoop', 'GelZoo', 'Macroben', 'Detritus',
                                   'PP', 'Bact', 'Microzoop', 'Zoop', 'GelZoo', 'AtlHer', 'SmPels', 
                                      'RivHer', 'Illex', 'Loligo', 'OthSquid', 'Detritus',
                                   'PP', 'Zoop', 'Micronekton', 'Detritus',
                                   'PP', 'Bact', 'Zoop', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'Discards', 'Detritus',
                                   'Bact', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'Discards', 'Detritus',
                                   'PP', 'Bact', 'Detritus',
                                   'PP', 'Bact', 'Detritus',
                                   'Bact', 'Macroben', 'Megaben', 'Discards', 'Detritus',
                                   'Bact', 'Macroben', 'Megaben', 'Discards', 'Detritus',
                                   'PP', 'Bact', 'Micronekton', 'Macroben', 'OthShrimp', 'Discards', 'Detritus',
                                   'Zoop', 'Micronekton', 'Macroben', 'OthShrimp', 'AtlHer', 'AtlMack', 'SmPels', 'OthSquid', 'RivHer', 'Juvs',
                                   'Zoop', 'Micronekton', 'Macroben', 'OthShrimp', 'AtlHer', 'AtlMack', 'SmPels', 'OthSquid', 'RivHer', 'Juvs',
                                   'Zoop', 'Micronekton', 'Macroben', 'OthShrimp', 'AtlHer', 'AtlMack', 'SmPels', 'OthSquid', 'RivHer', 'Juvs',
                                   'PP', 'Zoop', 'GelZoo', 'Micronekton', 'Macroben', 'Juvs',
                                   'PP', 'Zoop', 'GelZoo', 'Micronekton', 'Macroben', 'Juvs',
                                   'PP', 'Zoop', 'Micronekton', 'Macroben', 'Juvs',
                                   'GelZoo', 'Micronekton', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'OthShrimp', 'AtlHer', 
                                      'SmPels', 'Illex', 'Loligo', 'OthSquid', 'OthDem', 'Discards', 'Detritus',
                                   'GelZoo', 'Micronekton', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'OthShrimp', 'AtlHer', 
                                      'SmPels', 'Illex', 'Loligo', 'OthSquid', 'OthDem', 'Discards', 'Detritus',
                                   'GelZoo', 'Micronekton', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'OthShrimp', 'Juvs', 'AtlHer', 
                                      'SmPels', 'Illex', 'Loligo', 'OthSquid', 'OthDem', 'Discards', 'Detritus',
                                   'GelZoo', 'Micronekton', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'OthShrimp', 'Juvs', 'AtlHer', 
                                      'SmPels', 'Illex', 'Loligo', 'OthSquid', 'OthFlats', 'SmFlats', 'OthDem', 'AtlCod', 'Had', 
                                      'SilHake', 'RedHake', 
                                   'GelZoo', 'Micronekton', 'Macroben', 'Megaben', 'AtlScal', 'Clams', 'OthShrimp', 'AtlHer', 
                                      'SmPels', 'Illex', 'Loligo', 'OthSquid', 'OthDem', 'Discards', 'Detritus',
                                   'Zoop', 'Macroben', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 'OthSquid', 'RivHer',
                                      'Shad', 'OthPels', 'BSB', 'But', 'Cod', 'OthFlats', 'Had', 'OffHake', 'LgSkates', 'Pol', 'Rays', 
                                      'Red', 'RedHake', 'Scup', 'SilHake', 'SmDog', 'SmSkates', 'SpDog', 'StrBass', 'SumFl', 'Weak', 
                                      'WhHake', 'WinFl', 'YTFl', 'OthDem', 'OtherSci', 'Sharks', 'LgPels', 'Seals', 'BalWhale', 
                                      'ToothWhale', 'Dolphins', 'Seabirds', 'Detritus',
                                   'GelZoo', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 'OthSquid',
                                   'Micronekton', 'AtlHer', 'AtlMack', 'SmPels', 'RivHer', 'Shad', 'OthPels', 'Cod', 'BSB', 'But',  
                                      'OthFlats', 'Had', 'OffHake', 'LgSkates', 'Pol', 'Rays', 'Red', 'RedHake', 'Scup', 'SilHake', 
                                      'SmDog', 'SmSkates', 'SpDog', 'StrBass', 'SumFl', 'Weak', 'WhHake', 'WinFl', 'YTFl', 'OthDem', 
                                      'OtherSci',
                                   'Zoop', 'GelZoo', 'Micronekton', 'Macroben', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 
                                      'OthSquid', 'RivHer', 'Shad', 'Detritus',
                                   'GelZoo', 'Micronekton', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 'OthSquid', 'RivHer',
                                      'Shad', 'OthPels', 'Cod', 'But', 'Had', 'OffHake', 'LgSkates', 'Pol', 'Scup', 'SilHake', 'SmDog',
                                      'SpDog', 'ToothWhale', 
                                   'GelZoo', 'Micronekton', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 'OthSquid', 'RivHer',
                                      'Shad', 'OthPels', 'Cod', 'But', 'Had', 'OffHake', 'LgSkates', 'Pol', 'Scup', 'SilHake', 'SmDog',
                                      'SpDog',
                                   'Zoop', 'Micronekton', 'OthShrimp', 'AtlHer', 'AtlMack', 'SmPels', 'Illex', 'Loligo', 'OthSquid', 'RivHer',
                                      'Shad', 'OthPels', 'But', 'OffHake', 'Scup', 'SilHake', 'Discards'),
                        DC     = c(0.15, 0.85,
                                   0.15, 0.40, 0.10, 0.35,
                                   0.55, 0.15, 0.15, 0.05, 0.01, 0.09,
                                   0.082, 0.02, 0.05, 0.70, 0.02, 0.02, 0.004, 0.001, 0.001, 0.001, 0.001, 0.10, 
                                   0.114, 0.742, 0.029, 0.115,
                                   0.20, 0.155, 0.01, 0.30, 0.01, 0.005, 0.005, 0.005, 0.31,
                                   0.13, 0.52, 0.15, 0.01, 0.02, 0.03, 0.14,
                                   0.60, 0.20, 0.20,
                                   0.60, 0.20, 0.20,
                                   0.12, 0.55, 0.18, 0.03, 0.12,
                                   0.12, 0.55, 0.18, 0.03, 0.12,
                                   0.055, 0.32, 0.11, 0.12, 0.015, 0.05, 0.33,
                                   0.09, 0.45, 0.15, 0.07, 0.005, 0.005, 0.015, 0.06, 0.005, 0.15,     
                                   0.09, 0.45, 0.15, 0.07, 0.005, 0.005, 0.015, 0.06, 0.005, 0.15,
                                   0.09, 0.45, 0.15, 0.07, 0.005, 0.005, 0.015, 0.06, 0.005, 0.15,
                                   0.16, 0.75, 0.03, 0.04, 0.01, 0.01,
                                   0.16, 0.75, 0.03, 0.04, 0.01, 0.01,
                                   0.01, 0.95, 0.02, 0.01, 0.01,
                                   0.005, 0.06, 0.65, 0.05, 0.04, 0.06, 0.02, 0.05, 0.002, 0.001, 0.001, 0.001, 0.04, 0.01, 0.01,
                                   0.005, 0.06, 0.65, 0.05, 0.04, 0.06, 0.02, 0.05, 0.002, 0.001, 0.001, 0.001, 0.04, 0.01, 0.01,
                                   0.01, 0.035, 0.65, 0.04, 0.05, 0.06, 0.005, 0.01, 0.08, 0.002, 0.001, 0.001, 0.001, .025, 0.02, 0.01,
                                   0.001, 0.02, 0.179, 0.01, 0.01, 0.01, 0.15, 0.01, 0.25, 0.08, 0.02, 0.02, 0.02, 0.005, 0.005, 
                                      0.01, 0.05, 0.05, 0.05, 0.05,
                                   0.005, 0.06, 0.65, 0.05, 0.04, 0.06, 0.02, 0.05, 0.002, 0.001, 0.001, 0.001, 0.04, 0.01, 0.01,
                                   0.03, 0.02, 0.13, 0.12, 0.05, 0.05, 0.05, 0.05, 0.02, 0.01, 0.15, rep(0.005, 22), 0.02, 0.01, 
                                      0.03, 0.01, 0.03, 0.01, 0.01, 0.01, 0.03, 0.05,   
                                   0.07, 0.07, 0.07, 0.73, 0.02, 0.02, 0.02,                              
                                   0.08, 0.32, 0.155, 0.154, 0.09, 0.04, 0.001, 0.015, rep(0.005, 21), 0.02, 0.02,
                                   0.52, 0.001, 0.289, 0.10, 0.03, 0.02, 0.01, 0.003, 0.003, 0.004, 0.001, 0.001, 0.018,     
                                   0.003, 0.025, 0.20, 0.16, 0.19, 0.08, 0.08, 0.09, 0.03, 0.02, 0.001, rep(0.012, 10), 0.001,
                                   0.003, 0.025, 0.20, 0.16, 0.19, 0.08, 0.08, 0.09, 0.03, 0.02, 0.001, rep(0.0121, 10),
                                   0.03, 0.12, 0.06, 0.16, 0.12, 0.24, 0.023, 0.023, 0.023, 0.03, 0.01, 0.001, 0.01, 0.01, 0.01, 0.01, 0.12))
                                   
GB.diet <- rbindlist(list(GB.diet, emax.diet))
                                   
#Set up file for Rpath
Groups <- GB.mod[Type !=3, list(Group, Type)]
GB.diet[, RPATH := factor(RPATH, levels = levels(Groups[, Group]))]

GB.diet.output <- data.table(Group = Groups[, Group])
for(i in 1:length(Groups[, Group])){
  Group.dc <- GB.diet[RPATH == Groups[i, Group], list(PYPATH, DC)]
  if(nrow(Group.dc) == 0) Group.dc <- data.table(PYPATH = Groups[, Group], DC = 0)
  setnames(Group.dc, c("PYPATH", "DC"), c("Group", as.character(Groups[i, Group])))
  GB.diet.output <- merge(GB.diet.output, Group.dc, by = 'Group', all.x = T)
  }
  
GB.diet <- merge(GB.diet.output, Groups, by = 'Group')
setkey(GB.diet, Type, Group)

save(GB.diet, file = paste(out.dir, "GB_diet.RData", sep = ''))
write.csv(GB.diet, file = paste(out.dir, "GB_diet.csv", sep = ''), row.names = F)


#add in changing NAs to zero
