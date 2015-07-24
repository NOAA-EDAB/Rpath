#Rpath_assign.r
#Assign species to Rpath groups
#7/14
#SML

#User parameters
data.dir <- "L:\\EcoAP\\Data\\survey\\"
out.dir  <- "L:\\EcoAP\\Data\\survey\\"

#-------------------------------------------------------------------------------
#Required packages
library(RODBC); library(data.table)

#-------------------------------------------------------------------------------
#User created functions

#-------------------------------------------------------------------------------
#Grab SVSPP data
channel <- odbcDriverConnect()
path.spp <- as.data.table(sqlQuery(channel, "select COMNAME, SCINAME, SVSPP from SVSPECIES_LIST"))
path.spp <- path.spp[!is.na(SVSPP), ]

#Merge EMAX codes
load(paste(data.dir, 'EMAX_groups.RData', sep = ''))
emax[, SCINAME := NULL]
path.spp <- merge(path.spp, emax, by = 'SVSPP', all.x = T)

#Assign Rpath codes
#Drop dead, unknowns, and PR
path.spp <- path.spp[!SVSPP %in% c(0, 300, 400, 406, 408, 410, 412, 414, 519, 520, 605:607) & SVSPP < 950, ] 

#Gadids
path.spp[SVSPP == 69,  RPATH := 'OffHake']
path.spp[SVSPP == 72,  RPATH := 'SilHake']
path.spp[SVSPP == 73,  RPATH := 'Cod']
path.spp[SVSPP == 74,  RPATH := 'Had']
path.spp[SVSPP == 75,  RPATH := 'Pol']
path.spp[SVSPP == 76,  RPATH := 'WhHake']
path.spp[SVSPP == 77,  RPATH := 'RedHake']

#Flats
path.spp[SVSPP == 105, RPATH := 'YTFl']
path.spp[SVSPP == 103, RPATH := 'SumFl']
path.spp[SVSPP == 106, RPATH := 'WinFl']
path.spp[SVSPP == 101, RPATH := 'AtlHal']
path.spp[SVSPP %in% c(99, 102, 104, 107, 108), RPATH := 'OthFlats']
path.spp[SVSPP %in% c(98, 100, 109, 110, 117, 118, 221, 773:799, 
                      821, 825:827, 853, 866, 873, 882), RPATH := 'SmFlats']

#Pelagics
path.spp[SVSPP == 121, RPATH := 'AtlMack']
path.spp[SVSPP == 32,  RPATH := 'AtlHer']
path.spp[SVSPP == 36,  RPATH := 'AtlMen']
path.spp[SVSPP %in% c(30, 33, 34), RPATH := 'RivHer']
path.spp[SVSPP == 35,  RPATH := 'Shad']
path.spp[SVSPP %in% c(2, 37, 38, 112, 123:125, 127:129, 203, 204, 381, 382,
                      426, 463:471, 568:585, 701, 702, 744:746, 860, 877), RPATH := 'OthPels']
path.spp[SVSPP %in% c(31, 37, 39, 43:46, 66, 68, 113, 130, 132:134, 138, 175, 
                      181, 205, 208:212, 423, 427:430, 475, 476, 734, 851, 856, 
                      859, 865, 875, 889, 890, 892), RPATH := 'SmPels']
path.spp[SVSPP %in% c(700, 703:708, 747, 938:943), RPATH := 'LgPels']

#Other
path.spp[SVSPP == 197, RPATH := 'Goose']
path.spp[SVSPP == 155, RPATH := 'Red']
path.spp[SVSPP == 139, RPATH := 'StrBass']
path.spp[SVSPP == 135, RPATH := 'Blu']
path.spp[SVSPP == 141, RPATH := 'BSB']
path.spp[SVSPP == 143, RPATH := 'Scup']
path.spp[SVSPP == 131, RPATH := 'But']
path.spp[SVSPP == 136, RPATH := 'AtlCroak']
path.spp[SVSPP == 654, RPATH := 'RedDrum']
path.spp[SVSPP == 145, RPATH := 'Weak']
path.spp[SVSPP %in% c(140, 142, 146:149, 529, 530, 631:653, 655, 657, 858), RPATH := 'OtherSci']
path.spp[SVSPP %in% c(151, 621:624), RPATH := 'Tile']
path.spp[SVSPP %in% c(379, 380), RPATH := 'Sturg']
path.spp[SVSPP == 193, RPATH := 'Pout']

#Elasmobranchs
path.spp[SVSPP == 13, RPATH := 'SmDog']
path.spp[SVSPP == 15, RPATH := 'SpDog']
path.spp[SVSPP %in% 22:23, RPATH := 'LgSkates']
path.spp[SVSPP %in% c(20, 24:28, 368, 370:372, 924), RPATH := 'SmSkates']
path.spp[SVSPP %in% c(4, 5, 18, 19, 21, 29, 270:272, 367, 373:378), RPATH := 'Rays']
path.spp[SVSPP %in% c(3, 6:12, 14, 16, 17, 350:366, 600:603, 925:937), RPATH := 'Sharks']

#Mesopelagics/Inverts
path.spp[SVSPP %in% c(47, 51:56, 58, 59, 90:93, 111, 114, 126, 137, 144, 
                      150, 154, 158, 227:229, 231:268, 280),  RPATH := 'MesoPels']
path.spp[SVSPP == 301, RPATH := 'AmLob']
path.spp[SVSPP == 306, RPATH := 'Shrimp']
path.spp[SVSPP %in% c(285:287, 291:299, 304, 305, 307, 316, 910:915), RPATH := 'OthShrimp']
path.spp[SVSPP == 310, RPATH := 'RedCrab']
path.spp[SVSPP %in% c(302, 303, 308, 309, 311:315, 317:322, 324:329, 332:334, 339:343, 516:518, 604), RPATH := 'Megaben']
path.spp[SVSPP == 401, RPATH := 'AtlScal']
path.spp[SVSPP %in% c(403, 409), RPATH := 'Clams']
path.spp[SVSPP %in% c(323, 330, 331, 335:338, 344:349, 402, 404, 405, 407, 411, 413, 415:420), RPATH := 'Macroben']
path.spp[SVSPP == 502, RPATH := 'Illex']
path.spp[SVSPP == 503, RPATH := 'Loligo']
path.spp[SVSPP %in% c(501, 504:515), RPATH := 'OthSquid']

#Other demersal
path.spp[is.na(RPATH), RPATH := 'OthDem']

path.spp[, RPATH := as.factor(RPATH)]

#Add NESPP3
cfspp <- as.data.table(sqlQuery(channel, "select SVSPP, NESPP3, NAFOSPP from CFSPP"))
odbcClose(channel)

path.spp <- merge(path.spp, cfspp, by = 'SVSPP', all.x = T, allow.cartesian = T)

setcolorder(path.spp, c('RPATH', 'COMNAME', 'SCINAME', 'SVSPP', 'NESPP3', 'NAFOSPP', 'EMAX'))
setkey(path.spp, RPATH, SVSPP, NESPP3, NAFOSPP)
path.spp <- unique(path.spp)

save(path.spp, file = paste(out.dir, "Rpath_groups.RData", sep = ''))
write.csv(path.spp, file = paste(out.dir, "Rpath_groups.csv", sep = ''))