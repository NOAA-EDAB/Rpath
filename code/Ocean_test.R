#Ocean Test example
#Uses the parameters included with EwE version 5
#Fleets are unique to this test

if(Sys.info()['sysname']=="Windows"){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#To download Rpath
#This only needs to be done the first time you run the script
# library(devtools)
# devtools::install_github('slucey/Rpath/Rpath', 
#                          auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath); library(data.table)

#Generate the parameter files using create.rpath.param() 
groups <- c('Apex pred', 'Juv Apex', 'Mesopelagics', 'Epipelagics', 'Benthic fishes', 
            'Benthopelagics', 'Zooplankton Lg', 'Benthos', 'Microzooplankton',
            'Phytoplankton', 'Detritus', 'Longline', 'Trawlers')
types <- c(rep(0, 9), 1, 2, 3, 3)
model <- create.rpath.param(group = groups, type = types, parameter = 'model')
model[, Biomass := c(0.055, 0.055, 2.533, 0.516, 1.388, 0.022, 9.864, 4.772, 2.434, 0.9, 1.0, NA, NA)]
model[, PB := c(1.157, 2.434, 0.607, 1.991, 0.074, 0.104, 0.466, 0.108, 19.812, 393.435, NA, NA, NA)]
model[, QB := c(14.951, '', 2.748, 9.230, 0.324, 0.431, 2.684, 0.382, 96.561, NA, NA, NA, NA)]
model[Group == 'Juv Apex', ProdCons := 0.239]
model[, BioAcc := c(rep(0, 11), NA, NA)]
model[, Unassim := c(rep(0.2, 9), 0, 0, NA, NA)]
model[, Detritus := c(rep(1, 10), 0, 1, 1)]
model[Group == 'Apex pred',      Longline := .005]
model[Group == 'Mesopelagics',   Longline := .05]
model[Group == 'Epipelagics',    Longline := .005]
model[Group == 'Benthic fishes', Trawlers := .03]
model[Group == 'Apex pred',      Longline.disc := .0001]
model[Group == 'Juv Apex',       Longline.disc := .0001]
model[Group == 'Mesopelagics',   Longline.disc := .0001]
model[Group == 'Epipelagics',    Longline.disc := .0001]
model[Group == 'Benthic fishes', Trawlers.disc := .0001]
write.csv(model, file = paste(data.dir, 'Ocean_mod.csv', sep = ''), row.names = F)

diet <- create.rpath.param(group = groups, type = types, parameter = 'diet')
diet[, 'Apex pred'      := c(NA, 0.05, 0.1, 0.75, NA, NA, 0.1, rep(NA, 4))]
diet[, 'Juv Apex'       := c(NA, NA, 0.05, 0.05, NA, NA, 0.9, rep(NA, 4))]
diet[, Mesopelagics     := c(NA, NA, 0.1, 0.05, NA, NA, 0.25, NA, 0.6, NA, NA)]
diet[, Epipelagics      := c(NA, NA, 0.1, NA, NA, NA, 0.4, NA, 0.4, 0.1, NA)]
diet[, 'Benthic fishes' := c(rep(NA, 4), 0.15, NA, NA, 0.4, NA, NA, 0.45)]
diet[, Benthopelagics   := c(rep(NA, 5), 0.15, 0.2, 0.05, NA, NA, 0.6)]
diet[, 'Zooplankton Lg' := c(rep(NA, 8), 0.6, 0.4, NA)]
diet[, Benthos          := c(rep(NA, 7), 0.05, NA, NA, 0.95)]
diet[, Microzooplankton := c(rep(NA, 9), 1.0, NA)]
write.csv(diet, file = paste(data.dir, 'Ocean_diet.csv', sep = ''), row.names = F)

#Need an empty juvenile file so create with a dummy variable and then remove
juv <- create.rpath.param(group = 'Fish1', parameter = 'juvenile')
juv <- juv[StanzaName != 'Fish1', ]
write.csv(juv, file = paste(data.dir, 'Ocean_juv.csv', sep = ''), row.names = F)

ped <- create.rpath.param(group = groups, type = types, parameter = 'pedigree')
write.csv(ped, file = paste(data.dir, 'Ocean_ped.csv', sep = ''), row.names = F)

#file skeletons
modfile  <- paste(data.dir, 'Ocean_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'Ocean_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'Ocean_ped.csv',  sep = '')
juvfile  <- paste(data.dir, 'Ocean_juv.csv',  sep = '')

#Run ecosim on the R Ecosystem parameter files
ocean <- ecopath(modfile, dietfile, pedfile, 'Ocean Test System')
