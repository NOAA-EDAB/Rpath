#R Ecosystem
#Tutorial for using Rpath

#User parameters - file locations
#I have a windows machine and a linux machine, hence the windows toggle
if(Sys.info()['sysname']=="Windows"){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
}

#To download Rpath
# #This only needs to be done the first time you run the script
# library(devtools)
# devtools::install_github('slucey/RpathDev/Rpath')

library(Rpath); library(data.table)

#Groups and types for the R Ecosystem

groups <- c('Seabirds', 'Whales', 'Seals', 'JuvRoundfish1', 'AduRoundfish1', 
            'JuvRoundfish2', 'AduRoundfish2', 'JuvFlatfish1', 'AduFlatfish1',
            'JuvFlatfish2', 'AduFlatfish2', 'OtherGroundfish', 'Foragefish1',
            'Foragefish2', 'OtherForagefish', 'Megabenthos', 'Shellfish',
            'Macrobenthos', 'Zooplankton', 'Phytoplankton', 'Detritus', 
            'Discards', 'Trawlers', 'Midwater', 'Dredgers')
stgroups <- c(rep(NA, 3), rep('Roundfish1', 2), rep('Roundfish2', 2), 
              rep('Flatfish1', 2), rep('Flatfish2', 2), rep(NA, 14))

types  <- c(rep(0, 19), 1, rep(2, 2), rep(3, 3))

#-------------------------------------------------------------------------------
#Create Model File
REco.params <- create.rpath.param(group = groups, type = types, stgroup = stgroups,
                                  nstanzas = rep(2,4))

#Fill appropriate columns
#Model
biomass <- c(0.0149, 0.454, NA, NA, 1.39, NA, 5.553, NA, 5.766, NA,
             0.739, 7.4, 5.1, 4.7, 5.1, NA, 7, 17.4, 23, 10, rep(NA, 5))

pb <- c(0.098, 0.031, 0.100, 2.026, 0.42, 2.1, 0.425, 1.5, 0.26, 1.1, 0.18, 0.6,
        0.61, 0.65, 1.5, 0.9, 1.3, 7, 39, 240, rep(NA, 5))

qb <- c(76.750, 6.976, 34.455, NA, 2.19, NA, 3.78, NA, 1.44, NA, 1.69,
        1.764, 3.52, 5.65, 3.6, 2.984, rep (NA, 9))

REco.params$model[, Biomass := biomass]
REco.params$model[, PB := pb]
REco.params$model[, QB := qb]

#EE for groups w/o biomass
REco.params$model[Group %in% c('Seals', 'Megabenthos'), EE := 0.8]

#Production to Consumption for those groups without a QB
REco.params$model[Group %in% c('Shellfish', 'Zooplankton'), ProdCons:= 0.25]
REco.params$model[Group == 'Macrobenthos', ProdCons := 0.35]

#Biomass accumulation and unassimilated production
REco.params$model[, BioAcc  := c(rep(0, 22), rep(NA, 3))]
REco.params$model[, Unassim := c(rep(0.2, 18), 0.4, rep(0, 3), rep(NA, 3))]

#Detrital Fate
REco.params$model[, Detritus := c(rep(1, 20), rep(0, 5))]
REco.params$model[, Discards := c(rep(0, 22), rep(1, 3))]

#Fisheries
#Landings
trawl  <- c(rep(0, 4), 0.08, 0, 0.32, 0, 0.09, 0, 0.05, 0.2, rep(0, 10), rep(NA, 3))
mid    <- c(rep(0, 12), 0.3, 0.08, 0.02, rep(0, 7), rep(NA, 3))
dredge <- c(rep(0, 15), 0.1, 0.5, rep(0, 5), rep(NA, 3))
REco.params$model[, Trawlers := trawl]
REco.params$model[, Midwater := mid]
REco.params$model[, Dredgers := dredge]

#Discards
trawl.d  <- c(1e-5, 1e-7, 0.001, 0.001, 0.005, 0.001, 0.009, 0.001, 0.04, 0.001,
              0.01, 0.08, 0.001, 0.001, 0.001, rep(0, 7), rep(NA, 3))
mid.d    <- c(rep(0, 2), 0.001, 0.001, 0.01, 0.001, 0.01, rep(0, 4), 0.05, 0.05,
              0.01, 0.01, rep(0, 7), rep(NA, 3))
dredge.d <- c(rep(0, 3), 0.001, 0.05, 0.001, 0.05, 0.001, 0.05, 0.001, 0.01, 0.05,
              rep(0, 3), 0.09, 0.01, 1e-4, rep(0, 4), rep(NA, 3))
REco.params$model[, Trawlers.disc := trawl.d]
REco.params$model[, Midwater.disc := mid.d]
REco.params$model[, Dredgers.disc := dredge.d]

#Calculate the multistanza biomass/consumption
#Group parameters
REco.params$stanzas$stgroups[, VBGF_Ksp := c(0.145, 0.295, 0.0761, 0.112)]
REco.params$stanzas$stgroups[, Wmat     := c(0.0769, 0.561, 0.117,  0.321)]

#Individual stanza parameters
REco.params$stanzas$stanzas[, First   := c(rep(c(0, 24), 3), 0, 48)]
REco.params$stanzas$stanzas[, Last    := c(rep(c(23, 400), 3), 47, 400)]
REco.params$stanzas$stanzas[, Z       := c(2.026, 0.42, 2.1, 0.425, 1.5, 
                                           0.26, 1.1, 0.18)]
REco.params$stanzas$stanzas[, Leading := rep(c(F, T), 4)]

REco.params <- rpath.stanzas(REco.params)

#Plot multistanza plots
stanzaplot(REco.params, StanzaNum = 1)
stanzaplot(REco.params, StanzaNum = 2)
stanzaplot(REco.params, StanzaNum = 3)
stanzaplot(REco.params, StanzaNum = 4)

#-------------------------------------------------------------------------------
#Modify Diet File
REco.params$diet[, Seabirds        := c(rep(NA, 11), 0.1, 0.25, 0.2, 0.15, 
                                         rep(NA, 6), 0.3)]
REco.params$diet[, Whales          := c(rep(NA, 3), 0.01, NA, 0.01, NA, 0.01, 
                                         NA, 0.01, rep(NA, 4), 0.1, rep(NA, 3), 
                                         0.86, rep(NA, 3))]
REco.params$diet[, Seals           := c(rep(NA, 3), 0.05, 0.1, 0.05, 0.2, 0.005, 
                                         0.05, 0.005, 0.01, 0.24, rep(0.05, 4), 
                                         0.09, rep(NA, 5))]
REco.params$diet[, JuvRoundfish1   := c(rep(NA, 3), rep(c(1e-4, NA), 4), 1e-3, 
                                         rep(NA, 2), 0.05, 1e-4, NA, .02, 0.7785, 
                                         0.1, 0.05, NA)]
REco.params$diet[, AduRoundfish1   := c(rep(NA, 5), 1e-3, 0.01, 1e-3, 0.05, 1e-3, 
                                         0.01, 0.29, 0.1, 0.1, 0.347, 0.03, NA, 
                                         0.05, 0.01, rep(NA, 3))]
REco.params$diet[, JuvRoundfish2   := c(rep(NA, 3), rep(c(1e-4, NA), 4), 1e-3, 
                                         rep(NA, 2), 0.05, 1e-4, NA, .02, 0.7785, 
                                         0.1, .05, NA)]
REco.params$diet[, AduRoundfish2   := c(rep(NA, 3), 1e-4, NA, 1e-4, NA, rep(1e-4, 4), 
                                         0.1, rep(0.05, 3), 0.2684, 0.01, 0.37, 0.001, 
                                         NA, 0.1, NA)]
REco.params$diet[, JuvFlatfish1    := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 3), 
                                         rep(1e-4, 2), NA, 0.416, 0.4334, 0.1, 0.05, 
                                         NA)]
REco.params$diet[, AduFlatfish1    := c(rep(NA, 7), rep(1e-4, 5), rep(NA, 2), 0.001, 
                                         0.05, 0.001, 0.6, 0.2475, NA, 0.1, NA)]
REco.params$diet[, JuvFlatfish2    := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 3),
                                         rep(1e-4, 2), NA, 0.416, 0.4334, 0.1, 0.05, 
                                         NA)]
REco.params$diet[, AduFlatfish2    := c(rep(NA, 7), 1e-4, NA, 1e-4, rep(NA, 4), 
                                         rep(1e-4, 3), 0.44, 0.3895, NA, 0.17, NA)]
REco.params$diet[, OtherGroundfish := c(rep(NA, 3), rep(1e-4, 8), 0.05, 0.08, 0.0992, 
                                         0.3, 0.15, 0.01, 0.3, 0.01, rep(NA, 3))]
REco.params$diet[, Foragefish1     := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA)]
REco.params$diet[, Foragefish2     := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA)]
REco.params$diet[, OtherForagefish := c(rep(NA, 3), rep(c(1e-4, NA), 4), rep(NA, 7), 
                                         0.8196, 0.06, 0.12, NA)]
REco.params$diet[, Megabenthos     := c(rep(NA, 15), 0.1, 0.03, 0.55, rep(NA, 2), 0.32,
                                         NA)]
REco.params$diet[, Shellfish       := c(rep(NA, 18), 0.3, 0.5, 0.2, NA)]
REco.params$diet[, Macrobenthos    := c(rep(NA, 16), 0.01, rep(0.2, 2), NA, 0.59, NA)]
REco.params$diet[, Zooplankton     := c(rep(NA, 18), 0.2, 0.6, 0.2, NA)]

#-------------------------------------------------------------------------------
#Ecopath
REco <- rpath(REco.params, 'R Ecosystem')

#There is an rpath method for generic functions print and summary
REco
summary(REco)
 
#print also includes the mortalities toggle
print(REco, morts = T)

#Summary table of parameters or mortalities can be output
write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Parameters.csv', sep = ''))
write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Mortalities.csv', sep = ''), 
            morts = T)

#Webplot plots the resultant food web
set.seed(34)
#Set the plot order to avoid overlapping titles/points (optional)
my.order <- c(c(22, 20, 21), c(17, 19, 18), c(8, 4, 2, 14, 11, 13, 10, 15, 16, 9, 6),
              c(25, 1, 7, 12), c(3, 5, 24), 23)

tiff(file = paste(out.dir, 'R_Ecosystem.tif'), height = 1500, width = 1700, res = 300)
webplot(REco, labels = T, fleets = T, box.order = my.order, label.cex = 0.65)
dev.off()

#Highlight example
set.seed(34)
tiff(file = paste(out.dir, 'Highlight_AduRoundfish1.tif'), height = 1800, width = 2000, 
     res = 300)
webplot(REco, fleets = T, highlight = 5, box.order = my.order, label.cex = 0.65)
dev.off()

#----------------------------------------------------------------------------------------------------------------------
#Ecosim
#Use the Rpath object from above to run a 100 year Ecosim
#Steps are:
# 1 - rism.scenario
# 2 - modify scenario
# 3 - rsim.run
# 4 - view results with rsim.plot

#Test Adams-Bashforth
#A - base run
REco.sim <- rsim.scenario(REco, REco.params, 100)
REco.AB.1 <- rsim.run(REco.sim, method = 'AB', years = 100)
rsim.plot(REco.AB.1, groups[1:22])

#B - double trawling effort
REco.sim <- rsim.scenario(REco, REco.params, 100)
REco.sim$fishing$EFFORT[26:101, 2] <- 2
REco.AB.2 <- rsim.run(REco.sim, method = 'AB', years = 100)
rsim.plot(REco.AB.2, groups[1:22])


#Test Runge-Kutta 4
#A - base run
REco.sim <- rsim.scenario(REco, juvfile, 100)
REco.RK4.1 <- rsim.run(REco.sim, method = 'RK4', years = 100)
rsim.plot(REco.RK4.1, groups[1:22])

#B - double trawling effort
REco.sim <- rsim.scenario(REco, juvfile, 100)
REco.sim$fishing$EFFORT[26:101, 2] <- 2
REco.RK4.2 <- rsim.run(REco.sim, method = 'RK4', years = 100)
rsim.plot(REco.RK4.2, groups[1:22])









# #Write out the basic outputs from ecosim
# write.Rpath.sim(REco.s1, file = paste(out.dir, 'R_Ecosystem_Ecosim_s1.csv', sep = ''))

#Plot Relative biomass for sim run
ecosim.plot(REco.base)

#Compare with EwE
ewe.base       <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_base_Biomass.csv', sep = ''), skip = 9))
ewe.base.catch <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_base_Yield.csv',   sep = ''), skip = 9))
r.base         <- as.data.table(REco.base$out_BB)
r.base.catch   <- as.data.table(REco.base$out_CC * 12)

r.base       <- r.base[1:1200, ]
r.base.catch <- r.base.catch[1:1200, ]

groups <- copy(names(ewe.base)[1:20])

jpeg(file = paste(out.dir, 'R_Ecosystem_compare_EwE_R_Base.jpg', sep = ''),
     height = 1700, width = 1500, res = 200)
opar <- par(mfrow = c(5, 4), mar = c(2, 2, 2, 2), oma = c(2, 4, 0, 4))
for(i in 1:length(groups)){
  setnames(ewe.base,       groups[i], 'V1')
  setnames(r.base,         groups[i], 'V1')
  setnames(ewe.base.catch, groups[i], 'V1')
  setnames(r.base.catch,   groups[i], 'V1')
  
  bio.y.max <- max(ewe.base$V1, r.base$V1) + .1 * max(ewe.base$V1, r.base$V1)
  bio.y.min <- min(ewe.base$V1, r.base$V1) - .1 * min(ewe.base$V1, r.base$V1)
  cat.y.max <- max(ewe.base.catch$V1, r.base.catch$V1) + .12 * max(ewe.base.catch$V1, r.base.catch$V1)
  cat.y.min <- min(ewe.base.catch$V1, r.base.catch$V1) - .08 * min(ewe.base.catch$V1, r.base.catch$V1)
  
  plot(r.base$V1,    col = 'blue', typ = 'l', main = groups[i], xlab = '', ylab = '', ylim = c(bio.y.min, bio.y.max))
  lines(ewe.base$V1, col = 'blue', lty = 2)
  par(new = T)
  plot(r.base.catch$V1,    col = 'red', typ = 'l', xlab = '', ylab = '', axes = F, ylim = c(cat.y.min, cat.y.max))
  lines(ewe.base.catch$V1, col = 'red', lty = 2)
  axis(4)
  
  setnames(ewe.base,       'V1', groups[i])
  setnames(r.base,         'V1', groups[i])
  setnames(ewe.base.catch, 'V1', groups[i])
  setnames(r.base.catch,   'V1', groups[i])
}

mtext(1, text = 'Months', line = 1, outer = T)
mtext(2, text = 'Biomass (blue)', line = 1, outer = T)
mtext(4, text = 'Catch (red)', line = 1, outer = T)

dev.off()

#-------------------------------------------------------------------------------------------------------------------
#Scenario 2 - Increase F on non-stanza group
#Increase F on Shellfish
REco.shell <- copy(REco.init)
REco.shell <- ecosim.run(REco.shell, 0, 25)
REco.shell <- frate.adjust(REco.shell, 'Shellfish', 'Dredgers', 2)
#REco.shell <- frate.adjust(REco.shell, 'Shellfish', 'Dredgers', .14314, multiplier = F, discard = F)
#REco.shell <- frate.adjust(REco.shell, 'Shellfish', 'Dredgers', .00286, multiplier = F, target  = F)
frate.table(REco.shell)
REco.shell <- ecosim.run(REco.shell, 25, 100, init_run = F)

ecosim.plot(REco.shell)

#Compare with EwE
ewe.shell       <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_doubleF_shellfish_Biomass.csv', sep = ''), skip = 9))
ewe.shell.catch <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_doubleF_shellfish_Yield.csv',   sep = ''), skip = 9))
r.shell         <- as.data.table(REco.shell$out_BB)
r.shell.catch   <- as.data.table(REco.shell$out_CC * 12)

r.shell       <- r.shell[1:1200, ]
r.shell.catch <- r.shell.catch[1:1200, ]

groups <- copy(names(ewe.shell)[1:20])

jpeg(file = paste(out.dir, 'R_Ecosystem_compare_EwE_R_doubleF_shellfish.jpg', sep = ''),
     height = 1700, width = 1500, res = 200)
opar <- par(mfrow = c(5, 4), mar = c(2, 2, 2, 2), oma = c(2, 4, 0, 4))
for(i in 1:length(groups)){
  setnames(ewe.shell,       groups[i], 'V1')
  setnames(r.shell,         groups[i], 'V1')
  setnames(ewe.shell.catch, groups[i], 'V1')
  setnames(r.shell.catch,   groups[i], 'V1')
  
  bio.y.max <- max(ewe.shell$V1, r.shell$V1) + .1 * max(ewe.shell$V1, r.shell$V1)
  bio.y.min <- min(ewe.shell$V1, r.shell$V1) - .1 * min(ewe.shell$V1, r.shell$V1)
  cat.y.max <- max(ewe.shell.catch$V1, r.shell.catch$V1) + .12 * max(ewe.shell.catch$V1, r.shell.catch$V1)
  cat.y.min <- min(ewe.shell.catch$V1, r.shell.catch$V1) - .08 * min(ewe.shell.catch$V1, r.shell.catch$V1)
  
  plot(r.shell$V1,    col = 'blue', typ = 'l', main = groups[i], xlab = '', ylab = '', ylim = c(bio.y.min, bio.y.max))
  lines(ewe.shell$V1, col = 'blue', lty = 2)
  par(new = T)
  plot(r.shell.catch$V1,    col = 'red', typ = 'l', xlab = '', ylab = '', axes = F, ylim = c(cat.y.min, cat.y.max))
  lines(ewe.shell.catch$V1, col = 'red', lty = 2)
  axis(4)
  
  setnames(ewe.shell,       'V1', groups[i])
  setnames(r.shell,         'V1', groups[i])
  setnames(ewe.shell.catch, 'V1', groups[i])
  setnames(r.shell.catch,   'V1', groups[i])
}

mtext(1, text = 'Months', line = 1, outer = T)
mtext(2, text = 'Biomass (blue)', line = 1, outer = T)
mtext(4, text = 'Catch (red)', line = 1, outer = T)

dev.off()

#---------------------------------------------------------------------------------------------------------
#Scenario 3 - Increase F on stanza group
#Increase F on Adult roundfish 1
REco.rf1 <- copy(REco.init)
REco.rf1 <- ecosim.run(REco.rf1, 0, 25)
REco.rf1 <- frate.adjust(REco.rf1, 'AduRoundfish1', 'Trawlers', 2)
REco.rf1 <- frate.adjust(REco.rf1, 'AduRoundfish1', 'Midwater', 2)
REco.rf1 <- frate.adjust(REco.rf1, 'AduRoundfish1', 'Dredgers', 2)

frate.table(REco.rf1)
REco.rf1 <- ecosim.run(REco.rf1, 25, 100, init_run = F)

ecosim.plot(REco.rf1)

#Compare with EwE
ewe.rf1       <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_doubleF_rf1_Biomass.csv', sep = ''), skip = 9))
ewe.rf1.catch <- as.data.table(read.csv(paste(out.dir, 'EwE_R_Ecosystem_doubleF_rf1_Yield.csv',   sep = ''), skip = 9))
r.rf1         <- as.data.table(REco.rf1$out_BB)
r.rf1.catch   <- as.data.table(REco.rf1$out_CC * 12)

r.rf1       <- r.rf1[1:1200, ]
r.rf1.catch <- r.rf1.catch[1:1200, ]

groups <- copy(names(ewe.rf1)[1:20])

jpeg(file = paste(out.dir, 'R_Ecosystem_compare_EwE_R_doubleF_rf1_2.jpg', sep = ''),
     height = 1700, width = 1500, res = 200)
opar <- par(mfrow = c(5, 4), mar = c(2, 2, 2, 2), oma = c(2, 4, 0, 4))
for(i in 1:length(groups)){
  setnames(ewe.rf1,       groups[i], 'V1')
  setnames(r.rf1,         groups[i], 'V1')
  setnames(ewe.rf1.catch, groups[i], 'V1')
  setnames(r.rf1.catch,   groups[i], 'V1')
  
  bio.y.max <- max(ewe.rf1$V1, r.rf1$V1) + .1 * max(ewe.rf1$V1, r.rf1$V1)
  bio.y.min <- min(ewe.rf1$V1, r.rf1$V1) - .1 * min(ewe.rf1$V1, r.rf1$V1)
  cat.y.max <- max(ewe.rf1.catch$V1, r.rf1.catch$V1) + .12 * max(ewe.rf1.catch$V1, r.rf1.catch$V1)
  cat.y.min <- min(ewe.rf1.catch$V1, r.rf1.catch$V1) - .08 * min(ewe.rf1.catch$V1, r.rf1.catch$V1)
  
  plot(r.rf1$V1,    col = 'blue', typ = 'l', main = groups[i], xlab = '', ylab = '', ylim = c(bio.y.min, bio.y.max))
  lines(ewe.rf1$V1, col = 'blue', lty = 2)
  par(new = T)
  plot(r.rf1.catch$V1,    col = 'red', typ = 'l', xlab = '', ylab = '', axes = F, ylim = c(cat.y.min, cat.y.max))
  lines(ewe.rf1.catch$V1, col = 'red', lty = 2)
  axis(4)
  
  setnames(ewe.rf1,       'V1', groups[i])
  setnames(r.rf1,         'V1', groups[i])
  setnames(ewe.rf1.catch, 'V1', groups[i])
  setnames(r.rf1.catch,   'V1', groups[i])
}

mtext(1, text = 'Months', line = 1, outer = T)
mtext(2, text = 'Biomass (blue)', line = 1, outer = T)
mtext(4, text = 'Catch (red)', line = 1, outer = T)

dev.off()

#---------------------------------------------------------------------------------------------------------


