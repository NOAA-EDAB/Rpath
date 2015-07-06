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
# devtools::install_github('slucey/Rpath/Rpath', 
#                           auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath); library(data.table)

#file skeletons
modfile  <- paste(data.dir, 'REco_mod_2.csv',  sep = '')
dietfile <- paste(data.dir, 'REco_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'REco_ped.csv',  sep = '')
juvfile  <- paste(data.dir, 'REco_juv_2.csv',  sep = '')

#Run ecosim on the R Ecosystem parameter files
REco <- ecopath(modfile, dietfile, pedfile, 'R Ecosystem')

# #There is an rpath method for generic functions print and summary
# REco
# summary(REco)
# 
# #print also includes the mortalities toggle
# print(REco, morts = T)
# 
# #Summary table of parameters or mortalities can be output
# write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Parameters.csv', sep = ''))
# write.Rpath(REco, file = paste(out.dir, 'R_Ecosystem_Mortalities.csv', sep = ''), morts = T)
# 
# #Webplot plots the resultant food web
# set.seed(34)
# #Set the plot order to avoid overlapping titles/points (optional)
# my.order <- c(c(22, 20, 21), c(17, 19, 18), c(8, 4, 2, 14, 11, 13, 10, 15, 16, 9, 6),
#               c(25, 1, 7, 12), c(3, 5, 24), 23)
# 
# tiff(file = paste(out.dir, 'R_Ecosystem.tif'), height = 1500, width = 1700, res = 300)
# webplot(REco, labels = T, fleets = T, box.order = my.order, label.cex = 0.65)
# dev.off()
# 
# #Highlight example
# set.seed(34)
# tiff(file = paste(out.dir, 'Highlight_AduRoundfish1.tif'), height = 1800, width = 2000, res = 300)
# webplot(REco, fleets = T, highlight = 5, box.order = my.order, label.cex = 0.65)
# dev.off()

#----------------------------------------------------------------------------------------------------------------------
#Use the Rpath object from above to run a 100 year Ecosim
#Scenario 1 - Equilibrium
REco.init <- ecosim.init(REco, years = 100, juvfile)
REco.init$NoIntegrate[2]  <- 1
REco.init$NoIntegrate[4]  <- 3
REco.base   <- copy(REco.init)
REco.base   <- ecosim.run(REco.i1, 0, 100)

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

jpeg(file = paste(out.dir, 'R_Ecosystem_compare_EwE_R_doubleF_rf1.jpg', sep = ''),
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


