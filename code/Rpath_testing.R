#Rpath_testing.R
#Testing of new Rpath functions
#SML
#12/14

#User parameters

windows <- T
if(windows == T){
  r.dir    <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\code\\"
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
  stable   <- "L:\\PhD\\Rpath\\code\\"
}
if(windows == F){
  r.dir    <- "/home/slucey/slucey/Rpath/code/"
  data.dir <- "/home/slucey/slucey/Rpath/data/"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs/"
  stable   <- "/home/slucey/slucey/PhD/Rpath/code/"
}

#To download Rpath
# library(devtools)
# devtools::install_github('slucey/Rpath/Rpath', 
#                           auth_token = 'd95526d2fb3f6e34f9c8481b1740e0033ac1d623')

library(Rpath)

#This is a balanced (and much larger) example with stanzas
#GOA testing
modfile  <- paste(data.dir, 'GOA.csv',     sep = '')
dietfile <- paste(data.dir, 'GOAdiet.csv', sep = '')
pedfile  <- paste(data.dir, 'GOAped.csv',  sep = '')
GOA.juvfile  <- paste(data.dir, 'GOAjuv.csv',  sep = '')

GOA <- ecopath(modfile, dietfile, pedfile, eco.name = 'Gulf of Alaska')

GOA
summary(GOA)

png(file = paste(out.dir, 'Webplot_example_1.png'), height = 1500, width = 1700, res = 200)
webplot(GOA, fleets = T)
dev.off()

png(file = paste(out.dir, 'Webplot_example_2.png'), height = 1500, width = 1700, res = 200)
webplot(GOA, highlight = 33, fleets = T)
dev.off()

png(file = paste(out.dir, 'Webplot_example_3.png'), height = 1500, width = 1700, res = 200)
webplot(GOA, highlight = 33, fleets = T, labels = T, label.num = T, label.pos = 3)
dev.off()

#Ecosim
GOA.base <- ecosim.init(GOA, YEARS = 100, juvfile)
GOA.sim  <- ecosim.run(GOA.base, 0, 100)

#Plot Relative biomass
library(data.table)
biomass <- as.data.table(GOA.sim$out_BB[, 2:23])
rel.bio <- data.table(Month = 1:nrow(biomass))
for(i in 1:ncol(biomass)){
  sp.bio.start <- biomass[1, i, with = F]
  sp.rel.bio   <- biomass[ , i, with = F] / as.numeric(sp.bio.start)
  rel.bio <- cbind(rel.bio, sp.rel.bio)
}

ymax <- max(rel.bio[, 2:ncol(rel.bio), with = F])
xmax <- max(rel.bio[, Month])
plot(0, 0, ylim = c(0, ymax), xlim = c(0, xmax), xlab = '', ylab = '', type = 'n')
line.col <- heat.colors(23)
for(i in 1:(ncol(rel.bio) - 1)){
  lines(rel.bio[, c(1, i + 1), with = F], col = line.col[i])
}




#Microbenchmark verses old version
library(microbenchmark)
ecosim_run_new <- ecosim_run

#Run stable ecosim
source(paste(stable, "ecopath_r_v0_04.r", sep = ''))
source(paste(stable, "ecosim_r_v0_04.r", sep = ''))

# Load ecosim dll                                             
dyn.load(paste(r.dir, 'ecosim.dll', sep = ''))
BaseFile <- modfile
DietFile <- dietfile
PedFile  <- pedfile
JuvFile  <- juvfile

#Ecopath benchmark
path.bm <- microbenchmark(ecopathR(BaseFile, DietFile, PedFile),
                          ecopath(modfile, dietfile, pedfile, eco.name = 'Gulf of Alaska'))

#Ecosim benchmark
ecosim.init.bm <- microbenchmark(load_ecosim( Years = 100, BaseFile, DietFile, PedFile, JuvFile),
                                 ecosim_init(GOA, YEARS = 100, juvfile))

ecosim.bm <- microbenchmark(ecosim_run(base_sim,0,100),
                            ecosim_run_new(GOA.sim, 0, 100))


#R Ecosystem
modfile  <- paste(data.dir, 'REco_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'REco_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'REco_ped.csv',  sep = '')
juvfile  <- paste(data.dir, 'REco_juv.csv',  sep = '')

REco <- ecopath(modfile, dietfile, pedfile, 'R Ecosystem')
REco

summary(REco)

set.seed(34)
my.order <- c(c(3, 23), c(24, 5, 12), c(1, 25, 9, 7, 2), 
              c(4, 13, 6, 8, 14, 11, 15, 10, 16), c(17, 19, 18, 22, 20, 21))

png(file = paste(out.dir, 'R_Ecosystem.png'), height = 1500, width = 1700, res = 300)
webplot(REco, labels = T, fleets = T, box.order = my.order, label.cex = 0.65)
dev.off()

#Highlight example
set.seed(34)
tiff(file = paste(out.dir, 'Highlight_AduFlatfish1.tif'), height = 1800, width = 2000, res = 300)
webplot(REco, fleets = T, highlight = 9, box.order = my.order, label.cex = 0.65)
dev.off()

#Ecosim
REco.sim  <- ecosim.init(REco, YEARS = 100, juvfile)
REco.base <- ecosim.run(REco.sim, 0, 100)

#Plot Relative biomass
library(data.table)
biomass <- as.data.table(REco.base$out_BB[, 2:23])
rel.bio <- data.table(Month = 1:nrow(biomass))
for(i in 1:ncol(biomass)){
  sp.bio.start <- biomass[1, i, with = F]
  sp.rel.bio   <- biomass[ , i, with = F] / as.numeric(sp.bio.start)
  rel.bio <- cbind(rel.bio, sp.rel.bio)
}

ymax <- max(rel.bio[, 2:ncol(rel.bio), with = F])
xmax <- max(rel.bio[, Month])
plot(0, 0, ylim = c(0, ymax), xlim = c(0, xmax), xlab = '', ylab = '', type = 'n')
line.col <- heat.colors(23)
for(i in 1:(ncol(rel.bio) - 1)){
  lines(rel.bio[, c(1, i + 1), with = F], col = line.col[i])
}

#Roundfish 1
lines(rel.bio[, c(1, 5), with = F], col = 'black', lwd = 2)
lines(rel.bio[, c(1, 6), with = F], col = 'red', lwd = 2)

#Roundfish 2
lines(rel.bio[, c(1, 7), with = F], col = 'green3', lwd = 2)
lines(rel.bio[, c(1, 8), with = F], col = 'blue3', lwd = 2)

#Flatfish 1
lines(rel.bio[, c(1, 9),  with = F], col = 'orange3', lwd = 2)
lines(rel.bio[, c(1, 10), with = F], col = 'purple3', lwd = 2)

#Flatfish 2
lines(rel.bio[, c(1, 11), with = F], col = 'brown', lwd = 2)
lines(rel.bio[, c(1, 12), with = F], col = 'pink', lwd = 2)

#Seals
lines(rel.bio[, c(1, 4), with = F], lwd = 3)

#Whales
lines(rel.bio[, c(1, 3), with = F], lwd = 3)

#GOA.juvfile
GOA.juv <- read.csv(GOA.juvfile)
#REco juvfile
REco.juv <- read.csv(juvfile)
