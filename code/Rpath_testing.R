#Rpath_testing.R
#Testing of new Rpath functions
#SML
#12/14

#User parameters

windows <- F
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


#These groups are unbalanced for the time being
#D1F1 testing
modfile  <- paste(data.dir, 'D1F1_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'D1F1_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'AnchovyBay_ped.csv',  sep = '')

d1f1 <- ecopath(modfile, dietfile, pedfile)
summary(d1f1)



#D2F3 testing
modfile  <- paste(data.dir, 'D2F3_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'D2F3_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'D2F3_ped.csv',  sep = '')

d2f3 <- ecopath(modfile, dietfile, pedfile)
summary(d2f3)


#This is a balanced (and much larger) example with stanzas
#GOA testing
modfile  <- paste(data.dir, 'GOA.csv',     sep = '')
dietfile <- paste(data.dir, 'GOAdiet.csv', sep = '')
pedfile  <- paste(data.dir, 'GOAped.csv',  sep = '')
juvfile  <- paste(data.dir, 'GOAjuv.csv',  sep = '')

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
GOA.sim <- ecosim.init(GOA, YEARS = 100, juvfile)
GOA.sim <- ecosim.run(GOA.sim, 0, 100)


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
webplot(REco, labels = T, fleets = T, box.order = my.order)

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

