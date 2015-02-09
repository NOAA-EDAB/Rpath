#Rpath_testing.R
#Testing of new Rpath functions
#SML
#12/14

#User parameters

windows <- F
if(windows == T){
  r.dir    <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath_git\\Rpath\\code\\"
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath_git\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath_git\\Rpath\\outputs\\"
  stable   <- "L:\\PhD\\Rpath\\code\\"
}
if(windows == F){
  r.dir    <- "slucey/Rpath/code/"
  data.dir <- "slucey/Rpath/data/"
  out.dir  <- "slucey/Rpath/outputs/"
  stable   <- "slucey/PhD/Rpath/code/"
}

source(paste(r.dir, 'Rpath_methods.R', sep = ''))
source(paste(r.dir, 'ecopath_r.r',     sep = ''))
source(paste(r.dir, 'ecosim_r.r',      sep = ''))


#These groups are unbalanced for the time being
#D1F1 testing
modfile  <- paste(data.dir, 'D1F1_mod.csv',  sep = '')
dietfile <- paste(data.dir, 'D1F1_diet.csv', sep = '')
pedfile  <- paste(data.dir, 'D1F1_ped.csv',  sep = '')

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

#Ecosim
library(Rcpp)
sourceCpp(paste(r.dir, 'test_ecosim.cpp', sep = ''))

GOA.sim <- ecosim_init(GOA, YEARS = 100, juvfile)
GOA.sim <- ecosim_run(GOA.sim, 0, 100)


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


