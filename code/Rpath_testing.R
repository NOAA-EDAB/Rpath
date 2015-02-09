#Rpath_testing.R
#Testing of new Rpath functions
#SML
#12/14

#User parameters

windows <- T
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



#Run stable ecosim
source(paste(stable, "ecopath_r_v0_04.r", sep = ''))
source(paste(stable, "ecosim_r_v0_04.r", sep = ''))

# Load ecosim dll                                             
dyn.load(paste(r.dir, 'ecosim.dll', sep = ''))

# Basic Ecopath parameters (Biomass, Production, Consumption, Fishing)
BaseFile <- modfile
# Ecopath Diet Matrix
DietFile <- dietfile
# Data pedigree (if no pedigree, use file of shown format, with 1.0 for all values) 
PedFile  <- pedfile
# Juvenile adult split groups.  IF NO JUV/ADU IN YOUR MODEL, use file with the 
# headers in the example, but no rows.  This is required as a placeholder 
# (minor bug needs fixing)
JuvFile  <- juvfile


# EXAMPLES
path <- ecopathR(BaseFile, DietFile, PedFile)  
summary(path) 

# EXAMPLE 2:  Load Ecosim parameters, declaring enough vector length for 100 years.
#             Run once in baseline (equilibrium), run with fishing change.

# Prepare for run by loading from csv files (includes balancing Ecopath)
base_sim <- load_ecosim( Years = 100, BaseFile, DietFile, PedFile, JuvFile)

#long_sim <- load_ecosim( Years = 500, BaseFile, DietFile, PedFile, JuvFile)
# Base run from year 0 to year 100
run0 <- ecosim_run(base_sim,0,100)  

# resulting output biomass
summary(run0$out_BB[34])
summary(GOA.sim$out_BB[34])

plot(GOA.sim$out_BB[34])


test.sim <- list(NUM_LIVING = 2, NUM_DEAD = 1, juv_N = 1, 
                 state_BB = c(1, 3.97, 0.85, 2.456),
                 SpawnBio = c(0, 0.87),
                 JuvNum = c(0, 1), AduNum = c(0, 2), firstMoAdu = 24,
                 WageS = matrix(0, 1, 1),
                 NageS = matrix(0, 1, 1),
                 out_BB =  data.frame(Sp1 = rep(0, 12*5),
                                      Sp2 = rep(0, 12*5),
                                      Sp3 = rep(0, 12*5)),
                 out_SSB = data.frame(Sp1 = rep(0, 12*5),
                                      Sp2 = rep(0, 12*5),
                                      Sp3 = rep(0, 12*5)),
                 out_rec = data.frame(Sp1 = rep(0, 12*5),
                                      Sp2 = rep(0, 12*5),
                                      Sp3 = rep(0, 12*5)))

cppFunction('int ecosim_test(List mod, int StartYear, int EndYear){
            int y, m, dd;
            int sp, i;
            int STEPS_PER_YEAR = 12;

            // Parse out List mod
            int NUM_LIVING = as<int>(mod["NUM_LIVING"]);
            int NUM_DEAD   = as<int>(mod["NUM_DEAD"]);
            int juv_N      = as<int>(mod["juv_N"]);
              
            NumericVector state_BB     = as<NumericVector>(mod["state_BB"]);
            NumericVector JuvNum       = as<NumericVector>(mod["JuvNum"]);
            NumericVector AduNum       = as<NumericVector>(mod["AduNum"]);
            NumericVector SpawnBio     = as<NumericVector>(mod["SpawnBio"]);
            NumericVector firstMoAdu   = as<NumericVector>(mod["firstMoAdu"]);
            NumericMatrix WageS        = as<NumericMatrix>(mod["WageS"]);
            NumericMatrix NageS        = as<NumericMatrix>(mod["NageS"]);
            DataFrame out_BB           = as<DataFrame>(mod["out_BB"]);
            DataFrame out_SSB          = as<DataFrame>(mod["out_SSB"]);
            DataFrame out_rec          = as<DataFrame>(mod["out_rec"]);
         
            for (y = StartYear; y < EndYear; y++){                                  
            for (m = 0; m < STEPS_PER_YEAR; m++){
                     
           // dd is index for saving monthly output
  			      dd = y * STEPS_PER_YEAR + m;
           
           // save current state to output matrix
           // For non-split pools, SSB is output output as the same as B.
					 // For adult split pools, it is overwritten with "actual" SSB, 
					 // for juvs it is 0. 
              for (sp = 1; sp <= NUM_LIVING + NUM_DEAD; sp++){ 
                  NumericVector sp_out_BB  = out_BB[sp];
                  NumericVector sp_out_SSB = out_SSB[sp];
                  NumericVector sp_out_rec = out_rec[sp];

                  sp_out_BB[dd]  = state_BB[sp];
                  //if(sp = 34)  Rcout << sp_out_BB[dd] << std::endl;
                  sp_out_SSB[dd] = state_BB[sp];
                  sp_out_rec[dd] = state_BB[sp];
              }
              for (i = 1; i <= juv_N; i++){
                  NumericVector jv_out_SSB = out_SSB[JuvNum[i]];
                  NumericVector ad_out_SSB = out_SSB[AduNum[i]];
                  NumericVector ad_out_rec = out_rec[AduNum[i]];
                  
                  jv_out_SSB[dd] = 0.0;
									ad_out_SSB[dd] = SpawnBio[i];
									ad_out_rec[dd] = NageS(firstMoAdu[i], i) * WageS(firstMoAdu[i], i);
              }
}
}
return 0;
            }')        

ecosim_test(test.sim, 0, 3)

ecosim_test(GOA.sim, 0, 100)
as.data.table(GOA.sim$out_BB[34])
