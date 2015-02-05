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

GOA.sim <- ecosim_init(GOA, YEARS = 100, juvfile)

sourceCpp(paste(r.dir, 'test_ecosim.cpp', sep = ''))

deriv_master(GOA.sim, 0, 0, 0)

#test pointers
test <- list(x = data.frame(x = 1:10, y = rep(2, 10)), z = rep(0, 10))

cppFunction('NumericMatrix RmatT (DataFrame x){
              //Tranposes the column and rows
              int nc = x.nrows();
              int nr = x.size();
              NumericMatrix y(Dimension(nr, nc));
              
              for( int i = 0; i < x.size(); i++){
                NumericVector xcol = x[i];
                NumericMatrix::Row yy = y( i, _);
                yy = xcol;
              }
              
              return y;
            }')

cppFunction('int cpptest(List test){
            DataFrame x = as<DataFrame>(test["x"]);
            NumericVector z = as<NumericVector>(test["z"]);
            
            int n = x.nrows();
            for( int i = 0; i < n; i++){
              z[i] = x(i, 1) * x(i, 2);
              x[i, 2]++;
              }
            return 0;
            }')


cpptest(test)
test
cpptest(test)





cppFunction('NumericMatrix DFtoA (DataFrame x){
              int nr = x.nrows();
              int nc = x.size();
              NumericMatrix y(Dimension(nr, nc));
              for( int i = 0; i < x.size(); i++){
                NumericVector xcol = x[i];
                NumericMatrix::Column yy = y( _, i);
                yy = xcol;
              }
              return y;
              }')

cppFunction('NumericMatrix RmatT (DataFrame x){
              //Tranposes the column and rows
              int nr = x.nrows();
              int nc = x.size();
              NumericMatrix y(Dimension(nc, nr));
              for( int i = 0; i < x.size(); i++){
                NumericVector xcol = x[i];
                NumericMatrix::Row yy = y( i, _);
                yy = xcol;
              }
              return y;
              }')

cppFunction('List cppTest(List mod){
            #define STEPS_PER_YEAR 12
            NumericVector Q;
            int NumPredPreyLinks           = as<int>(mod["NumPredPreyLinks"]);
            NumericVector preyYY           = as<NumericVector>(mod["preyYY"]);
            NumericVector PreyFrom         = as<NumericVector>(mod["PreyFrom"]);
            NumericVector PreyTo           = as<NumericVector>(mod["PreyTo"]);
            NumericVector HandleSwitch     = as<NumericVector>(mod["HandleSwitch"]);
            NumericVector COUPLED          = as<NumericVector>(mod["COUPLED"]);

            for (int links=1; links<=NumPredPreyLinks; links++){
   	          int prey = PreyFrom[links];
 		          int pred = PreyTo[links];
              Q[links] = pow(preyYY[prey], COUPLED * HandleSwitch[links]);
            }
            List out = List::create(Named("predYY")         = predYY,
                                    Named("Q")              = Q);
            return out;
            }')

cppFunction('int deriv_master(List mod, int y, int m, int d){

            int STEPS_PER_YEAR = 12;
            int sp, i;

            // Parse out List mod
            int NUM_LIVING                 = as<int>(mod["NUM_LIVING"]);
            int NUM_DEAD                   = as<int>(mod["NUM_DEAD"]);
            int juv_N                      = as<int>(mod["juv_N"]);
            NumericVector spnum            = as<NumericVector>(mod["spnum"]);
            NumericVector B_BaseRef        = as<NumericVector>(mod["B_BaseRef"]);
            NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
            NumericVector state_Ftime      = as<NumericVector>(mod["state_Ftime"]);
            NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
            NumericVector stanzaBasePred   = as<NumericVector>(mod["stanzaBasePred"]);
            NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
            NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
            NumericVector preyYY           = as<NumericVector>(mod["preyYY"]);
            NumericVector predYY           = as<NumericVector>(mod["predYY"]);
            NumericMatrix force_bysearch   = as<NumericMatrix>(mod["force_bysearch"]);

            
            //  Set YY = B/B(ref) for functional response; note that foraging time is
            //  used here as a biomass modifier before the main functional response  
            for (sp = 1; sp <= NUM_LIVING + NUM_DEAD; sp++){
              preyYY[sp] = state_Ftime[sp] * state_BB[sp] / B_BaseRef[sp];
              predYY[sp] = state_Ftime[sp] * state_BB[sp] / B_BaseRef[sp];
              }    
            // The sun always has a biomass of 1 for Primary Production
            preyYY[0] = 1.0;  predYY[0] = 1.0;
            
            // Set functional response biomass for juvenile and adult groups (including foraging time) 
            for (i=1; i<=juv_N; i++){
              if (stanzaBasePred[JuvNum[i]]>0){
                predYY[JuvNum[i]] = state_Ftime[JuvNum[i]] * 
                stanzaPred[JuvNum[i]]/stanzaBasePred[JuvNum[i]];
                predYY[AduNum[i]] = state_Ftime[AduNum[i]] * 
                stanzaPred[AduNum[i]]/stanzaBasePred[AduNum[i]];
              }
            }
            // add "mediation by search rate" KYA 7/8/08
            for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
              predYY[sp] *= force_bysearch(y*STEPS_PER_YEAR+m, sp); 
            }
            }')

deriv_master(GOA.sim, 0, 0, 0)




test <- cpptest(GOA.sim)




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
summary(run0$out_BB)
