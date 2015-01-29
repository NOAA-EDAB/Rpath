#Rpath_testing.R
#Testing of new Rpath functions
#SML
#12/14

#User parameters

windows <- F
if(windows == T){
  r.dir    <- "L:\\Rpath\\code\\"
  data.dir <- "L:\\Rpath\\data\\"
  out.dir  <- "L:\\Rpath\\outputs\\"
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


cppFunction('List cpptest(List x){
            if (!x.inherits("Rpath.sim")) stop("Input must be a Rpath model");            
                       
            int year                       = as<int>(x["YEARS"]);
            CharacterVector spname         = as<Rcpp::CharacterVector>(x["spname"]);
            double *preyYY                 = as<double>(x["preyYY"]);
            List out = List::create(Named("Years")      = year,
                                    Named("Groups")     = spname,
                                    Named("preyYY")     = preyYY);
            
            out.attr("class") = "Rsim";
            return out;
            }')

cppFunction('List cpptest2 (List x){
            DataFrame test     = as<DataFrame>(x["force_bysearch"]);
            const int c = test.size(); 
            const int r = (as<int>(x["YEARS"]) * 12) + 1;
            NumericMatrix test2 = as<NumericMatrix>(test);
//            NumericMatrix test  = as<NumericMatrix>(x["NageS"]);
//            NumericVector test2 = as<NumericVector>(test);
            
//            test2.attr("dim") = Dimension(1, test2.size());
//            List out = List::create(Named("NageS")  = test,
//                                    Named("NageSv") = test2,
//                                    Named("WageS")  = test3);
            return List::create(c, r, test, test2);
}')

cppFunction('NumericVector DFtoV(DataFrame x){
            int nRows = x.nrows();
            NumericVector y(nRows * x.size());
            for( int i = 0; i < x.size(); i++){
              NumericVector xcol = x[i];
              for( int j = 0; j < nRows; j++){   
                y[(i * nRows) + j] = xcol[j];
              }
            }
            return y;
            }')

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

cppFunction('List powtest(NumericVector x, NumericVector y){
            int n = x.size();
            NumericVector powerout(n);
            for (int i = 0; i < x.size(); i++){
              powerout[i] = pow(x[i], y[i]);
            }
            List out = List::create(Named("size") = x.size(),
                                    Named("result") = powerout);
            return out;
            }')

cppTest(GOA.sim$predYY, RmatT(GOA.sim$force_bysearch), GOA.sim$NUM_LIVING, GOA.sim$NUM_DEAD)



test <- cpptest(GOA.sim)

sourceCpp(paste(r.dir, 'test_ecosim.cpp', sep = ''))


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
