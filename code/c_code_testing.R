#Testing
if(Sys.info()['sysname']=="Windows"){
  data.dir <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\data\\"
  out.dir  <- "C:\\Users\\Sean.Lucey\\Desktop\\Rpath\\outputs\\"
}

if(Sys.info()['sysname']=="Linux"){
  r.dir    <- "/home/slucey/slucey/Rpath/Rpath/R"
  data.dir <- "/home/slucey/slucey/Rpath/data"
  out.dir  <- "/home/slucey/slucey/Rpath/outputs"
}
#source RcppExports.R
source(file.path(r.dir, 'RcppExports.R'))

#Constants from header
MAX_MONTHS_STANZA <-  400      # Max number age pools within stanza calcs (months)
LO_DISCARD        <- 1e-4     # Discard Limit for sense
HI_DISCARD        <- 1e+4     # Discard Limit for sense
STEPS_PER_YEAR    <- 12                   # Steps per year between stanzas ("months")
DELTA_T           <- 0.083333333333333333 # 1/12 for steps per month
SORWT             <- 0.5                  # Fast Equilibrium smoother
EPSILON           <- 1E-8                 # Test threshold for "too close to zero"
BIGNUM            <- 1E+8 

load(file.path(data.dir, 'REco_params.RData'))

library('Rpath'); library('data.table')

REco <- rpath(REco.params, 'R Ecosystem')

REco.sim <- rsim.scenario(REco, REco.params, 100)

params  <- REco.sim$params
state   <- REco.sim$start_state
fishing <- REco.sim$fishing
forcing <- REco.sim$forcing
stanzas <- REco.sim$stanzas

NoIntegrate <- params$NoIntegrate

SplitSetPred(stanzas, state)

dyt <-  deriv_vector(params, state, forcing, fishing, stanzas, 0, 0, 0)

library(Rcpp)

cppFunction('
NumericVector test(int NUM_GROUPS){
  NumericVector dyt0(NUM_GROUPS);
  return(dyt0);
}
')
  
cppFunction('
int test2(int NUM_GROUPS, int y, int m){
  int i;
  NumericVector zero(NUM_GROUPS);
  NumericVector one(NUM_GROUPS);
  for(i = 0; i < NUM_GROUPS; i++){
    one[i] = 1;
  }
  Rcpp::Rcout << one << std::endl;

  return(0);
}
')

dd = 1;
StartYear <- 0
EndYear <- 3
#MAIN LOOP STARTS HERE
# ASSUMES STEPS_PER_MONTH will always be 1.0, took out divisions     
for (y in StartYear:EndYear){                                  
  # Monthly loop                                     
  for (m in 1:12){   
    dd = y * 12 + m; # dd is index for monthly output                
    # Load old state and old derivative
    old_BB    = state$BB
    old_Ftime = state$Ftime         
    dydt0     = dyt$DerivT
    
    # Calculate new derivative    
    dyt   = deriv_vector(params, state, forcing, fishing, stanzas, y, m, 0);
    dydt1     = dyt$DerivT
    
    # Now Update the new State Biomass using Adams-Basforth
    new_BB = 
      ifelse( NoIntegrate == 0,
              (1.0 - SORWT) * dyt$biomeq + SORWT * dyt$old_BB,
              ifelse( NoIntegrate > 0,
                      old_BB + (DELTA_T / 2.0) * (3.0 * dydt1 - dydt0),
                      old_BB)); 
    
    // Then Update Foraging Time 
    NumericVector pd = ifelse(NoIntegrate < 0, stanzaPred, old_BB);
    NumericVector new_Ftime = ifelse((FoodGain > 0) & (pd > 0),
                                     0.1 + 0.9 * old_Ftime * ((1.0 - FtimeAdj) + FtimeAdj * FtimeQBOpt / 
                                                                (FoodGain / pd)),
                                     old_Ftime);
    
    // Monthly Stanza (split pool) update
    if(Nsplit > 0){
      SplitUpdate(stanzas, state, forcing, dyt, y, m + 1);
      SplitSetPred(stanzas, state); 
    }
    new_BB = ifelse(NoIntegrate < 0, as<NumericVector>(state["BB"]), new_BB);
    
    // Calculate catch assuming fixed Frate and exponential biomass change.
    // kya 9/10/15 - replaced with linear, diff on monthly scale is minor
    NumericVector new_CC = (DELTA_T * FishingLoss / old_BB) * 
      (new_BB + old_BB) / 2.0;
