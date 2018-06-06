
#include "ecosim.h"
                       
//################################################################----------
// Runge-Kutta 4th order method for integrating Ecosim equations
// Currently does not contain aged-structured species.
// [[Rcpp::export]] 
List rk4_run (List params, List instate, List forcing, List fishing, List stanzas,
                 int StartYear, int EndYear){

   int y, m, dd, t; 
// Input rates are in units of years or years^-1.  Integration is wri tten so
// that integration timesteps always line up with months, for data reasons.
// STEPS_PER_YEAR should be 12 (for months), and STEPS_PER_MONTH sets the
// rk4 integration timestep.  So effective integration timestep with respect
// to input rates (years) is 1/(12*STEPS_PER_MONTH).
   const int STEPS_PER_MONTH = as<int>(params["RK4_STEPS"]);
   const double hh           = DELTA_T/(double)STEPS_PER_MONTH;

// Get some basic needed numbers from the params List
   const int NUM_BIO = as<int>(params["NUM_LIVING"]) + as<int>(params["NUM_DEAD"]);
   const int NumPredPreyLinks = as<int>(params["NumPredPreyLinks"]);
   const int NumFishingLinks  = as<int>(params["NumFishingLinks"]);
   
// Switches for run modes
   const int BURN_YEARS = as<int>(params["BURN_YEARS"]);
   int CRASH_YEAR = -1;
   int MEASURE_MONTH = 5;   
   
// Flag for group-sepcific Integration method (NoIntegrate=0 means Fast Eq)   
   const NumericVector NoIntegrate = as<NumericVector>(params["NoIntegrate"]);
   
// Parameters needed directly for foraging time adjustment
   const NumericVector B_BaseRef        = as<NumericVector>(params["B_BaseRef"]);
   const NumericVector FtimeAdj         = as<NumericVector>(params["FtimeAdj"]);
   const NumericVector FtimeQBOpt       = as<NumericVector>(params["FtimeQBOpt"]);
// FtimeAdj is monthly unit so adjust for sub-monthly integration
   const NumericVector FtimeStep     = FtimeAdj/STEPS_PER_MONTH;

// Number of split groups
   const int Nsplit = as<int>(stanzas["Nsplit"]);

// Parameter need to track catch by Gear
   const NumericVector FishFrom = as<NumericVector>(params["FishFrom"]);
   
// Monthly output matrices                     
   NumericMatrix out_BB(EndYear*12, NUM_BIO+1);           
   NumericMatrix out_CC(EndYear*12, NUM_BIO+1);          
   NumericMatrix out_SSB(EndYear*12, NUM_BIO+1);        
   NumericMatrix out_rec(EndYear*12, NUM_BIO+1);
   NumericMatrix out_Gear_CC(EndYear*12, NumFishingLinks+1);
// Annual output matrices
   NumericMatrix annual_CC(EndYear, NUM_BIO+1);
   NumericMatrix annual_BB(EndYear, NUM_BIO+1);
   NumericMatrix annual_QB(EndYear, NUM_BIO+1);
   NumericMatrix annual_Qlink(EndYear, NumPredPreyLinks+1);   
// Accumulator for monthly catch values
   NumericVector cum_CC(NUM_BIO+1);
   NumericVector cum_Gear_CC(NumFishingLinks +1);

//SML
// Update sums of split groups to total biomass for derivative calcs
   if(Nsplit > 0){ 
     SplitSetPred(stanzas, instate);
   }
//SML

// Load state, set some initial values.  Make sure state is COPY, not pointer   
   List state = clone(instate);
   dd =  StartYear * STEPS_PER_YEAR;  // dd is monthly index for data storage

// KYA 6/12/17 an initial derivative call just to declare deriv in right scope
   List dyt = deriv_vector(params,state,forcing,fishing,stanzas,1,0,0);
      NumericVector FoodGain = as<NumericVector>(dyt["FoodGain"]);
      NumericVector Qlink    = as<NumericVector>(dyt["Qlink"]);   

      // MAIN LOOP STARTS HERE with years loop
   for (y = StartYear; y <= EndYear; y++){
   if (y<1){stop("RK Year can't be less than 1");}
   // Monthly loop                                     
      for (m = 0; m < STEPS_PER_YEAR; m++){
         cum_CC = NumericVector(NUM_BIO+1);  // monthly catch to accumulate   
         cum_Gear_CC = NumericVector(NumFishingLinks+1);
         dd = (y-1) * STEPS_PER_YEAR + m;  
         // Sub-monthly integration loop
            for (t=0; t< STEPS_PER_MONTH; t++){
               double tt = (double)t*hh;  // sub monthly timestep in years     
            // Load state vars (note: this creates pointers, not new vals)
               NumericVector old_BB    = as<NumericVector>(state["BB"]);
               NumericVector old_Ftime = as<NumericVector>(state["Ftime"]);         

            // Calculate base derivative and RK-4 derivative steps (overwrites YY)   
               List YY = state;
               List k1 = deriv_vector(params,YY,forcing,fishing,stanzas,y,m,tt);
               NumericVector kk1 = as<NumericVector>(k1["DerivT"]);  
      
               YY["BB"] = old_BB + 0.5*kk1*hh;
               List k2 = deriv_vector(params,YY,forcing,fishing,stanzas,y,m,tt + 0.5*hh);
               NumericVector kk2 = as<NumericVector>(k2["DerivT"]);  

               YY["BB"] = old_BB + 0.5*kk2*hh;
               List k3 = deriv_vector(params,YY,forcing,fishing,stanzas,y,m,tt + 0.5*hh);
               NumericVector kk3 = as<NumericVector>(k3["DerivT"]);

               YY["BB"] = old_BB + kk3*hh; 
               List k4 = deriv_vector(params,YY,forcing,fishing,stanzas,y,m,tt + hh);
               NumericVector kk4 = as<NumericVector>(k4["DerivT"]);  

            // Take an rk4 step          
               NumericVector new_BB = old_BB + hh*(kk1 + 2*kk2 + 2*kk3 + kk4)/6.0;   

           // Update Foraging time state variable
           // pd term is used to indicate differrent values used for 
           // age-structured species, defaults to BB for non-aged structure
              NumericVector pd = old_BB;
              FoodGain         = as<NumericVector>(k1["FoodGain"]);
              Qlink            = as<NumericVector>(k1["Qlink"]);
              NumericVector new_Ftime =  ifelse((FoodGain>0)&(pd>0),
                 0.1 + 0.9*old_Ftime* 
                 ((1.0-FtimeStep) + FtimeStep*FtimeQBOpt/(FoodGain/pd)),
                 old_Ftime);

          // Accumulate Catch (small timestep, so linear average)
             NumericVector FishingLoss = as<NumericVector>(k1["FishingLoss"]);
             cum_CC += (hh * FishingLoss/old_BB) * (new_BB+old_BB)/2.0;
             
          // Track catch by gear
             NumericVector old_BB_flink = as<NumericVector>(old_BB[FishFrom]);
             NumericVector new_BB_flink = as<NumericVector>(new_BB[FishFrom]);
             NumericVector GearCatch = as<NumericVector>(k1["GearCatch"]);
             cum_Gear_CC += (hh * GearCatch / old_BB_flink) * (new_BB_flink + old_BB_flink)/2.0;
             
          // Set state to new values, including min/max traps
             state["BB"]    = pmax(pmin(new_BB, B_BaseRef*BIGNUM), B_BaseRef*EPSILON);
             state["Ftime"] = pmin(new_Ftime, 2.0);      
          }// end of sub-monthly (t-indexed) loop

      // Insert Monthly Stanza (split pool) update here
      //SML
         // Calculate new derivative    
        List dyt = deriv_vector(params, state, forcing, fishing, stanzas, y, m, 0);
        if(Nsplit > 0){ 
          SplitUpdate(stanzas, state, forcing, dyt, y, m + 1);
          SplitSetPred(stanzas, state);
        }
      
      // Make a copy of the current state for bounds testing
         NumericVector cur_BB = as<NumericVector>(state["BB"]);
        
        // KYA 8/9/17 one of the NA or NaN flags is reading back as a negative integer (-2^32)
        // Not sure why.  This sets any negative biomass (assuming this means NaN) to NA_REAL
        cur_BB = ifelse((cur_BB<0),NA_REAL,cur_BB);
        
        // If the run is during the "burn-in" years, and biomass goes
        // into the discard range, set flag to exit the loop.  Should set "bad"
        // biomass values to NA                 
        if (y < BURN_YEARS){ cur_BB = ifelse((cur_BB<B_BaseRef*LO_DISCARD)|
            (cur_BB>B_BaseRef*HI_DISCARD),
            NA_REAL,cur_BB);
        }
        
        // If biomass goes crazy or hits NA, exit loop with crash signal.  Note it 
        // should still write the NA or INF values back to the output.
        if ( any(is_na(cur_BB)) | any(is_infinite(cur_BB)) | any(is_nan(cur_BB)) )  {
          CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;
        }
        
        
      // Write to monthly output matricies (vector write)     				          									                    
         out_BB( dd, _) = cur_BB;
         out_SSB(dd, _) = cur_BB;
         out_rec(dd, _) = cur_BB;
         out_CC( dd, _) = cum_CC;
         out_Gear_CC(dd, _) = cum_Gear_CC;
         annual_CC(y-1, _) = annual_CC(y-1, _) + cum_CC;
         if (m==MEASURE_MONTH){
           annual_BB(y-1, _)    = cur_BB;
           annual_QB(y-1, _)    = FoodGain/cur_BB;
           annual_Qlink(y-1, _) = Qlink;
         }      
    }  // End of main months loop
    
  }// End of years loop
  
// Write Last timestep (note:  SSB and rec only used for age-structure species) 
   //out_BB( dd+1, _) = as<NumericVector>(state["BB"]);
   //out_SSB(dd+1, _) = as<NumericVector>(state["BB"]);
   //out_rec(dd+1, _) = as<NumericVector>(state["BB"]);
   //out_CC( dd+1, _) = out_CC( dd, _); 
  
// Create Rcpp list to output  
List outdat = List::create(
  _["out_BB"]=out_BB,
  _["out_CC"]=out_CC,
  _["out_Gear_CC"]=out_Gear_CC,
  _["annual_CC"]=annual_CC,
  _["annual_BB"]=annual_BB,
  _["annual_QB"]=annual_QB,
  _["annual_Qlink"]=annual_Qlink,
  _["end_state"]=state,
  _["crash_year"]=CRASH_YEAR);
  
// Return is an Rcpp List
   return(outdat);
} 

//-----#################################################################----
// Adams_Basforth two-step method for integrating ecosim equations.
// Fixed monthly timestep with "fast integration" Walters et al. '97 method
// for faster-turnover species.
// Currently does not contain aged-structured species.
// [[Rcpp::export]] 
List Adams_run (List params, List instate, List forcing, List fishing, List stanzas,
                 int StartYear, int EndYear, List InitDeriv){
     
int y, m, dd; 

// Get some basic needed numbers from the params List
   const int NUM_BIO = as<int>(params["NUM_LIVING"]) + as<int>(params["NUM_DEAD"]);
   const int NumPredPreyLinks          = as<int>(params["NumPredPreyLinks"]);
   const int NumFishingLinks  = as<int>(params["NumFishingLinks"]);

// Forcing Biomass
   //const int NumForcedBio        = as<int>(fitting["NumForcedBio"]);
   //const IntegerVector BforceNum = as<IntegerVector>(fitting["BforceNum"]);
   // KYA TODO 11/1/2017 put force_bybio in RK!
   NumericMatrix force_bybio       = as<NumericMatrix>(forcing["bybio"]);

// Switches for run modes
   const int BURN_YEARS = as<int>(params["BURN_YEARS"]);
   int CRASH_YEAR = -1;
   int MEASURE_MONTH = 5;
   
// Flag for group-sepcific Integration method (NoIntegrate=0 means Fast Eq)   
   const NumericVector NoIntegrate = as<NumericVector>(params["NoIntegrate"]);

// Parameters needed directly for foraging time adjustment
   const NumericVector B_BaseRef  = as<NumericVector>(params["B_BaseRef"]);
   const NumericVector FtimeAdj   = as<NumericVector>(params["FtimeAdj"]);
   const NumericVector FtimeQBOpt = as<NumericVector>(params["FtimeQBOpt"]);
   
// Parameters from stanzas
   const int Nsplit         = as<int>(stanzas["Nsplit"]);
   NumericVector stanzaPred = as<NumericVector>(stanzas["stanzaPred"]);
   
// Parameter need to track catch by Gear
   const NumericVector FishFrom = as<NumericVector>(params["FishFrom"]);
   
// Monthly output matrices                     
   NumericMatrix out_BB( EndYear * 12, NUM_BIO + 1);           
   NumericMatrix out_CC( EndYear * 12, NUM_BIO + 1);          
   NumericMatrix out_SSB(EndYear * 12, NUM_BIO + 1);        
   NumericMatrix out_rec(EndYear * 12, NUM_BIO + 1);
   NumericMatrix out_Gear_CC(EndYear*12, NumFishingLinks+1);
// Annual output matrices
   NumericMatrix annual_CC(EndYear, NUM_BIO+1);
   NumericMatrix annual_BB(EndYear, NUM_BIO+1);
   NumericMatrix annual_QB(EndYear, NUM_BIO+1);
   NumericMatrix annual_Qlink(EndYear, NumPredPreyLinks+1); 

// Load state and call initial derivative (todo: allow other start times)
// Use Clone to make sure state/stanzas are copies of instate/instanzas, not pointers   
   List state = clone(instate);
   //List stanzas = clone(instanzas);
   
   // Update sums of split groups to total biomass for derivative calcs
   if(Nsplit > 0){
     SplitSetPred(stanzas, state); 
   } 
   
   // Use the initial derivative calculated outside of function
      List dyt = InitDeriv;

   dd = StartYear * STEPS_PER_YEAR;

// MAIN LOOP STARTS HERE
// ASSUMES STEPS_PER_MONTH will always be 1.0, took out divisions     
   for (y = StartYear; y <= EndYear; y++){
      if (y<1){stop("Adams Year can't be less than 1");}
   // Monthly loop                                     
      for (m = 0; m < STEPS_PER_YEAR; m++){   
				 dd = (y-1) * STEPS_PER_YEAR + m; // dd is index for monthly output                
      // Load old state and old derivative
         NumericVector old_BB    = as<NumericVector>(state["BB"]);
         NumericVector old_Ftime = as<NumericVector>(state["Ftime"]);         
 	       NumericVector dydt0     = as<NumericVector>(dyt["DerivT"]);
         
      // Calculate new derivative    
 	       dyt   = deriv_vector(params, state, forcing, fishing, stanzas, y, m, 0);
      
      // Extract needed parts of the derivative
         NumericVector dydt1       = as<NumericVector>(dyt["DerivT"]); 
 				 NumericVector FoodGain    = as<NumericVector>(dyt["FoodGain"]);					
         NumericVector biomeq      = as<NumericVector>(dyt["biomeq"]);
         NumericVector FishingLoss = as<NumericVector>(dyt["FishingLoss"]);  
         NumericVector Qlink       = as<NumericVector>(dyt["Qlink"]);
               
      // Now Update the new State Biomass using Adams-Basforth
         NumericVector new_BB = 
                       ifelse( NoIntegrate == 0,
                         (1.0 - SORWT) * biomeq + SORWT * old_BB,
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
        
        // Track catch by gear
        NumericVector old_BB_flink = as<NumericVector>(old_BB[FishFrom]);
        NumericVector new_BB_flink = as<NumericVector>(new_BB[FishFrom]);
        NumericVector GearCatch = as<NumericVector>(dyt["GearCatch"]);
        NumericVector new_Gear_CC = (DELTA_T * GearCatch / old_BB_flink) * 
                                    (new_BB_flink + old_BB_flink)/2.0;
        
        // NumericVector new_CC = 
        //               ifelse( new_BB==old_BB,
        //                 FishingLoss*DELTA_T,
        //                 (FishingLoss*DELTA_T/old_BB) *
        //                   (new_BB-old_BB)/log(new_BB/old_BB) 
        //               );       

     // KYA 9/13/17 - noticed that cur_BB is never copied back into state["BB"] after
     // bounds testing.  Error?  Moved state["BB"] setting, hopefully nothing breaks...
     // Set state to new values, including min/max traps
        //state["BB"]    = pmax(pmin(new_BB, B_BaseRef * BIGNUM), B_BaseRef * EPSILON);
        //state["Ftime"] = pmin(new_Ftime, 2.0);    
                                                                 										                         
     // Make a copy of the current state for bounds testing
        NumericVector cur_BB = pmax(pmin(new_BB, B_BaseRef * BIGNUM), B_BaseRef * EPSILON); // as<NumericVector>(state["BB"]);

        NumericVector bforce = force_bybio((y-1) * STEPS_PER_YEAR + m, _);
        cur_BB = ifelse(bforce>B_BaseRef * EPSILON, bforce, cur_BB);
        
     // insert forced biomass levels    
        //for (i=0; i<NumForcedBio; i++){
         // sp = BforceNum[i]; /* if (NoIntegrate[sp]>=0){cur_BB[sp]=BforceVal(sp,m);} */
        //}            
        
     // KYA 8/9/17 one of the NA or NaN flags is reading back as a negative integer (-2^32)
     // Not sure why.  This sets any negative biomass (assuming this means NaN) to NA_REAL
        cur_BB = ifelse((cur_BB<0),NA_REAL,cur_BB);

     // If the run is during the "burn-in" years, and biomass goes
     // into the discard range, set flag to exit the loop.  Should set "bad"
     // biomass values to NA                 
        if (y < BURN_YEARS){ cur_BB = ifelse((cur_BB<B_BaseRef*LO_DISCARD)|
                                             (cur_BB>B_BaseRef*HI_DISCARD),
                                      NA_REAL,cur_BB);
        }
      
    // If biomass goes crazy or hits NA, exit loop with crash signal.  Note it 
    // should still write the NA or INF values back to the output.
    //NOJUV make sure crash tests work for juveniles.
    
           if ( any(is_na(cur_BB)) | any(is_infinite(cur_BB)) | any(is_nan(cur_BB)) )  {
          CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;
       }
  
    // KYA 9/13/17 - Now copy cur_BB into state    
       state["BB"]    = cur_BB;
       state["Ftime"] = pmin(new_Ftime, 2.0);
  
     // Write to output matricies     				          									                    
        out_BB( dd, _) = old_BB;
        out_SSB(dd, _) = old_BB;
        out_rec(dd, _) = old_BB;
        out_CC( dd, _) = new_CC;
        out_Gear_CC(dd, _) = new_Gear_CC;
        annual_CC(y-1, _) = annual_CC(y-1, _) + new_CC;
        if (m==MEASURE_MONTH){
          annual_BB(y-1, _)    = old_BB;
          annual_QB(y-1, _)    = FoodGain/old_BB;
          annual_Qlink(y-1, _) = Qlink;
        }
        
     //NOJUV    for (i = 1; i <= juv_N; i++){
     //NOJUV    out_SSB(dd, JuvNum[i]) = 0.0;
     //NOJUV 		out_SSB(dd, AduNum[i]) = SpawnBio[i];
     //NOJUV 		out_rec(dd, AduNum[i]) = NageS(firstMoAdu[i], i) * WageS(firstMoAdu[i], i);

     }  // End of main months loop
     
   }// End of years loop

// Write Last timestep 
   //out_BB( dd+1, _) = as<NumericVector>(state["BB"]);
   //out_SSB(dd+1, _) = as<NumericVector>(state["BB"]);
   //out_rec(dd+1, _) = as<NumericVector>(state["BB"]);
   //out_CC( dd+1, _) = out_CC( dd, _); // the "next" time interval
   //out_Gear_CC(dd+1, _) = out_Gear_CC(dd, _);        
//NOJUV         for (i = 1; i <= juv_N; i++){
//NOJUV              out_SSB(dd, JuvNum[i]) = 0.0;
//NOJUV						  out_SSB(dd, AduNum[i]) = SpawnBio[i];
//NOJUV						  out_rec(dd, AduNum[i]) = NageS(firstMoAdu[i], i) * WageS(firstMoAdu[i], i);
//NOJUV         }
 
// Create Rcpp list to output  
   List outdat = List::create(
     _["out_BB"]=out_BB,
     _["out_CC"]=out_CC,
     _["out_Gear_CC"]=out_Gear_CC,
     _["annual_CC"]=annual_CC,
     _["annual_BB"]=annual_BB,
     _["annual_QB"]=annual_QB,
     _["annual_Qlink"]=annual_Qlink,     
     _["end_state"]=state,
     _["crash_year"]=CRASH_YEAR,
     _["dyt"]=dyt);
  
// Return is an Rcpp List
   return(outdat);
                     
} 

//################################################################----------
// Main derivative (dB/dt) calculations for ecosim model, using Rcpp vector
// package.
// [[Rcpp::export]] 
List deriv_vector(List params, List state, List forcing, List fishing, List stanzas, 
                  int inyear, int m, double tt){

int sp, links, prey, pred, gr, egr, dest, isp, ist, ieco;
   
   //Rcout << inyear <<" "<< m << std::endl;
   if (inyear<1){stop("Derivative Year can't be less than 1");}
   
// forcing time index (in months)
   // KYA 11/1/17 - added offset to deal with 0 vs 1 indexing in forcing files
   // (matters when trying to line up with "acutal" years)
   const int y  = inyear - 1;
   const int dd = y*STEPS_PER_YEAR+m;

// Base model size - number of groups and number of links (flows) by type
   const int NUM_GROUPS                = as<int>(params["NUM_GROUPS"]);
   const int NUM_LIVING                = as<int>(params["NUM_LIVING"]);
   const int NUM_DEAD                  = as<int>(params["NUM_DEAD"]);
   const int NumPredPreyLinks          = as<int>(params["NumPredPreyLinks"]);
   const int NumFishingLinks           = as<int>(params["NumFishingLinks"]);
   const int NumDetLinks               = as<int>(params["NumDetLinks"]);
   const int COUPLED                   = as<int>(params["COUPLED"]);
 
// NUM_GROUPS length input vectors
   const NumericVector B_BaseRef       = as<NumericVector>(params["B_BaseRef"]);
   const NumericVector MzeroMort       = as<NumericVector>(params["MzeroMort"]);
   const NumericVector UnassimRespFrac = as<NumericVector>(params["UnassimRespFrac"]);
   const NumericVector ActiveRespFrac  = as<NumericVector>(params["ActiveRespFrac"]);
   const NumericVector HandleSelf      = as<NumericVector>(params["HandleSelf"]);
   const NumericVector ScrambleSelf    = as<NumericVector>(params["ScrambleSelf"]);
   const NumericVector fish_Effort     = as<NumericVector>(params["fish_Effort"]);

// NumPredPreyLinks Length vectors
   const IntegerVector PreyFrom        = as<IntegerVector>(params["PreyFrom"]);
   const IntegerVector PreyTo          = as<IntegerVector>(params["PreyTo"]);
   const NumericVector QQ              = as<NumericVector>(params["QQ"]);
   const NumericVector DD              = as<NumericVector>(params["DD"]);
   const NumericVector VV              = as<NumericVector>(params["VV"]);
   const NumericVector HandleSwitch    = as<NumericVector>(params["HandleSwitch"]);
   const NumericVector PredPredWeight  = as<NumericVector>(params["PredPredWeight"]);
   const NumericVector PreyPreyWeight  = as<NumericVector>(params["PreyPreyWeight"]);

// NumFishingLinks lenghted vectors
   const IntegerVector FishFrom        = as<IntegerVector>(params["FishFrom"]);
   const IntegerVector FishThrough     = as<IntegerVector>(params["FishThrough"]);
   const IntegerVector FishTo          = as<IntegerVector>(params["FishTo"]);
   const NumericVector FishQ           = as<NumericVector>(params["FishQ"]);
   const IntegerVector DetFrom         = as<IntegerVector>(params["DetFrom"]);
   const IntegerVector DetTo           = as<IntegerVector>(params["DetTo"]);
   const NumericVector DetFrac         = as<NumericVector>(params["DetFrac"]);

// Age-structured parameters
   const int Nsplit                   = as<int>(stanzas["Nsplit"]);
   const NumericVector Nstanzas       = as<NumericVector>(stanzas["Nstanzas"]);
   const NumericVector stanzaPred     = as<NumericVector>(stanzas["stanzaPred"]);
   const NumericVector stanzaBasePred = as<NumericVector>(stanzas["stanzaBasePred"]);
   NumericMatrix EcopathCode          = as<NumericMatrix>(stanzas["EcopathCode"]);

// State vectors
   const NumericVector state_BB        = as<NumericVector>(state["BB"]);
   const NumericVector state_Ftime     = as<NumericVector>(state["Ftime"]);

//FISHING  NumericVector TerminalF        = as<NumericVector>(params["TerminalF"]);
//FISHING    NumericVector TARGET_BIO       = as<NumericVector>(params["TARGET_BIO"]);
//FISHING    NumericVector TARGET_F         = as<NumericVector>(params["TARGET_F"]);
//FISHING    NumericVector ALPHA            = as<NumericVector>(params["ALPHA"]);

// "Environmental" forcing matrices (dd-indexed month x species)
// SHOULD BE CONST, but no row extraction for CONST (per Rcpp issues wiki)
   NumericMatrix force_byprey     = as<NumericMatrix>(forcing["byprey"]);
   NumericMatrix force_bymort     = as<NumericMatrix>(forcing["bymort"]);
   NumericMatrix force_bysearch   = as<NumericMatrix>(forcing["bysearch"]);
   NumericMatrix force_bymigrate  = as<NumericMatrix>(forcing["bymigrate"]);
   
// Fishing forcing matrices (indexed year x species)  
// SHOULD BE CONST, but no row extraction for CONST (per Rcpp issues wiki)
   NumericMatrix FORCED_FRATE     = as<NumericMatrix>(fishing["FRATE"]);
   NumericMatrix FORCED_CATCH     = as<NumericMatrix>(fishing["CATCH"]);
   NumericMatrix EffortMat        = as<NumericMatrix>(fishing["EFFORT"]); 

// Components of derivative calculated here  
   NumericVector TotGain(NUM_GROUPS+1);       
   NumericVector TotLoss(NUM_GROUPS+1);         
   NumericVector LossPropToB(NUM_GROUPS+1);     
   NumericVector LossPropToQ(NUM_GROUPS+1);     
   NumericVector DerivT(NUM_GROUPS+1);           
   NumericVector biomeq(NUM_GROUPS+1);            
   NumericVector FoodLoss(NUM_GROUPS+1);       
   NumericVector FoodGain(NUM_GROUPS+1);       
   NumericVector UnAssimLoss(NUM_GROUPS+1);    
   NumericVector ActiveRespLoss(NUM_GROUPS+1);    
   NumericVector DetritalGain(NUM_GROUPS+1);   
   NumericVector FishingGain(NUM_GROUPS+1);    
   NumericVector MzeroLoss(NUM_GROUPS+1);
   NumericVector FishingLoss(NUM_GROUPS+1);
   NumericVector DetritalLoss(NUM_GROUPS+1);
   NumericVector FishingThru(NUM_GROUPS+1);
   NumericVector PredSuite(NUM_GROUPS+1);
   NumericVector HandleSuite(NUM_GROUPS+1); 
   NumericVector GearCatch(NumFishingLinks+1);
   NumericVector MigrateLoss(NUM_GROUPS+1);
   
// Set effective biomass for pred/prey response
// default is B/Bref
   NumericVector preyYY = state_Ftime * state_BB/B_BaseRef * force_byprey(dd,_);;
   NumericVector predYY = state_Ftime * state_BB/B_BaseRef * force_bysearch(dd,_);

// Set functional response biomass for juvenile and adult groups (including foraging time) 
   if(Nsplit > 0){
     for (isp = 1; isp <=Nsplit; isp++){
       for(ist = 1; ist <= Nstanzas[isp]; ist++){
         ieco = EcopathCode(isp, ist);
         if (stanzaBasePred[ieco] > 0){
           predYY[ieco] = state_Ftime[ieco] * stanzaPred[ieco] /
                          stanzaBasePred[ieco];
         }
       }
     }
   }
   
// Unroll Biomass Vectors (match pred, prey biomass for all links)
   NumericVector PYY = preyYY[PreyFrom];
   NumericVector PDY = predYY[PreyTo];

// Summed predator and prey suites for joint handling time and/or scramble functional response
   for (links=1; links<=NumPredPreyLinks; links++){
     PredSuite[PreyFrom[links]] += predYY[PreyTo[links]  ] * PredPredWeight[links];
     HandleSuite[PreyTo[links]] += preyYY[PreyFrom[links]] * PreyPreyWeight[links];
   }
   
// Unroll the suites into a NumPredPrey length vector
   NumericVector PdSuite = PredSuite[PreyFrom];
   NumericVector PySuite = HandleSuite[PreyTo];
   NumericVector Hself   = HandleSelf[PreyTo]; 
   NumericVector Sself   = ScrambleSelf[PreyTo];

//   // Main VECTOR CALC to calculate functional response for each predator/prey link
//     // (3) Additive version: primary used and published in Aydin (2004) 
//     // KYA 3/2/2012 setting "COUPLED" to zero means species are density dependent
//     // (based on their own modul) but don't interact otherwise.  This can magically
//     // create and destroy energy in the system but makes them act like a set
//     // of independent surplus production models for comparison purposes 
// NON-VECTOR VERSION:
//     Q =   QQ[links] * predYY[pred] * pow(preyYY[prey], COUPLED * HandleSwitch[links]) *
//       ( DD[links] / ( DD[links] - 1.0 + 
//       pow(HandleSelf[pred] * preyYY[prey]   + 
//       (1. - HandleSelf[pred]) * HandleSuite[pred],
//                                            COUPLED * HandleSwitch[links])) )*
//                                              ( VV[links] / ( VV[links] - 1.0 + 
//                                              ScrambleSelf[pred] * predYY[pred] + 
//                                              (1. - ScrambleSelf[pred]) * PredSuite[prey]) );
// Rcpp VECTOR VERSION
   NumericVector Q1 = 
             QQ * PDY * vpow(PYY, HandleSwitch * COUPLED) *
           ( DD / ( DD-1.0 + vpow((1.-Hself)*PYY + Hself*PySuite, COUPLED*HandleSwitch)) ) *
           ( VV / ( VV-1.0 +      (1.-Sself)*PDY + Sself*PdSuite) );
   Q1[0] = 1.0; // get rid of NaN - moved from KYA's code 6/12/17

// No vector solution here as we need to sum by both links and species 
   for (links=1; links<=NumPredPreyLinks; links++){
      prey = PreyFrom[links];
      pred = PreyTo[links];
   // If model is uncoupled, food loss doesn't change with prey or predator levels.
      if (COUPLED){  FoodLoss[prey]  += Q1[links]; }
      else{  FoodLoss[prey]  += state_BB[prey] * QQ[links]/B_BaseRef[prey]; }
      FoodGain[pred]         += Q1[links];
   }

// By Species Rates 
   UnAssimLoss    = FoodGain  * UnassimRespFrac; 
   ActiveRespLoss = FoodGain  * ActiveRespFrac;  												 
   MzeroLoss      = MzeroMort * state_BB;
   
   NumericVector NetProd = FoodGain - UnAssimLoss - ActiveRespLoss - MzeroLoss - FoodLoss;
   
// FISHING FUNCTIONS (multiple options depending on fishing method)
   
//   // MOST OF THE FOLLOWING FISHING SPECIFICATION METHODS ARE NOT SUPPORTED
//   // BY THE R-CODE, only fishing by effort (for gear) or by F-rate (for
//   // species) is supported at the end.
//   //
//   // BY CURRENT R-CODE.  ONLY    
//   // RFISH for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
//   // RFISH This sets EFFORT by time series of gear-target combinations
//   // RFISH 		   for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
//   // RFISH if -1 is an input value, uses TERMINAL F (last non-negative F) 		
//   //RFISH 		   for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
//   //RFISH 		       if (y+m+d == 0){fish_Effort[gr]=1.0;}
//   //RFISH 		       else           {fish_Effort[gr]=1.0;} // NOTE DEFAULT!  THIS CAN BE CHANGED TO 1.0
//   //RFISH        // Added 7/8/08 for forced effort
//   //RFISH            if (FORCED_EFFORT[gr][y] > -0.001) 
//   //RFISH 					    {fish_Effort[gr]=FORCED_EFFORT[gr][y];}
//   //RFISH 
//   //RFISH 			     if ((FORCED_TARGET[gr]>0) && (FORCED_CATCH[gr][y]>-EPSILON)){
//   //RFISH 			        totQ = 0.0;
//   //RFISH 			        sp   = FORCED_TARGET[gr];
//   //RFISH 			        for (links=1; links<=NumFishingLinks; links++){
//   //RFISH     					    if ((FishingThrough[links] == gr) && 
//   //RFISH 		    					    (FishingFrom[links]) == sp){
//   //RFISH 				    					totQ += FishingQ[links];
//   //RFISH 									}
//   //RFISH 						  }
//   //RFISH 						  fish_Effort[gr] = FORCED_CATCH[gr][y]/ 
//   //RFISH 							                  (totQ * state_BB[sp]);
//   //RFISH 							if (FORCED_CATCH[gr][y] >= state_BB[sp])
//   //RFISH 							   {fish_Effort[gr] = (1.0-EPSILON)*(state_BB[sp])/ 
//   //RFISH 				    			                  (totQ * state_BB[sp]);}
//   //RFISH 					 }
//   //RFISH 					 // By putting F after catch, Frates override absolute catch
//   //RFISH 			     if ((FORCED_FTARGET[gr]>0) && (FORCED_FRATE[gr][y]>-EPSILON)){
//   //RFISH 			        totQ = 0.0;
//   //RFISH 			        sp   = FORCED_FTARGET[gr];
//   //RFISH 			        for (links=1; links<=NumFishingLinks; links++){
//   //RFISH     					    if ((FishingThrough[links] == gr) && 
//   //RFISH 		    					    (FishingFrom[links]) == sp){
//   //RFISH 				    					totQ += FishingQ[links];
//   //RFISH 									}
//   //RFISH 						  }
//   //RFISH 						  fish_Effort[gr] = FORCED_FRATE[gr][y]/totQ;
//   //RFISH 							//if (FORCED_CATCH[gr][y] >= state_BB[sp])
//   //RFISH 							//   {fish_Effort[gr] = (1.0-EPSILON)*(state_BB[sp])/ 
//   //RFISH 				    //			                  (totQ * state_BB[sp]);}
//   //RFISH 					 }					 
//   //RFISH 					 
//   //RFISH 				   //if ((y==0) && (m==0) && (d==0)){
//   //RFISH 					 //    cout << path_species[gr] << " " << FORCED_TARGET[gr] << " " << path_species[sp] << " " 
//   //RFISH 					 //	      << state_BB[sp] << " " << FORCED_CATCH[gr][y] << " "
//   //RFISH 					 //		      << fish_Effort[gr] << endl;
//   //RFISH 					 //} 					 
//   //RFISH 			 }
//   
 double caught;   
   // Apply specified Effort by Gear to catch (using Ecopath-set Q)
      NumericVector EFFORT = (NumericVector)EffortMat(dd,_);
      for (links=1; links<=NumFishingLinks; links++){
 				 prey = FishFrom[links];
 				 gr   = FishThrough[links];
 				 dest = FishTo[links];
         egr  = FishThrough[links] - (NUM_LIVING + NUM_DEAD);
 				 caught =  FishQ[links] * EFFORT[egr] * state_BB[prey]; 
          FishingLoss[prey] += caught;
          FishingThru[gr]   += caught;
          FishingGain[dest] += caught;
          GearCatch[links] = caught;
 		}		
    NumericVector FORCE_F = (NumericVector)FORCED_FRATE(y,_);
    //  Special "CLEAN" fisheries assuming q=1, so specified input is Frate
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
             caught = FORCED_CATCH(y, sp) + FORCE_F[sp] * state_BB[sp];
             // KYA Aug 2011 removed terminal effort option to allow negative fishing pressure 
                // if (caught <= -EPSILON) {caught = TerminalF[sp] * state_BB[sp];}
             // KYA 10/6/17 Added productivity to BB limit for F>1 species (salmon inspired)
                if (caught >= state_BB[sp] + NetProd[sp]){caught = (1.0 - EPSILON) * (state_BB[sp] + NetProd[sp]);}
             FishingLoss[sp] += caught;
             FishingThru[0]  += caught;
             FishingGain[0]  += caught;
             //if(sp==1){Rprintf("%d %g %g %g \n",sp,FORCED_FRATE(y,sp),caught,FishingLoss[sp]);}
             //TerminalF[sp] = caught/state_BB[sp];
        }
            
// KINKED CONTROL RULE - NEEDS INPUT of TARGET BIOMASS and TARGET CATCH
  //FISHING      double RefBio, maxcaught;
  //FISHING      for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
  //FISHING        if (TARGET_BIO[sp] > EPSILON){
  //FISHING          RefBio    = state_BB[sp] / TARGET_BIO[sp];
  //FISHING          maxcaught = TARGET_F[sp] * state_BB[sp];         
  //FISHING          if      (RefBio > 1.0)           {caught = maxcaught;}
  //FISHING          else if (RefBio >= ALPHA[sp]) {caught = maxcaught * (RefBio - ALPHA[sp])/(1.0 - ALPHA[sp]);}
  //FISHING          else                             {caught = 0.0;}
  //FISHING          FishingLoss[sp] += caught;
  //FISHING          FishingThru[0]  += caught;
  //FISHING          FishingGain[0]  += caught;
  //FISHING          TerminalF[sp] = caught/state_BB[sp];
  //FISHING        }          
  //FISHING      }
    
// DETRITUS  - note: check interdetrital flow carefully, have had some issues
// (check by ensuring equlibrium run stays in equilibrium)
   int liv, det;
   double flow;
   for (links=1; links<=NumDetLinks; links++){
      liv  = DetFrom[links];
      det  = DetTo[links];
      flow = DetFrac[links] * (MzeroLoss[liv] + UnAssimLoss[liv]);
      DetritalGain[det] += flow;
      if (liv > NUM_LIVING) {DetritalLoss[liv] += flow; }
   }
   for (sp=NUM_LIVING+1; sp<=NUM_LIVING+NUM_DEAD; sp++){
      MzeroLoss[sp] = 0.0;
   }
    
// Add mortality forcing
   for (int i=1; i<=NUM_DEAD+NUM_LIVING; i++){
     FoodLoss[i]  *= force_bymort(y * STEPS_PER_YEAR + m, i);
     MzeroLoss[i] *= force_bymort(y * STEPS_PER_YEAR + m, i);
   }
   
// Add migration forcing
   MigrateLoss = clone(state_BB);
   for (int i=1; i<=NUM_DEAD+NUM_LIVING; i++){
     MigrateLoss[i]  *= force_bymigrate(y * STEPS_PER_YEAR + m, i);
   }
   
// Sum up derivitive parts (vector sums)
// Override for group 0 (considered "the sun", never changing)        
   TotGain = FoodGain + DetritalGain + FishingGain;      
   LossPropToQ = UnAssimLoss + ActiveRespLoss;
   LossPropToB = FoodLoss + MzeroLoss + FishingLoss + MigrateLoss + DetritalLoss; 
   TotGain[0]     = 0;
   LossPropToB[0] = 0;  
   LossPropToQ[0] = 0;
   TotLoss = LossPropToQ + LossPropToB;      
   biomeq  = TotGain/(TotLoss/state_BB);
   biomeq[0] = 1.0;
   DerivT  = TotGain - TotLoss; 

// Rcpp List structure to return
// KYA 6/20/17 Rcpp bug (known) is max 18 items on List::create
// had to add Q1 (qlink), so removed LossPropToQ (used nowhere?)
// SML 8/8/17 - actually don't need most of these...commenting
// out to track catch by gear
   List deriv = List::create(
     //_["preyYY"]=preyYY,
     //_["predYY"]=predYY,
     _["TotGain"]=TotGain,
     _["TotLoss"]=TotLoss,
     _["DerivT"]=DerivT,
     _["biomeq"]=biomeq,
     _["LossPropToB"]=LossPropToB, 
     //_["LossPropToQ"]=LossPropToQ,                           
     _["FoodLoss"]=FoodLoss,
     _["FoodGain"]=FoodGain,
     _["UnAssimLoss"]=UnAssimLoss,
     _["ActiveRespLoss"]=ActiveRespLoss,
     _["DetritalGain"]=DetritalGain,
     _["FishingGain"]=FishingGain,
     _["MzeroLoss"]=MzeroLoss,
     _["FishingLoss"]=FishingLoss,
     _["DetritalLoss"]=DetritalLoss,
     _["FishingThru"]=FishingThru,
     //_["PredSuite"]=PredSuite,
     //_["HandleSuite"]=HandleSuite,
     _["Qlink"]=Q1,
     _["GearCatch"]=GearCatch
     );

// Return is an Rcpp List     
   return(deriv);
}

// SplitSetPred function called in sim stanza initialize and update
// This function simply sums up across stanzas to get population-level
// Biomass, Numbers, and Consumption
// [[Rcpp::export]]
int SplitSetPred(List stanzas, List state){
  int isp, ist, ia, ieco;
  double Bt, pt, Nt;
  
  //stanza parameters
  const int Nsplit             = as<int>(stanzas["Nsplit"]);     
  const NumericVector Nstanzas = as<NumericVector>(stanzas["Nstanzas"]);
  NumericMatrix NageS          = as<NumericMatrix>(stanzas["NageS"]);
  NumericMatrix WageS          = as<NumericMatrix>(stanzas["WageS"]);
  NumericMatrix WWa            = as<NumericMatrix>(stanzas["WWa"]);
  NumericMatrix Age1           = as<NumericMatrix>(stanzas["Age1"]);
  NumericMatrix Age2           = as<NumericMatrix>(stanzas["Age2"]);
  NumericMatrix EcopathCode    = as<NumericMatrix>(stanzas["EcopathCode"]);
  NumericVector stanzaPred     = as<NumericVector>(stanzas["stanzaPred"]);

  //state parameters
  NumericVector state_BB = as<NumericVector>(state["BB"]);
  NumericVector state_NN = as<NumericVector>(state["NN"]);

  for (isp = 1; isp <= Nsplit; isp++){
    for (ist = 1; ist <= Nstanzas[isp]; ist++){
      ieco = EcopathCode(isp, ist);
      Bt = 1e-30;
      pt = 1e-30;
      Nt = 1e-30;
      for (ia = Age1(isp, ist); ia <= Age2(isp, ist); ia++){
        Bt = Bt + NageS(ia, isp) * WageS(ia, isp);
        pt = pt + NageS(ia, isp) * WWa(ia, isp);
        Nt = Nt + NageS(ia, isp);
      }
      state_BB[ieco] = Bt;
      state_NN[ieco] = Nt;
      stanzaPred[ieco] = pt;
    }
  }
  return(0);
}

// SplitUpdate function
// Update numbers, weight, and biomass for multistanza groups
// [[Rcpp::export]]
int SplitUpdate(List stanzas, List state, List forcing, List deriv, int yr, int mon){
  int isp, ist, ia, ieco=0, last, first;  //KYA 6/12/17 ieco=0 to stop annoying warning
  double Su, Gf, Nt;

  //stanza parameters
  const int Nsplit                   = as<int>(stanzas["Nsplit"]);     
  const NumericVector Nstanzas       = as<NumericVector>(stanzas["Nstanzas"]);
  const NumericVector vBM            = as<NumericVector>(stanzas["vBM"]);
  const NumericVector Wmat           = as<NumericVector>(stanzas["Wmat"]);
  //const NumericVector Wmat001        = as<NumericVector>(stanzas["Wmat001"]);  
  //const NumericVector Wmat50         = as<NumericVector>(stanzas["Wmat50"]);  
  //const NumericVector WmatSpread     = as<NumericVector>(stanzas["WmatSpread"]);  
  //const NumericVector Amat001        = as<NumericVector>(stanzas["Amat001"]);  
  //const NumericVector Amat50         = as<NumericVector>(stanzas["Amat50"]);  
  //const NumericVector AmatSpread     = as<NumericVector>(stanzas["AmatSpread"]);  
  const NumericVector baseEggsStanza = as<NumericVector>(stanzas["baseEggsStanza"]);
  const NumericVector RscaleSplit    = as<NumericVector>(stanzas["RscaleSplit"]);
  const NumericVector RzeroS         = as<NumericVector>(stanzas["RzeroS"]);
  const NumericVector RecPower       = as<NumericVector>(stanzas["RecPower"]);
  const NumericVector vBGFd          = as<NumericVector>(stanzas["vBGFd"]);
  const NumericVector SpawnEnergy    = as<NumericVector>(stanzas["SpawnEnergy"]);
  const NumericVector SpawnX         = as<NumericVector>(stanzas["SpawnX"]);
  const NumericVector baseSpawnBio   = as<NumericVector>(stanzas["baseSpawnBio"]);
  NumericVector SpawnBio             = as<NumericVector>(stanzas["SpawnBio"]);  
  NumericVector EggsStanza           = as<NumericVector>(stanzas["EggsStanza"]);  
  NumericMatrix NageS                = as<NumericMatrix>(stanzas["NageS"]);
  NumericMatrix WageS                = as<NumericMatrix>(stanzas["WageS"]);
  NumericMatrix SplitAlpha           = as<NumericMatrix>(stanzas["SplitAlpha"]);
  NumericMatrix WWa                  = as<NumericMatrix>(stanzas["WWa"]);
  NumericMatrix Age1                 = as<NumericMatrix>(stanzas["Age1"]);
  NumericMatrix Age2                 = as<NumericMatrix>(stanzas["Age2"]);
  NumericMatrix EcopathCode          = as<NumericMatrix>(stanzas["EcopathCode"]);
  NumericVector stanzaPred           = as<NumericVector>(stanzas["stanzaPred"]);

  //state parameters
  const NumericVector state_BB = as<NumericVector>(state["BB"]);
  
  //forcing parameters
  NumericMatrix force_byrecs   = as<NumericMatrix>(forcing["byrecs"]);

  //derivatives
  const NumericVector LossPropToB = as<NumericVector>(deriv["LossPropToB"]);
  const NumericVector FoodGain    = as<NumericVector>(deriv["FoodGain"]);

  for (isp = 1; isp <= Nsplit; isp++){
    // Update numbers and body weights
    SpawnBio[isp] = 0;
    for(ist = 1; ist <= Nstanzas[isp]; ist++){
      ieco = EcopathCode(isp, ist);
      Su = exp(-LossPropToB[ieco] / STEPS_PER_YEAR / state_BB[ieco]);
      Gf = FoodGain[ieco] / stanzaPred[ieco];
      for(ia = Age1(isp, ist); ia <= Age2(isp, ist); ia++){
        NageS(ia, isp) = NageS(ia, isp) * Su;
        WageS(ia, isp) = vBM[isp] * WageS(ia, isp) + Gf * SplitAlpha(ia, isp);
      // KYA 5/7/18 - started to add alternate maturity method, commented out for now
        //if(Wmat[isp]<0.0){
        //  if ((WageS(ia, isp)>Wmat001[isp])&&(ia>Amat001[isp])){
        //    SpawnBio[isp] += NageS(ia, isp) * WageS(ia, isp)/(1. + exp( 
        //                   -((WageS(ia, isp) - Wmat50[isp]) / WmatSpread[isp])   
        //                   -((ia             - Amat50[isp]) / AmatSpread[isp]) ));                   
        //  }
        //}
        //else{
          if(WageS(ia, isp) > Wmat[isp]){
            SpawnBio[isp] += NageS(ia, isp) * (WageS(ia, isp) - Wmat[isp]);
          };
        //}
      }
    }
    EggsStanza[isp] = SpawnBio[isp] * SpawnEnergy[isp] * SpawnX[isp] /
                      (SpawnX[isp] - 1.0 + (SpawnBio[isp] / baseSpawnBio[isp]));
    EggsStanza[isp] *= force_byrecs(yr * STEPS_PER_YEAR + mon, ieco);

    //Rprintf("%g %g %g %g %g",isp,SpawnBio[isp],EggsStanza[isp],)
    // Need to add monthly recruitment

    // now update n and wt looping backward over age
    last  = Age2(isp, Nstanzas[isp]);
    first = Age1(isp, 1);

    Nt = NageS(last, isp) + NageS(last - 1, isp);    
    if(Nt == 0){Nt = 1e-30;}
    
    WageS(last, isp) = (WageS(last, isp) * NageS(last, isp) + WageS(last - 1, isp) * 
                        NageS(last - 1, isp)) / Nt;
    NageS(last, isp) = Nt;

    for(ia = last - 1; ia > first; ia--){
      NageS(ia, isp) = NageS(ia - 1, isp);
      WageS(ia, isp) = WageS(ia - 1, isp);
    }

    //Apply number of eggs to youngest slot. Includes Walter recruit power
    if(baseEggsStanza[isp] > 0){
      NageS(first, isp) = RscaleSplit[isp] * RzeroS[isp] * pow(double(EggsStanza[isp] / 
                          baseEggsStanza[isp]), double(RecPower[isp]));
    }
    WageS(first, isp) = 0;

    //Uses generalized vonB (exponent is d)
    //Added for stability 4/13/07 (Unlucky Friday)
    for(ia = 0; ia <= last; ia++){
      WWa(ia, isp) = pow(double(WageS(ia, isp)), double(vBGFd[isp]));
    }
  }

return(0);
}