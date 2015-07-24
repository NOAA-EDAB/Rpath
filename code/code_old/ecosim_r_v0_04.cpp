//------------------------------------------------------------------------------
//Version 0.04 on 07-March-2013
//
// Ecosim_r is an implemetation of Ecosim written in C by Kerim Aydin
// (Kerim.Aydin@noaa.gov) and Sarah Gaichas (Sarah.Gaichas@noaa.gov), designed 
// for calling from R functions, based on published Ecosim algorithms.  
//
// This program consists of 5 functions:
//
// ecosim_run():     Main function called from R; it sets pointers of c variables
//                   to data structures created in R, calls the Adams-Basforth 
//                   routine, and returns run results to R.
// Adams-Basforth(): Strandard 1st order Adams-Basforth numerical integration using 
//                   current and previous step derivatives (e.g. Atkinson (1989) An 
//                   Introduction to Numerical Analysis).  
//                   Timestep is fixed to monthly to accord with delay-difference
//                   dynamics of aged groups; can be changed with some tweaks.  
//                   Organisms with production rates that are unstable at monthly
//                   timesteps are integrated with "fast equilibrium" integration
//                   (Walters and Christensen, pers. comm.)
// deriv_master():   Calculates derivative of ecosystem state based on curent
//                   biomass.  Includes predator/prey functional responses,
//                   fishing, and forcing.
// update_stanzas(): Runs delay-difference equations for age-pooled groups (monthly
//                   age pools of juveniles and adults).  Weight growth of age
//                   pools is based on consumption rates calculated by deriv_master.
// SplitSetPred():   Simple function to sum age-structure into total biomass for
//                   feeding responses. 
//
// VERSION 0.04 NOTES
// - Fishing is limited to forcing by Effort anomaly from Ecopath (by gear) or
//   by absolute F rate (on individual species).  Fishing
//   by absolute catch or effort is not currently implemented (the code for
//   absolute catch and by-target effort calculations is commented out in
//   deriv_master); no imput vectors exist in the calling R code for these yet.
// - Split pools are limited to 2 stages (Juvenile and Adult).
// - Functional response is that modified by Aydin (2004) and Gaichas et. al
//   (2012), and differs from EwE functional response in handling time and
//   Mediation implementation.  Foraging Time is implemented; handling time can 
//   be done by individual prey or by summed prey; predators can be by 
//   individual or "scramble". 
// - Fitting functions are not here.
// - Sensitivity (sense) routines and random generation of ecosystems should
//   be handled in R.
//
// Compilation:  This code has been compiled using gcc and interfacing with R
// in a Linux and Windows environment.  dll included (ecosim.dll) has been 
// tested on multiple windows machines and shouldn't need additional resources.
// Compile command, including linkages to libraries in R, contained in makefile.
//
// IMPORTANT: In Windows, this should be compiled by MinGW gcc compiler 
// (recommended by R).  In particular, Cygwin compilations do not seem to work.
// Compile by going to containing directory in MinGW and typing 'make' (make
// utility is part of MinGW).  In some cases, directories in makefile pointing
// to the compiler may need to be changed.
//
// In Linux, the created ecosim.dll will still work as a linux library, if you
// are particular you can change the dll extension to a linux library extension
// in the makefile and in the R calling code.  
// -----------------------------------------------------------------------------

// Load main headers (standard C packages) AND global declarations
#include "ecosim_r_v0_04.h"
#include <iostream>

extern "C" void ecosim_run(
   // List of pointers to vectors created in R, must be in the same order
   // as in the calling R script. 
      int *YEARS,
      int *BEGIN_YEAR,
      int *END_YEAR,
      int *numlist,
      int *flags,
      double *ratelist,
      int *ppind,
      double *pplist,
      int *fishind,
      double *fishlist,
      int *detind,
      double *detlist,
      int *juvind,
      double *juvlist,
      double *targlist,
      int *outflags,
      double *statelist,
      double *NageS,
      double *WageS,
      double *WWa,
      double *SplitAlpha,
      double *derivlist,
      double *force_byprey,
      double *force_bymort,
      double *force_byrecs,
      double *force_bysearch,
      double *FORCED_FRATE,
      double *FORCED_CATCH,
      double *out_BB,
      double *out_CC,
      double *out_SSB,
      double *out_rec
){
 struct newsim v;
 int r, t, byy, eyy;

 // Almost all of the Main routine is loading R vectors into the newsim structure v

 // Basic values, size of vectors, and single-valued flags
    v.init_run         = 0;  
    v.YEARS            = *YEARS;
    byy                = *BEGIN_YEAR;
    eyy                = *END_YEAR;
    v.NUM_GROUPS       = numlist[0];
    v.NUM_LIVING       = numlist[1];
    v.NUM_DEAD         = numlist[2];
    v.NUM_GEARS        = numlist[3];
    v.juv_N            = numlist[4];
    v.NumPredPreyLinks = numlist[5];
    v.NumFishingLinks  = numlist[6];
    v.NumDetLinks      = numlist[7];
    v.BURN_YEARS       = flags[0];
    v.COUPLED          = flags[1];
    v.CRASH_YEAR       = -1;

 // Indexing is by pointer arithmatic, for example ratelist + 2 * r is the same as ratelist[2*r]   
 
 //  Basic group values
     r = v.NUM_GROUPS + 1;
     v.B_BaseRef        = ratelist + 0 * r;
     v.MzeroMort        = ratelist + 1 * r;
     v.UnassimRespFrac  = ratelist + 2 * r;
     v.ActiveRespFrac   = ratelist + 3 * r;
     v.FtimeAdj         = ratelist + 4 * r;
     v.FtimeQBOpt       = ratelist + 5 * r;
     v.PBopt            = ratelist + 6 * r;
     v.NoIntegrate      = ratelist + 7 * r;
     v.HandleSelf       = ratelist + 8 * r;
     v.ScrambleSelf     = ratelist + 9 * r;
     v.PredTotWeight    = ratelist + 10 * r;
     v.PreyTotWeight    = ratelist + 11 * r;
     v.fish_Effort      = ratelist + 12 * r;

 // Predator-prey link parameters  
    r = v.NumPredPreyLinks + 1;
    v.PreyFrom         = ppind + 0 * r;
    v.PreyTo           = ppind + 1 * r;
    v.QQ               = pplist + 0 * r;
    v.DD               = pplist + 1 * r;
    v.VV               = pplist + 2 * r;
    v.HandleSwitch     = pplist + 3 * r;
    v.PredPredWeight   = pplist + 4 * r;
    v.PreyPreyWeight   = pplist + 5 * r;

 // fishing link parameters  
    r = v.NumFishingLinks + 1;
    v.FishingFrom         = fishind + 0 * r;
    v.FishingThrough      = fishind + 1 * r;
    v.FishingTo           = fishind + 2 * r;
    v.FishingQ            = fishlist + 0 * r;

 // detritus flow link parameters   
    r = v.NumDetLinks + 1;
    v.DetFrom          = detind + 0 * r;
    v.DetTo            = detind + 1 * r;
    v.DetFrac          = detlist + 0 * r;

 // juvenile/adult parameters  
    r = v.juv_N + 1; 
    v.JuvNum           = juvind + 0 * r;
    v.AduNum           = juvind + 1 * r;
    v.firstMoJuv       = juvind + 2 * r;
    v.lastMoJuv        = juvind + 3 * r;
    v.firstMoAdu       = juvind + 4 * r;
    v.lastMoAdu        = juvind + 5 * r;
    v.SpawnBio         = juvlist + 0 * r;
    v.EggsStanza       = juvlist + 1 * r;
    v.stanzaGGJuv      = juvlist + 2 * r;
    v.stanzaGGAdu      = juvlist + 3 * r;
    v.SpawnEnergy      = juvlist + 4 * r;
    v.SpawnX           = juvlist + 5 * r;
    v.SpawnAllocR      = juvlist + 6 * r;
    v.SpawnAllocG      = juvlist + 7 * r;
    v.recruits         = juvlist + 8 * r;
    v.RzeroS           = juvlist + 9 * r;
    v.baseEggsStanza   = juvlist +10 * r;
    v.baseSpawnBio     = juvlist +11 * r;
    v.Rbase            = juvlist +12 * r;
    v.RecMonth         = juvlist +13 * r;
    v.VonBD            = juvlist +14 * r;
    v.aduEqAgeZ        = juvlist +15 * r;
    v.juvEqAgeZ        = juvlist +16 * r;
    v.RecPower         = juvlist +17 * r;
    v.Wmat001          = juvlist +18 * r;
    v.Wmat50           = juvlist +19 * r;
    v.Amat001          = juvlist +20 * r;
    v.Amat50           = juvlist +21 * r;
    v.WmatSpread       = juvlist +22 * r;
    v.AmatSpread       = juvlist +23 * r;
    v.vBM              = juvlist +24 * r;
    v.RscaleSplit      = juvlist +25 * r;

 // state variables  
    r = v.NUM_GROUPS + 1;
    v.state_BB         = statelist + 0 * r;
    v.state_Ftime      = statelist + 1 * r;
    v.state_NN         = statelist + 2 * r;
    v.stanzaPred       = statelist + 3 * r;
    v.stanzaBasePred   = statelist + 4 * r;
 
 // state matrices for age structure.  rmat_trans switches from R order to C order 
    v.NageS            = rmat_trans(v.juv_N+1, MAX_MONTHS_STANZA+1, NageS);
    v.WageS            = rmat_trans(v.juv_N+1, MAX_MONTHS_STANZA+1,WageS);
    v.WWa              = rmat_trans(v.juv_N+1, MAX_MONTHS_STANZA+1, WWa);
    v.SplitAlpha       = rmat_trans(v.juv_N+1, MAX_MONTHS_STANZA+1, SplitAlpha);

 // Instantaneous state of the derivative for analysis
    r = v.NUM_GROUPS + 1;
    v.deriv_TotGain        = derivlist + 0 * r;
    v.deriv_TotLoss        = derivlist + 1 * r;
    v.deriv_LossPropToB    = derivlist + 2 * r;
    v.deriv_LossPropToQ    = derivlist + 3 * r;
 
    v.deriv_DerivT         = derivlist + 4 * r;
    v.deriv_dyt            = derivlist + 5 * r;
    v.deriv_biomeq         = derivlist + 6 * r;
 
    v.deriv_FoodGain       = derivlist + 7 * r;
    v.deriv_DetritalGain   = derivlist + 8 * r;
    v.deriv_FoodLoss       = derivlist + 9 * r;
    v.deriv_UnAssimLoss    = derivlist +10 * r;
    v.deriv_ActiveRespLoss = derivlist +11 * r;
 
    v.deriv_MzeroLoss      = derivlist +12 * r;
    v.deriv_DetritalLoss   = derivlist +13 * r;
    v.deriv_TotDetOut      = derivlist +14 * r;
    
    v.deriv_FishingLoss    = derivlist +15 * r;
    v.deriv_FishingGain    = derivlist +16 * r;
    v.deriv_FishingThru    = derivlist +17 * r;
 
    v.deriv_preyYY         = derivlist +18 * r;
    v.deriv_predYY         = derivlist +19 * r;
    cout << "v.deriv_preyYY" << v.deriv_preyYY;
    v.deriv_PredSuite      = derivlist +20 * r;
    v.deriv_HandleSuite    = derivlist +21 * r;
    v.deriv_TerminalF      = derivlist +22 * r;

  // fixed "forcing" such as non time-dependent control rules
     v.TARGET_BIO          = targlist + 0 * r;
     v.TARGET_FRATE        = targlist + 1 * r;
     v.ALPHA               = targlist + 2 * r;
     
	// Monthly forcing matrices
		 v.force_byprey   = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, force_byprey  );
     v.force_bymort   = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, force_bymort  );
     v.force_byrecs   = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, force_byrecs  );
     v.force_bysearch = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, force_bysearch);

  // Annual forcing matrices
     v.FORCED_FRATE   = rmat_trans(v.NUM_GROUPS+1, v.YEARS+1,    FORCED_FRATE );
     v.FORCED_CATCH   = rmat_trans(v.NUM_GROUPS+1, v.YEARS+1,    FORCED_CATCH );

	// monthly matrices for storing all outputs
     v.out_BB  = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, out_BB);
     v.out_CC  = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, out_CC);
     v.out_SSB = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, out_SSB);
     v.out_rec = rmat_trans(v.NUM_GROUPS+1, v.YEARS*12+1, out_rec);
     
  // Making sure range of years is in bounds  
	   if (byy <         0) byy = 0;
	   if (byy > v.YEARS-1) byy = v.YEARS-1;
	   if (eyy <=      byy) eyy = byy + 1;
	   if (eyy > v.YEARS  ) eyy = v.YEARS;

  // Run Ecosim model contained in structure v from year byy to year eyy	 
     Adams_Basforth(&v, byy, eyy);
  
  // Output status flags
     outflags[0] = v.CRASH_YEAR;

	// Return to R at this point
}

// ---------------------------------------------------------------------------- 
int Adams_Basforth (struct newsim *v, int StartYear, int EndYear){
     
int y, m, c, j, dd;
int sp, t, i, ageMo, s, link,prey,pred,links,gr;
int inbound;
double old_B, new_B, nn, ww, bb, pd;
double newbio; float ystep;
double bioanom, bioratio;
unsigned int LL;

   // Update sums of split groups to total biomass for derivative calcs
      SplitSetPred(v); 
      
 	 // If not starting from previous timestep, call derivative twice to
 	 // set deriv and deriv(t-1)
      if (v->init_run){
         deriv_master(v,0,0,0);
         deriv_master(v,0,0,0);
      }

    // MAIN LOOP STARTS HERE     
       for (y=StartYear; y<EndYear; y++){                                  
         for (m=0; m<STEPS_PER_YEAR; m++){
                     
           // dd is index for saving monthly output
				      dd = y*STEPS_PER_YEAR+m;
           
           // save current state to output matrix
           // For non-split pools, SSB is output output as the same as B.
					 // For adult split pools, it is overwritten with "actual" SSB, 
					 // for juvs it is 0. 
              for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){ 
                  v->out_BB[sp][dd]  = v->state_BB[sp];
                  v->out_SSB[sp][dd] = v->state_BB[sp];
                  v->out_rec[sp][dd] = v->state_BB[sp];
              }
              for (i = 1; i<= v->juv_N; i++){
                  v->out_SSB[v->JuvNum[i]][dd] = 0.0;
									v->out_SSB[v->AduNum[i]][dd] = v->SpawnBio[i];
									v->out_rec[v->AduNum[i]][dd] = v->NageS[i][v->firstMoAdu[i]] * v->WageS[i][v->firstMoAdu[i]];
              }
       
			    // note: code for sub-monthly steps has been removed    
 		      // for (d=0; d<STEPS_PER_MONTH; d++){

 		          // Calculate Derivative for a given timestep
 	               deriv_master(v,y,m,0);

 	            // Loop through species, applying Adams-Basforth or fast
 	            // equilibrium depending on species.
 	            
 								 for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){                    

                  // Adjust feeding time 
 	                   if (v->NoIntegrate[sp] < 0){pd = v->stanzaPred[sp];}
 	                   else                         {pd = v->state_BB[sp];}	                      					       
                     if ((pd > 0) && (v->deriv_FoodGain[sp]>0)){
 										      v->state_Ftime[sp] = 0.1 + 0.9 * v->state_Ftime[sp] * 
                                 ((1.0 - v->FtimeAdj[sp]/(double)STEPS_PER_MONTH) 
 																+ v->FtimeAdj[sp]/(double)STEPS_PER_MONTH * 
 										            v->FtimeQBOpt[sp] / (v->deriv_FoodGain[sp] / pd));
                          }	
                     MAX_THRESHOLD(v->state_Ftime[sp], 2.0);
                      
                  // Biomass update for non-split groups
 								     old_B=v->state_BB[sp];
 								     new_B=v->state_BB[sp];
 								     if (v->NoIntegrate[sp] == 0){
 										 // 'Fast equilibrium'
												new_B        = (1.0 - SORWT) * v->deriv_biomeq[sp] 
 												               + SORWT * v->state_BB[sp];
 										 }
 										 else if (v->NoIntegrate[sp] == sp){
                         // Adams-Basforth timestep 
                            new_B = v->state_BB[sp] + (DELTA_T / 2.0) * 
                                   (3.0 * (v->deriv_DerivT[sp]) - v->deriv_dyt[sp]); 
 										       }
 										       
						      // If the new biomass goes to infinity or something, set a
									// flag to exit the loop and return the problem. 
						         if (isnan(new_B) || isinf(new_B)) {
                        v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;
                     }
                  
									// If the run is during the "burn-in" years, and biomass goes
									// into the discard range, set flag to exit the loop.  
                     if (y<v->BURN_YEARS){
                       if  ( (new_B < v->B_BaseRef[sp]*LO_DISCARD) ||
                             (new_B > v->B_BaseRef[sp]*HI_DISCARD) ){
                           v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;                         
                       }                          
                     }
                    
                 // Otherwise if biomass goes very large or very small, set to 
                 // min or max.
								    MIN_THRESHOLD(new_B,v->B_BaseRef[sp]*EPSILON);
                    MAX_THRESHOLD(new_B,v->B_BaseRef[sp]*BIGNUM);                    
                    
 										v->state_BB[sp] = new_B;                                
 						        
                 // sum up Catch at every time step (catch is cumulative)
 								    if (fabs(v->state_BB[sp]/old_B-1.0) > EPSILON){
 								       v->out_CC[sp][dd] = ((v->deriv_FishingLoss[sp] * DELTA_T)/old_B) *
 									                       (v->state_BB[sp] - old_B) /
 										 	    							 log(v->state_BB[sp]/old_B);
 									  }
 									  else {v->out_CC[sp][dd] = (v->deriv_FishingLoss[sp] * DELTA_T);
 								    }		  
  									         
 								 } // End of species loop
		    		       
 					//}  // End of days loop
           
          // Monthly Stanza (split pool) update
 					   update_stanzas(v,y,m+1);
             SplitSetPred(v);
           
				 // As above for non-split groups, check for out-of-bounds, numerical
				 // errors, or burn-in model stoppages for juvs and then adults.   
            for (i = 1; i<= v->juv_N; i++){
                sp = v->JuvNum[i];
                new_B = v->state_BB[sp];
						    if (isnan(new_B) || isinf(new_B)) {
                   v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;
                }                  
                if (y<v->BURN_YEARS){
                    if  ( (new_B < v->B_BaseRef[sp]*LO_DISCARD) ||
                          (new_B > v->B_BaseRef[sp]*HI_DISCARD) ){
                         v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;                           
                    }                          
                }                
                sp    = v->AduNum[i];
                new_B = v->state_BB[sp];
						    if (isnan(new_B) || isinf(new_B)) {
                    v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;
                }     
                if (y<v->BURN_YEARS){
                   if  ( (new_B < v->B_BaseRef[sp]*LO_DISCARD) ||
                         (new_B > v->B_BaseRef[sp]*HI_DISCARD) ){
                       v->CRASH_YEAR=y; y=EndYear; m=STEPS_PER_YEAR;                           
                   }                          
                }                
            }
                           
 		    }  // End of main months loop
     
   }// End of years loop
       
   // If you get to the end and haven't crashed, save biomass in the
   // final time slot before exiting.
      if (v->CRASH_YEAR<0){
         dd = EndYear*STEPS_PER_YEAR;
         for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){ 
              v->out_BB[sp][dd]  = v->state_BB[sp];
              v->out_SSB[sp][dd] = v->state_BB[sp];
              v->out_rec[sp][dd] = v->state_BB[sp];
         }
         for (i = 1; i<= v->juv_N; i++){
              v->out_SSB[v->JuvNum[i]][dd] = 0.0;
						  v->out_SSB[v->AduNum[i]][dd] = v->SpawnBio[i];
						  v->out_rec[v->AduNum[i]][dd] = v->NageS[i][v->firstMoAdu[i]] * v->WageS[i][v->firstMoAdu[i]];
         }
      }  
      
return(0);
} 
//------------------------------------------------------------------------------
// SplitSetPred function called in sim stanza initialize and update
// This function simply sums up across juvenile and adult age structure to get 
// population-level Biomass, Numbers, and Consumption 
int SplitSetPred(struct newsim *v){

int ageMo, i;  
double Bt, pt, Nt;

	// loop over split pools
     for (i = 1; i<= v->juv_N; i++){
         Bt = 1e-30;
         pt = 1e-30;
         Nt = 1e-30;
      // loop over juv monthly age classes
         for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){
             Bt = Bt + v->NageS[i][ageMo] * v->WageS[i][ageMo];
             pt = pt + v->NageS[i][ageMo] * v->WWa[i][ageMo];
             Nt = Nt + v->NageS[i][ageMo];
             
         }
         v->state_BB[v->JuvNum[i]]   = Bt;
         v->stanzaPred[v->JuvNum[i]] = pt;
         v->state_NN[v->JuvNum[i]]   = Nt;
    
      // loop over adult monthly age classes
         Bt = 1e-30;
         pt = 1e-30;
         Nt = 1e-30;       
         for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++){
             Bt = Bt + v->NageS[i][ageMo] * v->WageS[i][ageMo];
             pt = pt + v->NageS[i][ageMo] * v->WWa[i][ageMo];
             Nt = Nt + v->NageS[i][ageMo];
         }
         v->state_BB[v->AduNum[i]] = Bt;
         v->stanzaPred[v->AduNum[i]] = pt;
         v->state_NN[v->AduNum[i]] = Nt;    
         v->stanzaGGJuv[i] =  v->deriv_FoodGain[v->JuvNum[i]]/v->stanzaPred[v->JuvNum[i]]; 
         v->stanzaGGAdu[i] =  v->deriv_FoodGain[v->AduNum[i]]/v->stanzaPred[v->AduNum[i]]; 
     }  

return(0);
}

// ----------------------------------------------------------------------------
// Update juvenile adult or "stanza" age structure during sim run 
// on monthly timesteps (not hardwiring months, but recommended)
 
 int update_stanzas(struct newsim *v, int yr, int mon){
 
int ageMo, i;  
double Su, Gf, Nt; 
double propRepro;
     
  // loop over split species groups to update n, wt, biomass in sim  
     for (i = 1; i<= v->juv_N; i++){
         v->SpawnBio[i] = 0;
      // loop over juv timesteps with juv survival and growth, also calculate Spawning Biomass
         Su = exp(-v->deriv_LossPropToB[v->JuvNum[i]]/STEPS_PER_YEAR/v->state_BB[v->JuvNum[i]]);    
         Gf = v->deriv_FoodGain[v->JuvNum[i]]/v->stanzaPred[v->JuvNum[i]];  // BB should be stanzaPred            
 				 v->stanzaGGJuv[i] = Gf; 						         
         for (ageMo=v->firstMoJuv[i]; ageMo<=v->lastMoJuv[i]; ageMo++){   
             v->NageS[i][ageMo] = v->NageS[i][ageMo] * Su;
             v->WageS[i][ageMo] = v->vBM[i] * v->WageS[i][ageMo] + Gf * v->SplitAlpha[i][ageMo];
             if ( (v->WageS[i][ageMo] > v->Wmat001[i]) && (ageMo > v->Amat001[i])) {
                 v->SpawnBio[i]  += v->WageS[i][ageMo] * v->NageS[i][ageMo] /
                                 (1. + exp(- ((v->WageS[i][ageMo]-v->Wmat50[i])/v->WmatSpread[i]) 
 																		      - ((ageMo          -v->Amat50[i])/v->AmatSpread[i])));                             
             }            
         }
 				      
      // loop over adult timesteps with adult survival and growth, also calculate Spawning Biomass
         Su = exp(-v->deriv_LossPropToB[v->AduNum[i]]/STEPS_PER_YEAR/v->state_BB[v->AduNum[i]]);
         Gf = v->deriv_FoodGain[v->AduNum[i]]/v->stanzaPred[v->AduNum[i]];   //BB should be stanzaPred         
         v->stanzaGGAdu[i] = Gf; 
         for (ageMo=v->firstMoAdu[i]; ageMo<=v->lastMoAdu[i]; ageMo++){
             v->NageS[i][ageMo] = v->NageS[i][ageMo] * Su;
             v->WageS[i][ageMo] = v->vBM[i] * v->WageS[i][ageMo] + Gf * v->SplitAlpha[i][ageMo];
 						if ( (v->WageS[i][ageMo] > v->Wmat001[i]) && (ageMo > v->Amat001[i])) {
                 v->SpawnBio[i]  += v->WageS[i][ageMo] * v->NageS[i][ageMo] /
                                 (1. + exp(- ((v->WageS[i][ageMo]-v->Wmat50[i])/v->WmatSpread[i]) 
 																		      - ((ageMo          -v->Amat50[i])/v->AmatSpread[i])));                             
             }     	       
         }
         
      // KYA This is a Beverton-Holt curve between Spawning Biomass and
      // Number of Eggs Produced; generally near-linear with default parameters        
         v->EggsStanza[i] = v->SpawnBio[i] * v->SpawnEnergy[i] * v->SpawnX[i] /
 				                             (v->SpawnX[i] - 1.0 + 
 			  															 (v->SpawnBio[i] / v->baseSpawnBio[i]));
 		 // Forcing for eggs applied here															 
 			   v->EggsStanza[i] *= v->force_byrecs[v->JuvNum[i]][yr*STEPS_PER_YEAR+mon];		

		 // If recruitment happens all at once (in a single month) that happens here)
     // Otherwise, recruitment is throughout the year.
 			  if (v->RecMonth[i]>0){
 				   if (mon == v->RecMonth[i]) {v->EggsStanza[i] *=1. ;}
 				 else                         {v->EggsStanza[i] *=0. ;}
 				}
     
     // now update n and wt looping backward over age
         Nt = v->NageS[i][v->lastMoAdu[i]] + v->NageS[i][v->lastMoAdu[i]-1];
         if (Nt == 0){Nt = 1e-30;}
         v->WageS[i][v->lastMoAdu[i]] = (v->WageS[i][v->lastMoAdu[i]]*v->NageS[i][v->lastMoAdu[i]] + 
                            v->WageS[i][v->lastMoAdu[i]-1]*v->NageS[i][v->lastMoAdu[i]-1])/ Nt;
         v->NageS[i][v->lastMoAdu[i]] = Nt;
         
         for (ageMo=v->lastMoAdu[i]-1; ageMo>v->firstMoJuv[i]; ageMo--) {
             v->NageS[i][ageMo] = v->NageS[i][ageMo-1];
             v->WageS[i][ageMo] = v->WageS[i][ageMo-1];   
         }
 
      // finally apply number of eggs to youngest juvenile slot. including Carl's recruit power        
         if (v->baseEggsStanza[i] > 0){
             v->NageS[i][v->firstMoJuv[i]] = v->RscaleSplit[i] * v->RzeroS[i] * 
                                       pow(double(v->EggsStanza[i] / v->baseEggsStanza[i]), 
                                       double(v->RecPower[i]));
         }
         v->WageS[i][v->firstMoJuv[i]] = 0;
 
     // Uses generalized vonB (exponent is d).    
     // ADDED FOR STABILITY 4/13/07 (Unlucky Friday)
        for (ageMo = 0; ageMo <= v->lastMoAdu[i]; ageMo++){
             v->WWa[i][ageMo] = pow(double(v->WageS[i][ageMo]), double(v->VonBD[i])); 
         }           
     }  // end of Juvenile Loop  

return(0); 
}

// ----------------------------------------------------------------------------
// Deriv master calculates the biomass dynamics derivative

int deriv_master(struct newsim *v, int y, int m, int d){

// Functional response vars     
   int sp, links, prey, pred, i;
   double Master_Density, Q;
   unsigned int LL;
// Fishing Vars
   int gr,dest;
   double caught, totQ;


     // Some derivative parts need to be set to zero
         LL=(v->NUM_GROUPS+1)*sizeof(double);
         memset(v->deriv_FoodLoss,		  0,LL);
         memset(v->deriv_FoodGain,		  0,LL);
         memset(v->deriv_UnAssimLoss,	  0,LL);
         memset(v->deriv_ActiveRespLoss,0,LL);   
         memset(v->deriv_DetritalGain,  0,LL);
         memset(v->deriv_FishingGain,   0,LL);
         memset(v->deriv_MzeroLoss,     0,LL);
         memset(v->deriv_FishingLoss,   0,LL);
         memset(v->deriv_DetritalLoss,  0,LL);
         memset(v->deriv_FishingThru,   0,LL);
 	       memset(v->deriv_PredSuite,     0,LL);
 	       memset(v->deriv_HandleSuite,   0,LL);
 				     			  
     //  Set YY = B/B(ref) for functional response; note that foraging time is
		 //  used here as a biomass modifier before the main functional response  
         for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
 		        v->deriv_preyYY[sp]=v->state_Ftime[sp] * v->state_BB[sp]/v->B_BaseRef[sp];
 		        v->deriv_predYY[sp]=v->state_Ftime[sp] * v->state_BB[sp]/v->B_BaseRef[sp];
 		     }    
 		     // The sun always has a biomass of 1 for Primary Production
					  v->deriv_preyYY[0]=1.0;  v->deriv_predYY[0]=1.0;
 		 		
 		// Set functional response biomass for juvenile and adult groups (including foraging time) 
 		   for (i=1; i<=v->juv_N; i++){
 		     if (v->stanzaBasePred[v->JuvNum[i]]>0){
 		        v->deriv_predYY[v->JuvNum[i]] = v->state_Ftime[v->JuvNum[i]] * 
 				       v->stanzaPred[v->JuvNum[i]]/v->stanzaBasePred[v->JuvNum[i]];
 		        v->deriv_predYY[v->AduNum[i]] = v->state_Ftime[v->AduNum[i]] * 
 				        v->stanzaPred[v->AduNum[i]]/v->stanzaBasePred[v->AduNum[i]];
 		     }
 		   }
    // add "mediation by search rate" KYA 7/8/08
        for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
            v->deriv_predYY[sp] *= v->force_bysearch[sp][y*STEPS_PER_YEAR+m]; 
        }
        
 	  // Summed predator and prey for joint handling time and/or scramble functional response
 	     for (links=1; links<=v->NumPredPreyLinks; links++){
 		          prey = v->PreyFrom[links];
 		          pred = v->PreyTo[links];
               v->deriv_PredSuite[prey]   += v->deriv_predYY[pred] * v->PredPredWeight[links];
               v->deriv_HandleSuite[pred] += v->deriv_preyYY[prey] * v->PreyPreyWeight[links];
 		  }

	 // Main loop to calculate functional response for each predator/prey link
      for (links=1; links<=v->NumPredPreyLinks; links++){
 		      prey = v->PreyFrom[links];
 		      pred = v->PreyTo[links];
 
      // MAIN FUNCTIONAL RESPONSE SET HERE TO OUTPUT TOTAL CONSUMPTION (Q); use 1 version only
         // (1) This is EwE Classic (c)(r)(tm)(4.0beta)(all rights reserved)
 				 //  Master_Density = v->deriv_predYY[pred];	 
     		 //   Q = v->QQ[links] * v->XX[links] * v->deriv_predYY[pred] * v->deriv_preyYY[prey] / 
         //       (Master_Density + ((v->XX[links] - 1.0) * 
 				 //	                   (1.0 + v->state_Ftime[prey]) / 2.0) );
  
         // (2) Multiplicative version which is too weird for me.
         // Q =   v->QQ[links] * 
         //     ( v->XX[links] * v->deriv_predYY[pred] / ((v->XX[links] - 1.0) + v->deriv_predYY[pred])     )*
         //     ( v->HH[links] * v->deriv_preyYY[prey] / ((v->HH[links] - 1.0) + v->deriv_preyYY[prey])     )*
 				 // 		( v->DD[links]                 / ((v->DD[links] - 1.0) + v->deriv_HandleSuite[pred]) )*
 				 // 		( v->VV[links]                 / ((v->VV[links] - 1.0) + v->deriv_PredSuite[prey])   );
         
 				 // (3) Additive version: primary used and published in Aydin (2004) 
            // KYA 3/2/2012 setting "COUPLED" to zero means species are density dependent
            // (based on their own vul) but don't interact otherwise.  This can magically
            // create and destroy energy in the system but makes them act like a set
            // of independent surplus production models for comparison purposes 
         Q =   v->QQ[links] * v->deriv_predYY[pred] * pow(v->deriv_preyYY[prey], v->COUPLED * v->HandleSwitch[links]) *
 				      ( v->DD[links] / ( v->DD[links] - 1.0 + 
 						                     pow(v->HandleSelf[pred] * v->deriv_preyYY[prey]   + 
 																 (1. - v->HandleSelf[pred]) * v->deriv_HandleSuite[pred],
                                                      v->COUPLED * v->HandleSwitch[links])) )*
             ( v->VV[links] / ( v->VV[links] - 1.0 + 
                                  v->ScrambleSelf[pred] * v->deriv_predYY[pred] + 
 						                     (1. - v->ScrambleSelf[pred]) * v->deriv_PredSuite[prey]) );
        
			 // Include any Forcing by prey   
 				  Q *= v->force_byprey[prey][y*STEPS_PER_YEAR+m]; 

			 // If model is uncoupled, food loss doesn't change with prey or predator levels.
				if (v->COUPLED){  v->deriv_FoodLoss[prey]  += Q; }
				           else{  v->deriv_FoodLoss[prey]  += v->state_BB[prey] * v->QQ[links]/v->B_BaseRef[prey]; }
        
			 // Energy Accounting
				  v->deriv_FoodGain[pred]           += Q;
          v->deriv_UnAssimLoss[pred]        += Q * v->UnassimRespFrac[pred]; 
          v->deriv_ActiveRespLoss[pred]     += Q * v->ActiveRespFrac[pred];  												 
 		  }
 		
     // Mzero Mortality
     for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
         v->deriv_MzeroLoss[sp] = v->MzeroMort[sp] * v->state_BB[sp];
     }	   
 		
		 
		// MOST OF THE FOLLOWING FISHING SPECIFICATION METHODS ARE NOT SUPPORTED
		// BY THE R-CODE, only fishing by effort (for gear) or by F-rate (for
		// species) is supported at the end.
		//
		// BY CURRENT R-CODE.  ONLY    
 		// RFISH for (gr=v->NUM_LIVING+v->NUM_DEAD+1; gr<=v->NUM_GROUPS; gr++){
    // RFISH This sets EFFORT by time series of gear-target combinations
 		// RFISH 		   for (gr=v->NUM_LIVING+v->NUM_DEAD+1; gr<=v->NUM_GROUPS; gr++){
    // RFISH if -1 is an input value, uses TERMINAL F (last non-negative F) 		
     //RFISH 		   for (gr=v->NUM_LIVING+v->NUM_DEAD+1; gr<=v->NUM_GROUPS; gr++){
     //RFISH 		       if (y+m+d == 0){v->fish_Effort[gr]=1.0;}
     //RFISH 		       else           {v->fish_Effort[gr]=1.0;} // NOTE DEFAULT!  THIS CAN BE CHANGED TO 1.0
     //RFISH        // Added 7/8/08 for forced effort
     //RFISH            if (v->FORCED_EFFORT[gr][y] > -0.001) 
     //RFISH 					    {v->fish_Effort[gr]=v->FORCED_EFFORT[gr][y];}
     //RFISH 
     //RFISH 			     if ((v->FORCED_TARGET[gr]>0) && (v->FORCED_CATCH[gr][y]>-EPSILON)){
     //RFISH 			        totQ = 0.0;
     //RFISH 			        sp   = v->FORCED_TARGET[gr];
     //RFISH 			        for (links=1; links<=v->NumFishingLinks; links++){
     //RFISH     					    if ((v->FishingThrough[links] == gr) && 
     //RFISH 		    					    (v->FishingFrom[links]) == sp){
     //RFISH 				    					totQ += v->FishingQ[links];
     //RFISH 									}
     //RFISH 						  }
     //RFISH 						  v->fish_Effort[gr] = v->FORCED_CATCH[gr][y]/ 
     //RFISH 							                  (totQ * v->state_BB[sp]);
     //RFISH 							if (v->FORCED_CATCH[gr][y] >= v->state_BB[sp])
     //RFISH 							   {v->fish_Effort[gr] = (1.0-EPSILON)*(v->state_BB[sp])/ 
     //RFISH 				    			                  (totQ * v->state_BB[sp]);}
     //RFISH 					 }
     //RFISH 					 // By putting F after catch, Frates override absolute catch
     //RFISH 			     if ((v->FORCED_FTARGET[gr]>0) && (v->FORCED_FRATE[gr][y]>-EPSILON)){
     //RFISH 			        totQ = 0.0;
     //RFISH 			        sp   = v->FORCED_FTARGET[gr];
     //RFISH 			        for (links=1; links<=v->NumFishingLinks; links++){
     //RFISH     					    if ((v->FishingThrough[links] == gr) && 
     //RFISH 		    					    (v->FishingFrom[links]) == sp){
     //RFISH 				    					totQ += v->FishingQ[links];
     //RFISH 									}
     //RFISH 						  }
     //RFISH 						  v->fish_Effort[gr] = v->FORCED_FRATE[gr][y]/totQ;
     //RFISH 							//if (FORCED_CATCH[gr][y] >= v->state_BB[sp])
     //RFISH 							//   {v->fish_Effort[gr] = (1.0-EPSILON)*(v->state_BB[sp])/ 
     //RFISH 				    //			                  (totQ * v->state_BB[sp]);}
     //RFISH 					 }					 
     //RFISH 					 
     //RFISH 				   //if ((y==0) && (m==0) && (d==0)){
     //RFISH 					 //    cout << path_species[gr] << " " << FORCED_TARGET[gr] << " " << path_species[sp] << " " 
     //RFISH 					 //	      << v->state_BB[sp] << " " << FORCED_CATCH[gr][y] << " "
     //RFISH 					 //		      << v->fish_Effort[gr] << endl;
     //RFISH 					 //} 					 
     //RFISH 			 }
 			 					 			 					 
   // Apply specified Effort by Gear to catch (using Ecopath-set Q)
      for (links=1; links<=v->NumFishingLinks; links++){
 				 prey = v->FishingFrom[links];
 				 gr   = v->FishingThrough[links];
 				 dest = v->FishingTo[links];
 				 caught = v->FishingQ[links] * v->fish_Effort[gr] * v->state_BB[prey]; 
          v->deriv_FishingLoss[prey] += caught;
          v->deriv_FishingThru[gr]   += caught;
          v->deriv_FishingGain[dest] += caught;
 		}		
 
    //  Special "CLEAN" fisheries assuming q=1, so specified input is Frate
        for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
             caught = v->FORCED_CATCH[sp][y] + v->FORCED_FRATE[sp][y] * v->state_BB[sp];
             // KYA Aug 2011 removed terminal effort option to allow negative fishing pressure 
                // if (caught <= -EPSILON) {caught = v->deriv_TerminalF[sp] * v->state_BB[sp];}
             if (caught>=v->state_BB[sp]){caught=(1.0-EPSILON)*(v->state_BB[sp]);}
             v->deriv_FishingLoss[sp] += caught;
             v->deriv_FishingThru[0]  += caught;
             v->deriv_FishingGain[0]  += caught;
             v->deriv_TerminalF[sp] = caught/v->state_BB[sp];
        }
    
    // KINKED CONTROL RULE - NEEDS INPUT of TARGET BIOMASS and TARGET CATCH
        double RefBio, maxcaught;
        for (sp=1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
            if (v->TARGET_BIO[sp] > EPSILON){
               RefBio    = v->state_BB[sp]/v->TARGET_BIO[sp];
               maxcaught = v->TARGET_FRATE[sp] * v->state_BB[sp];         
               if      (RefBio > 1.0)           {caught = maxcaught;}
               else if (RefBio >= v->ALPHA[sp]) {caught = maxcaught * (RefBio - v->ALPHA[sp])/(1.0 - v->ALPHA[sp]);}
               else                             {caught = 0.0;}
               v->deriv_FishingLoss[sp] += caught;
               v->deriv_FishingThru[0]  += caught;
               v->deriv_FishingGain[0]  += caught;
               v->deriv_TerminalF[sp] = caught/v->state_BB[sp];
            }          
        }
        
   // DETRITUS  - note: check interdetrital flow carefully, have had some issues
 	 // (check by ensuring equlibrium run stays in equilibrium)
	 	  int liv, det;
 		  double flow;
      for (links=1; links<=v->NumDetLinks; links++){
          liv  = v->DetFrom[links];
          det  = v->DetTo[links];
          flow = v->DetFrac[links] * (v->deriv_MzeroLoss[liv] + v->deriv_UnAssimLoss[liv]);
          v->deriv_DetritalGain[det] += flow;
          if (liv > v->NUM_LIVING) {v->deriv_DetritalLoss[liv] += flow; }
      }
      for (sp=v->NUM_LIVING+1; sp<=v->NUM_LIVING+v->NUM_DEAD; sp++){
         v->deriv_MzeroLoss[sp] = 0.0;
      }  
    
  // Add mortality forcing
     for (i=1; i<=v->NUM_DEAD+v->NUM_LIVING; i++){
        v->deriv_FoodLoss[i]  *= v->force_bymort[i][y*STEPS_PER_YEAR+m];
        v->deriv_MzeroLoss[i] *= v->force_bymort[i][y*STEPS_PER_YEAR+m];
     }
  
	// Sum up derivitive parts; move previous derivative to dyt        
     for (i=1; i<=v->NUM_DEAD+v->NUM_LIVING; i++){
         v->deriv_dyt[i]=v->deriv_DerivT[i];
         v->deriv_TotGain[i] = v->deriv_FoodGain[i] + v->deriv_DetritalGain[i] + 
				                                              v->deriv_FishingGain[i];      
         v->deriv_LossPropToQ[i] = v->deriv_UnAssimLoss[i]  + v->deriv_ActiveRespLoss[i];
         v->deriv_LossPropToB[i] = v->deriv_FoodLoss[i]     + v->deriv_MzeroLoss[i] +
                                v->deriv_FishingLoss[i]  + v->deriv_DetritalLoss[i]; 
                  
         v->deriv_TotLoss[i] = v->deriv_LossPropToQ[i] + v->deriv_LossPropToB[i];    
         v->deriv_DerivT[i]  = v->deriv_TotGain[i] - v->deriv_TotLoss[i]; 
         
 	   // Set biomeq for "fast equilibrium" of fast variables
  	    if (v->state_BB[i] > 0) {
            v->deriv_biomeq[i] = v->deriv_TotGain[i] / 
 					                 (v->deriv_TotLoss[i] / v->state_BB[i]);
        }          
     }
     
return 0;
}
// ----------------------------------------------------------------------------