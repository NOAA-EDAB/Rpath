#include <Rcpp.h>
using namespace Rcpp;
 
// Various Constants
#define STEPS_PER_MONTH    1        // Integration Steps per month (speed (5) vs. accuracy (10 to 30))
#define MAX_MONTHS_STANZA  400      // Max number age pools within stanza calcs (months)
#define LO_DISCARD         1e-4     // Discard Limit for sense
#define HI_DISCARD         1e+4     // Discard Limit for sense

// Some fixed definitions, shouldn't need to tweak
#define STEPS_PER_YEAR     12                   // Steps per year between stanzas ("months")
#define MIN_THRESHOLD(x,y) if ((x)<(y))(x)=(y)  // Set x to min y if too small
#define MAX_THRESHOLD(x,y) if ((x)>(y))(x)=(y)  // Set x to max y if too big
#define DELTA_T            0.083333333333333333
#define SORWT              0.5                  // Fast Equilibrium smoother
#define EPSILON            1E-8                 // Test threshold for "too close to zero"
#define BIGNUM             1E+8                 // Test threshold for "too big"

// [[Rcpp::export]] 
int deriv_test(List mod){
   return(3);
}

// Deriv master calculates the biomass dynamics derivative
// [[Rcpp::export]]
int deriv_master(List mod, int y, int m, int d){
 //if (!mod.inherits("Rpath.sim")) stop("Input must be a Rpath model");

  // Functional response vars     
  int sp, links, prey, pred, i;
  double caught, Q;
  //unsigned int LL;
  // Fishing vars
  int gr, dest;
  // double Master_Density, totQ;
  
  // Parse out List mod
  int NUM_GROUPS                 = as<int>(mod["NUM_GROUPS"]);
  int NUM_LIVING                 = as<int>(mod["NUM_LIVING"]);
  int NUM_DEAD                   = as<int>(mod["NUM_DEAD"]);
  int NumPredPreyLinks           = as<int>(mod["NumPredPreyLinks"]);
  int NumFishingLinks            = as<int>(mod["NumFishingLinks"]);
  int NumDetLinks                = as<int>(mod["NumDetLinks"]);
  int juv_N                      = as<int>(mod["juv_N"]);
  int COUPLED                    = as<int>(mod["COUPLED"]);
  
  NumericVector B_BaseRef        = as<NumericVector>(mod["B_BaseRef"]);
  NumericVector MzeroMort        = as<NumericVector>(mod["MzeroMort"]);
  NumericVector UnassimRespFrac  = as<NumericVector>(mod["UnassimRespFrac"]);
  NumericVector ActiveRespFrac   = as<NumericVector>(mod["ActiveRespFrac"]);
  NumericVector fish_Effort      = as<NumericVector>(mod["fish_Effort"]);
  NumericVector HandleSelf       = as<NumericVector>(mod["HandleSelf"]);
  NumericVector ScrambleSelf     = as<NumericVector>(mod["ScrambleSelf"]);
  NumericVector PreyFrom         = as<NumericVector>(mod["PreyFrom"]);
  NumericVector PreyTo           = as<NumericVector>(mod["PreyTo"]);
  NumericVector QQ               = as<NumericVector>(mod["QQ"]);
  NumericVector DD               = as<NumericVector>(mod["DD"]);
  NumericVector VV               = as<NumericVector>(mod["VV"]);
  NumericVector HandleSwitch     = as<NumericVector>(mod["HandleSwitch"]);
  NumericVector PredPredWeight   = as<NumericVector>(mod["PredPredWeight"]);
  NumericVector PreyPreyWeight   = as<NumericVector>(mod["PreyPreyWeight"]);
  NumericVector FishFrom         = as<NumericVector>(mod["FishFrom"]);
  NumericVector FishThrough      = as<NumericVector>(mod["FishThrough"]);
  NumericVector FishQ            = as<NumericVector>(mod["FishQ"]);
  NumericVector FishTo           = as<NumericVector>(mod["FishTo"]);
  NumericVector DetFrac          = as<NumericVector>(mod["DetFrac"]);
  NumericVector DetFrom          = as<NumericVector>(mod["DetFrom"]);
  NumericVector DetTo            = as<NumericVector>(mod["DetTo"]);
  NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
  NumericVector state_Ftime      = as<NumericVector>(mod["state_Ftime"]);
  NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
  NumericVector stanzaBasePred   = as<NumericVector>(mod["stanzaBasePred"]);
  NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
  NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
  NumericVector TotGain          = as<NumericVector>(mod["TotGain"]);
  NumericVector TotLoss          = as<NumericVector>(mod["TotLoss"]);
  NumericVector LossPropToB      = as<NumericVector>(mod["LossPropToB"]);
  NumericVector LossPropToQ      = as<NumericVector>(mod["LossPropToQ"]);
  NumericVector DerivT           = as<NumericVector>(mod["DerivT"]);
  NumericVector dyt              = as<NumericVector>(mod["dyt"]);
  NumericVector biomeq           = as<NumericVector>(mod["biomeq"]);
  NumericVector FoodGain         = as<NumericVector>(mod["FoodGain"]);
  NumericVector DetritalGain     = as<NumericVector>(mod["DetritalGain"]);     
  NumericVector FoodLoss         = as<NumericVector>(mod["FoodLoss"]);
  NumericVector UnAssimLoss      = as<NumericVector>(mod["UnAssimLoss"]);
  NumericVector ActiveRespLoss   = as<NumericVector>(mod["ActiveRespLoss"]);   
  NumericVector FishingGain      = as<NumericVector>(mod["FishingGain"]);
  NumericVector MzeroLoss        = as<NumericVector>(mod["MzeroLoss"]);
  NumericVector FishingLoss      = as<NumericVector>(mod["FishingLoss"]);
  NumericVector DetritalLoss     = as<NumericVector>(mod["DetritalLoss"]);
  NumericVector FishingThru      = as<NumericVector>(mod["FishingThru"]);
  NumericVector PredSuite        = as<NumericVector>(mod["PredSuite"]);
  NumericVector HandleSuite      = as<NumericVector>(mod["HandleSuite"]);
  NumericVector preyYY           = as<NumericVector>(mod["preyYY"]);
  NumericVector predYY           = as<NumericVector>(mod["predYY"]);
  NumericVector TerminalF        = as<NumericVector>(mod["TerminalF"]);
  NumericVector TARGET_BIO       = as<NumericVector>(mod["TARGET_BIO"]);
  NumericVector TARGET_F         = as<NumericVector>(mod["TARGET_F"]);
  NumericVector ALPHA            = as<NumericVector>(mod["ALPHA"]);
  
  NumericMatrix force_byprey     = as<NumericMatrix>(mod["force_byprey"]);
  NumericMatrix force_bymort     = as<NumericMatrix>(mod["force_bymort"]);
  NumericMatrix force_bysearch   = as<NumericMatrix>(mod["force_bysearch"]);
  NumericMatrix FORCED_FRATE     = as<NumericMatrix>(mod["FORCED_FRATE"]);
  NumericMatrix FORCED_CATCH     = as<NumericMatrix>(mod["FORCED_CATCH"]);
  
  // Some derivative parts need to be set to zero
  for(i = 0; i < NUM_GROUPS + 1; i++){
    FoodLoss[i]       = 0;
    FoodGain[i]       = 0;
    UnAssimLoss[i]    = 0;
    ActiveRespLoss[i] = 0;   
    DetritalGain[i]   = 0;
    FishingGain[i]    = 0;
    MzeroLoss[i]      = 0;
    FishingLoss[i]    = 0;
    DetritalLoss[i]   = 0;
    FishingThru[i]    = 0;
    PredSuite[i]      = 0;
    HandleSuite[i]    = 0;
  }
			     			  
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
            predYY[sp] *= force_bysearch(y * STEPS_PER_YEAR + m, sp); 
        }
        
 	  // Summed predator and prey for joint handling time and/or scramble functional response
 	     for (links=1; links<=NumPredPreyLinks; links++){
 		          prey = PreyFrom[links];
 		          pred = PreyTo[links];
               PredSuite[prey]   += predYY[pred] * PredPredWeight[links];
               HandleSuite[pred] += preyYY[prey] * PreyPreyWeight[links];
 		  }


	 // Main loop to calculate functional response for each predator/prey link
      for (links=1; links<=NumPredPreyLinks; links++){
 		      prey = PreyFrom[links];
 		      pred = PreyTo[links];
 
      // MAIN FUNCTIONAL RESPONSE SET HERE TO OUTPUT TOTAL CONSUMPTION (Q); use 1 version only
         // (1) This is EwE Classic (c)(r)(tm)(4.0beta)(all rights reserved)
 				 //  Master_Density = predYY[pred];	 
     		 //   Q = QQ[links] * XX[links] * predYY[pred] * preyYY[prey] / 
         //       (Master_Density + ((XX[links] - 1.0) * 
 				 //	                   (1.0 + state_Ftime[prey]) / 2.0) );
  
         // (2) Multiplicative version which is too weird for me.
         // Q =   QQ[links] * 
         //     ( XX[links] * predYY[pred] / ((XX[links] - 1.0) + predYY[pred])     )*
         //     ( HH[links] * preyYY[prey] / ((HH[links] - 1.0) + preyYY[prey])     )*
 				 // 		( DD[links]                 / ((DD[links] - 1.0) + HandleSuite[pred]) )*
 				 // 		( VV[links]                 / ((VV[links] - 1.0) + PredSuite[prey])   );
         
 				 // (3) Additive version: primary used and published in Aydin (2004) 
            // KYA 3/2/2012 setting "COUPLED" to zero means species are density dependent
            // (based on their own modul) but don't interact otherwise.  This can magically
            // create and destroy energy in the system but makes them act like a set
            // of independent surplus production models for comparison purposes 
         Q =   QQ[links] * predYY[pred] * pow(preyYY[prey], COUPLED * HandleSwitch[links]) *
 				      ( DD[links] / ( DD[links] - 1.0 + 
 						                     pow(HandleSelf[pred] * preyYY[prey]   + 
 																 (1. - HandleSelf[pred]) * HandleSuite[pred],
                                                      COUPLED * HandleSwitch[links])) )*
             ( VV[links] / ( VV[links] - 1.0 + 
                                  ScrambleSelf[pred] * predYY[pred] + 
 						                     (1. - ScrambleSelf[pred]) * PredSuite[prey]) );
        
			 // Include any Forcing by prey   
           Q *= force_byprey(y * STEPS_PER_YEAR + m, prey); 

			 // If model is uncoupled, food loss doesn't change with prey or predator levels.
				if (COUPLED){  FoodLoss[prey]  += Q; }
				           else{  FoodLoss[prey]  += state_BB[prey] * QQ[links]/B_BaseRef[prey]; }
        
			 // Energy Accounting
				  FoodGain[pred]           += Q;
          UnAssimLoss[pred]        += Q * UnassimRespFrac[pred]; 
          ActiveRespLoss[pred]     += Q * ActiveRespFrac[pred];  												 
 		  }
 		
     // Mzero Mortality
     for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
         MzeroLoss[sp] = MzeroMort[sp] * state_BB[sp];
     }	   
 		
		 
		// MOST OF THE FOLLOWING FISHING SPECIFICATION METHODS ARE NOT SUPPORTED
		// BY THE R-CODE, only fishing by effort (for gear) or by F-rate (for
		// species) is supported at the end.
		//
		// BY CURRENT R-CODE.  ONLY    
 		// RFISH for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
    // RFISH This sets EFFORT by time series of gear-target combinations
 		// RFISH 		   for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
    // RFISH if -1 is an input value, uses TERMINAL F (last non-negative F) 		
     //RFISH 		   for (gr=NUM_LIVING+NUM_DEAD+1; gr<=NUM_GROUPS; gr++){
     //RFISH 		       if (y+m+d == 0){fish_Effort[gr]=1.0;}
     //RFISH 		       else           {fish_Effort[gr]=1.0;} // NOTE DEFAULT!  THIS CAN BE CHANGED TO 1.0
     //RFISH        // Added 7/8/08 for forced effort
     //RFISH            if (FORCED_EFFORT[gr][y] > -0.001) 
     //RFISH 					    {fish_Effort[gr]=FORCED_EFFORT[gr][y];}
     //RFISH 
     //RFISH 			     if ((FORCED_TARGET[gr]>0) && (FORCED_CATCH[gr][y]>-EPSILON)){
     //RFISH 			        totQ = 0.0;
     //RFISH 			        sp   = FORCED_TARGET[gr];
     //RFISH 			        for (links=1; links<=NumFishingLinks; links++){
     //RFISH     					    if ((FishingThrough[links] == gr) && 
     //RFISH 		    					    (FishingFrom[links]) == sp){
     //RFISH 				    					totQ += FishingQ[links];
     //RFISH 									}
     //RFISH 						  }
     //RFISH 						  fish_Effort[gr] = FORCED_CATCH[gr][y]/ 
     //RFISH 							                  (totQ * state_BB[sp]);
     //RFISH 							if (FORCED_CATCH[gr][y] >= state_BB[sp])
     //RFISH 							   {fish_Effort[gr] = (1.0-EPSILON)*(state_BB[sp])/ 
     //RFISH 				    			                  (totQ * state_BB[sp]);}
     //RFISH 					 }
     //RFISH 					 // By putting F after catch, Frates override absolute catch
     //RFISH 			     if ((FORCED_FTARGET[gr]>0) && (FORCED_FRATE[gr][y]>-EPSILON)){
     //RFISH 			        totQ = 0.0;
     //RFISH 			        sp   = FORCED_FTARGET[gr];
     //RFISH 			        for (links=1; links<=NumFishingLinks; links++){
     //RFISH     					    if ((FishingThrough[links] == gr) && 
     //RFISH 		    					    (FishingFrom[links]) == sp){
     //RFISH 				    					totQ += FishingQ[links];
     //RFISH 									}
     //RFISH 						  }
     //RFISH 						  fish_Effort[gr] = FORCED_FRATE[gr][y]/totQ;
     //RFISH 							//if (FORCED_CATCH[gr][y] >= state_BB[sp])
     //RFISH 							//   {fish_Effort[gr] = (1.0-EPSILON)*(state_BB[sp])/ 
     //RFISH 				    //			                  (totQ * state_BB[sp]);}
     //RFISH 					 }					 
     //RFISH 					 
     //RFISH 				   //if ((y==0) && (m==0) && (d==0)){
     //RFISH 					 //    cout << path_species[gr] << " " << FORCED_TARGET[gr] << " " << path_species[sp] << " " 
     //RFISH 					 //	      << state_BB[sp] << " " << FORCED_CATCH[gr][y] << " "
     //RFISH 					 //		      << fish_Effort[gr] << endl;
     //RFISH 					 //} 					 
     //RFISH 			 }
 			 					 			 					 
   // Apply specified Effort by Gear to catch (using Ecopath-set Q)
      for (links=1; links<=NumFishingLinks; links++){
 				 prey = FishFrom[links];
 				 gr   = FishThrough[links];
 				 dest = FishTo[links];
 				 caught = FishQ[links] * fish_Effort[gr] * state_BB[prey]; 
          FishingLoss[prey] += caught;
          FishingThru[gr]   += caught;
          FishingGain[dest] += caught;
 		}		
 
    //  Special "CLEAN" fisheries assuming q=1, so specified input is Frate
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
             caught = FORCED_CATCH(y, sp) + FORCED_FRATE(y, sp) * state_BB[sp];
             // KYA Aug 2011 removed terminal effort option to allow negative fishing pressure 
                // if (caught <= -EPSILON) {caught = TerminalF[sp] * state_BB[sp];}
             if (caught >= state_BB[sp]){caught = (1.0 - EPSILON) * (state_BB[sp]);}
             FishingLoss[sp] += caught;
             FishingThru[0]  += caught;
             FishingGain[0]  += caught;
             TerminalF[sp] = caught/state_BB[sp];
        }
    
    // KINKED CONTROL RULE - NEEDS INPUT of TARGET BIOMASS and TARGET CATCH
        double RefBio, maxcaught;
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
            if (TARGET_BIO[sp] > EPSILON){
               RefBio    = state_BB[sp] / TARGET_BIO[sp];
               maxcaught = TARGET_F[sp] * state_BB[sp];         
               if      (RefBio > 1.0)           {caught = maxcaught;}
               else if (RefBio >= ALPHA[sp]) {caught = maxcaught * (RefBio - ALPHA[sp])/(1.0 - ALPHA[sp]);}
               else                             {caught = 0.0;}
               FishingLoss[sp] += caught;
               FishingThru[0]  += caught;
               FishingGain[0]  += caught;
               TerminalF[sp] = caught/state_BB[sp];
            }          
        }
        
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
     for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
        FoodLoss[i]  *= force_bymort(y * STEPS_PER_YEAR + m, i);
        MzeroLoss[i] *= force_bymort(y * STEPS_PER_YEAR + m, i);
     }
  
	// Sum up derivitive parts; move previous derivative to dyt        
     for (i=1; i<=NUM_DEAD+NUM_LIVING; i++){
         dyt[i]=DerivT[i];
         TotGain[i] = FoodGain[i] + DetritalGain[i] + FishingGain[i];      
         LossPropToQ[i] = UnAssimLoss[i]  + ActiveRespLoss[i];
         LossPropToB[i] = FoodLoss[i]     + MzeroLoss[i] +
                                FishingLoss[i]  + DetritalLoss[i]; 
                  
         TotLoss[i] = LossPropToQ[i] + LossPropToB[i];    
         DerivT[i]  = TotGain[i] - TotLoss[i]; 
         
 	   // Set biomeq for "fast equilibrium" of fast variables
  	    if (state_BB[i] > 0) {
            biomeq[i] = TotGain[i] / 
 					                 (TotLoss[i] / state_BB[i]);
        }          
     }     
return 0;
}



// SplitSetPred function called in sim stanza initialize and update
// This function simply sums up across juvenile and adult age structure to get 
// population-level Biomass, Numbers, and Consumption 
// [[Rcpp::export]]
int SplitSetPred(List mod){

int ageMo, i;  
double Bt, pt, Nt;

  // Parse out List mod
  int juv_N                      = as<int>(mod["juv_N"]);
  NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
  NumericMatrix WageS            = as<NumericMatrix>(mod["WageS"]);
  NumericVector WWa              = as<NumericVector>(mod["WWa"]);
  NumericMatrix NageS            = as<NumericMatrix>(mod["NageS"]);
  NumericVector state_NN         = as<NumericVector>(mod["state_NN"]);
  NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
  NumericVector stanzaGGJuv      = as<NumericVector>(mod["stanzaGGJuv"]);
  NumericVector stanzaGGAdu      = as<NumericVector>(mod["stanzaGGAdu"]);
  NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
  NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
  NumericVector firstMoJuv       = as<NumericVector>(mod["firstMoJuv"]);
  NumericVector lastMoJuv        = as<NumericVector>(mod["lastMoJuv"]);
  NumericVector firstMoAdu       = as<NumericVector>(mod["firstMoAdu"]);
  NumericVector lastMoAdu        = as<NumericVector>(mod["lastMoAdu"]);
  NumericVector FoodGain         = as<NumericVector>(mod["FoodGain"]);
  
  // loop over split pools
     for (i = 1; i<= juv_N; i++){
         Bt = 1e-30;
         pt = 1e-30;
         Nt = 1e-30;
      // loop over juv monthly age classes
         for (ageMo = firstMoJuv[i]; ageMo <= lastMoJuv[i]; ageMo++){
             Bt = Bt + NageS(ageMo, i) * WageS(ageMo, i);
             pt = pt + NageS(ageMo, i) * WWa(ageMo, i);
             Nt = Nt + NageS(ageMo, i);
             
         }
         state_BB[JuvNum[i]]   = Bt;
         stanzaPred[JuvNum[i]] = pt;
         state_NN[JuvNum[i]]   = Nt;
    
      // loop over adult monthly age classes
         Bt = 1e-30;
         pt = 1e-30;
         Nt = 1e-30;       
         for (ageMo = firstMoAdu[i]; ageMo <= lastMoAdu[i]; ageMo++){
             Bt = Bt + NageS(ageMo, i) * WageS(ageMo, i);
             pt = pt + NageS(ageMo, i) * WWa(ageMo, i);
             Nt = Nt + NageS(ageMo, i);
         }
         state_BB[AduNum[i]]   = Bt;
         stanzaPred[AduNum[i]] = pt;
         state_NN[AduNum[i]]   = Nt;    
         stanzaGGJuv[i] =  FoodGain[JuvNum[i]] / stanzaPred[JuvNum[i]]; 
         stanzaGGAdu[i] =  FoodGain[AduNum[i]] / stanzaPred[AduNum[i]]; 
     }  

return(0);
}



// Update juvenile adult or "stanza" age structure during sim run 
// on monthly timesteps (not hardwiring months, but recommended)
// [[Rcpp::export]] 
 int update_stanzas(List mod, int yr, int mon){
 
int ageMo, i;  
double Su, Gf, Nt; 
//double propRepro;

// Parse out List mod
  int juv_N                      = as<int>(mod["juv_N"]);
  NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
  NumericMatrix WageS            = as<NumericMatrix>(mod["WageS"]);
  NumericVector WWa              = as<NumericVector>(mod["WWa"]);
  NumericMatrix NageS            = as<NumericMatrix>(mod["NageS"]);
  NumericVector state_NN         = as<NumericVector>(mod["state_NN"]);
  NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
  NumericVector stanzaGGJuv      = as<NumericVector>(mod["stanzaGGJuv"]);
  NumericVector stanzaGGAdu      = as<NumericVector>(mod["stanzaGGAdu"]);
  NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
  NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
  NumericVector firstMoJuv       = as<NumericVector>(mod["firstMoJuv"]);
  NumericVector lastMoJuv        = as<NumericVector>(mod["lastMoJuv"]);
  NumericVector firstMoAdu       = as<NumericVector>(mod["firstMoAdu"]);
  NumericVector lastMoAdu        = as<NumericVector>(mod["lastMoAdu"]);
  NumericVector FoodGain         = as<NumericVector>(mod["FoodGain"]);
  NumericVector SpawnBio         = as<NumericVector>(mod["SpawnBio"]);
  NumericVector LossPropToB      = as<NumericVector>(mod["LossPropToB"]);
  NumericVector vBM              = as<NumericVector>(mod["vBM"]);
  NumericMatrix SplitAlpha       = as<NumericMatrix>(mod["SplitAlpha"]);
  NumericVector Wmat001          = as<NumericVector>(mod["Wmat001"]);
  NumericVector Wmat50           = as<NumericVector>(mod["Wmat50"]);
  NumericVector Amat001          = as<NumericVector>(mod["Amat001"]);
  NumericVector Amat50           = as<NumericVector>(mod["Amat50"]);
  NumericVector WmatSpread       = as<NumericVector>(mod["WmatSpread"]);
  NumericVector AmatSpread       = as<NumericVector>(mod["AmatSpread"]);
  NumericVector EggsStanza       = as<NumericVector>(mod["EggsStanza"]);
  NumericVector SpawnEnergy      = as<NumericVector>(mod["SpawnEnergy"]);
  NumericVector SpawnX           = as<NumericVector>(mod["SpawnX"]);
  NumericVector baseSpawnBio     = as<NumericVector>(mod["baseSpawnBio"]);
  NumericVector RecMonth         = as<NumericVector>(mod["RecMonth"]);
  NumericVector baseEggsStanza   = as<NumericVector>(mod["baseEggsStanza"]);
  NumericVector RscaleSplit      = as<NumericVector>(mod["RscaleSplit"]);
  NumericVector RzeroS           = as<NumericVector>(mod["RzeroS"]);
  NumericVector RecPower         = as<NumericVector>(mod["RecPower"]);
  NumericVector VonBD            = as<NumericVector>(mod["VonBD"]);
  NumericMatrix force_byrecs     = as<NumericMatrix>(mod["force_byrecs"]);
   
  // loop over split species groups to update n, wt, biomass in sim  
     for (i = 1; i <= juv_N; i++){
         SpawnBio[i] = 0;
      // loop over juv timesteps with juv survival and growth, also calculate Spawning Biomass
         Su = exp(-LossPropToB[JuvNum[i]] / STEPS_PER_YEAR / state_BB[JuvNum[i]]);    
         Gf = FoodGain[JuvNum[i]] / stanzaPred[JuvNum[i]];  // BB should be stanzaPred            
   			 stanzaGGJuv[i] = Gf; 						         
         for (ageMo = firstMoJuv[i]; ageMo <= lastMoJuv[i]; ageMo++){   
             NageS(ageMo, i) = NageS(ageMo, i) * Su;
             WageS(ageMo, i) = vBM[i] * WageS(ageMo, i) + Gf * SplitAlpha(ageMo, i);
             if ( (WageS(ageMo, i) > Wmat001[i]) && (ageMo > Amat001[i])) {
                 SpawnBio[i]  += WageS(ageMo, i) * NageS(ageMo, i) /
                                 (1. + exp(- ((WageS(ageMo, i) - Wmat50[i]) / WmatSpread[i]) 
 																		      - ((ageMo - Amat50[i]) /AmatSpread[i])));                             
             }            
         }
 				      
      // loop over adult timesteps with adult survival and growth, also calculate Spawning Biomass
         Su = exp(-LossPropToB[AduNum[i]] / STEPS_PER_YEAR / state_BB[AduNum[i]]);
         Gf = FoodGain[AduNum[i]] / stanzaPred[AduNum[i]];   //BB should be stanzaPred         
         stanzaGGAdu[i] = Gf; 
         for (ageMo = firstMoAdu[i]; ageMo <= lastMoAdu[i]; ageMo++){
             NageS(ageMo, i) = NageS(ageMo, i) * Su;
             WageS(ageMo, i) = vBM[i] * WageS(ageMo, i) + Gf * SplitAlpha(ageMo, i);
 						if ( (WageS(ageMo, i) > Wmat001[i]) && (ageMo > Amat001[i])) {
                 SpawnBio[i]  += WageS(ageMo, i) * NageS(ageMo, i) /
                                 (1. + exp(- ((WageS(ageMo, i) - Wmat50[i]) / WmatSpread[i]) 
 																		      - ((ageMo - Amat50[i]) / AmatSpread[i])));                             
             }     	       
         }
         
      // KYA This is a Beverton-Holt curve between Spawning Biomass and
      // Number of Eggs Produced; generally near-linear with default parameters        
         EggsStanza[i] = SpawnBio[i] * SpawnEnergy[i] * SpawnX[i] /
 				                             (SpawnX[i] - 1.0 + 
 			  															 (SpawnBio[i] / baseSpawnBio[i]));
 		 // Forcing for eggs applied here															 
          EggsStanza[i] *= force_byrecs(yr * STEPS_PER_YEAR + mon, JuvNum[i]);		

		 // If recruitment happens all at once (in a single month) that happens here)
     // Otherwise, recruitment is throughout the year.
 			  if (RecMonth[i] > 0){
 				   if (mon == RecMonth[i]) {EggsStanza[i] *=1. ;}
 				 else                         {EggsStanza[i] *=0. ;}
 				}
     
     // now update n and wt looping backward over age
         Nt = NageS(lastMoAdu[i], i) + NageS(lastMoAdu[i] - 1, i);
         if (Nt == 0){Nt = 1e-30;}
         WageS(lastMoAdu[i], i) = (WageS(lastMoAdu[i], i) * NageS(lastMoAdu[i], i) + 
                            WageS(lastMoAdu[i] - 1, i) * NageS(lastMoAdu[i] - 1, i)) / Nt;
         NageS(lastMoAdu[i], i) = Nt;
         
         for (ageMo = lastMoAdu[i] - 1; ageMo > firstMoJuv[i]; ageMo--) {
             NageS(ageMo, i) = NageS(ageMo - 1, i);
             WageS(ageMo, i) = WageS(ageMo - 1, i);   
         }
 
      // finally apply number of eggs to youngest juvenile slot. including Carl's recruit power        
         if (baseEggsStanza[i] > 0){
             NageS(firstMoJuv[i], i) = RscaleSplit[i] * RzeroS[i] * 
                                       pow(double(EggsStanza[i] / baseEggsStanza[i]), 
                                       double(RecPower[i]));
         }
         WageS(firstMoJuv[i], i) = 0;
 
     // Uses generalized vonB (exponent is d).    
     // ADDED FOR STABILITY 4/13/07 (Unlucky Friday)
        for (ageMo = 0; ageMo <= lastMoAdu[i]; ageMo++){
             WWa(ageMo, i) = pow(double(WageS(ageMo, i)), double(VonBD[i])); 
         }           
     }  // end of Juvenile Loop  

return(0); 
}


// [[Rcpp::export]] 
int Adams_Basforth (List mod, int StartYear, int EndYear){
     
int y, m, dd;//c, j, 
int sp, i; //t, ageMo, s, link, prey,pred,links,gr
//int inbound;
double old_B, new_B, pd; //nn, ww, bb, 
//double newbio; float ystep;
//double bioanom, bioratio;
//unsigned int LL;

// Parse out List mod
  int NUM_LIVING = as<int>(mod["NUM_LIVING"]);
  int NUM_DEAD   = as<int>(mod["NUM_DEAD"]);
  int juv_N      = as<int>(mod["juv_N"]);
  int init_run   = as<int>(mod["init_run"]);
  int BURN_YEARS = as<int>(mod["BURN_YEARS"]);
  int CRASH_YEAR = as<int>(mod["CRASH_YEAR"]);
  
  NumericVector B_BaseRef        = as<NumericVector>(mod["B_BaseRef"]);
  NumericVector NoIntegrate      = as<NumericVector>(mod["NoIntegrate"]);
  NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
  NumericVector state_Ftime      = as<NumericVector>(mod["state_Ftime"]);
  NumericVector FtimeAdj         = as<NumericVector>(mod["FtimeAdj"]);
  NumericVector FtimeQBOpt       = as<NumericVector>(mod["FtimeQBOpt"]);
  NumericMatrix WageS            = as<NumericMatrix>(mod["WageS"]);
  NumericMatrix NageS            = as<NumericMatrix>(mod["NageS"]);
  NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
  NumericVector SpawnBio         = as<NumericVector>(mod["SpawnBio"]);
  NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
  NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
  NumericVector firstMoAdu       = as<NumericVector>(mod["firstMoAdu"]);
  NumericVector DerivT           = as<NumericVector>(mod["DerivT"]);
  NumericVector dyt              = as<NumericVector>(mod["dyt"]);
  NumericVector biomeq           = as<NumericVector>(mod["biomeq"]);
  NumericVector FoodGain         = as<NumericVector>(mod["FoodGain"]);
  NumericVector FishingLoss      = as<NumericVector>(mod["FishingLoss"]);

  NumericMatrix out_BB           = as<NumericMatrix>(mod["out_BB"]);
  NumericMatrix out_CC           = as<NumericMatrix>(mod["out_CC"]);
  NumericMatrix out_SSB          = as<NumericMatrix>(mod["out_SSB"]);
  NumericMatrix out_rec          = as<NumericMatrix>(mod["out_rec"]);
  
   // Update sums of split groups to total biomass for derivative calcs
      SplitSetPred(mod); 
      
    // If not starting from previous timestep, call derivative twice to
 	 // set deriv and deriv(t-1)
      if (init_run){
         deriv_master(mod, 0, 0, 0);
         deriv_master(mod, 0, 0, 0);
      }

    // MAIN LOOP STARTS HERE     
       for (y = StartYear; y < EndYear; y++){                                  
         for (m = 0; m < STEPS_PER_YEAR; m++){
                     
           // dd is index for saving monthly output
				      dd = y * STEPS_PER_YEAR + m;
           
           // save current state to output matrix
           // For non-split pools, SSB is output output as the same as B.
					 // For adult split pools, it is overwritten with "actual" SSB, 
					 // for juvs it is 0. 
              for (sp = 1; sp <= NUM_LIVING + NUM_DEAD; sp++){ 
                  out_BB(dd, sp)  = state_BB[sp];
                  out_SSB(dd, sp) = state_BB[sp];
                  out_rec(dd, sp) = state_BB[sp];
              }
              for (i = 1; i <= juv_N; i++){
                  out_SSB(dd, JuvNum[i]) = 0.0;
									out_SSB(dd, AduNum[i]) = SpawnBio[i];
									out_rec(dd, AduNum[i]) = NageS(firstMoAdu[i], i) * WageS(firstMoAdu[i], i);
              }
       
			    // note: code for sub-monthly steps has been removed    
 		      // for (d=0; d<STEPS_PER_MONTH; d++){

 		          // Calculate Derivative for a given timestep
 	               deriv_master(mod, y, m, 0);

 	            // Loop through species, applying Adams-Basforth or fast
 	            // equilibrium depending on species.
 	            
 								 for (sp=1; sp <= NUM_LIVING + NUM_DEAD; sp++){                    

                  // Adjust feeding time 
 	                   if (NoIntegrate[sp] < 0){pd = stanzaPred[sp];}
 	                   else                         {pd = state_BB[sp];}	                      					       
                     if ((pd > 0) && (FoodGain[sp] > 0)){
 										      state_Ftime[sp] = 0.1 + 0.9 * state_Ftime[sp] * 
                                 ((1.0 - FtimeAdj[sp] / (double)STEPS_PER_MONTH) 
 																+ FtimeAdj[sp] / (double)STEPS_PER_MONTH * 
 										             FtimeQBOpt[sp] / (FoodGain[sp] / pd));
                          }	
                     MAX_THRESHOLD(state_Ftime[sp], 2.0);
                      
                  // Biomass update for non-split groups
 								     old_B = state_BB[sp];
 								     new_B = state_BB[sp];
 								     if (NoIntegrate[sp] == 0){
 										 // 'Fast equilibrium'
												new_B        = (1.0 - SORWT) * biomeq[sp] 
 												               + SORWT * state_BB[sp];
 										 }
 										 else if (NoIntegrate[sp] == sp){
                         // Adams-Basforth timestep 
                            new_B = state_BB[sp] + (DELTA_T / 2.0) * 
                                   (3.0 * (DerivT[sp]) - dyt[sp]); 
 										       }
 										       
						      // If the new biomass goes to infinity or something, set a
									// flag to exit the loop and return the problem. 
						         if (std::isnan(new_B) || std::isinf(new_B)) {
                        CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;
                     }
                  
									// If the run is during the "burn-in" years, and biomass goes
									// into the discard range, set flag to exit the loop.  
                     if (y < BURN_YEARS){
                       if  ( (new_B < B_BaseRef[sp] * LO_DISCARD) ||
                             (new_B > B_BaseRef[sp] * HI_DISCARD) ){
                           CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;                         
                       }                          
                     }
                    
                 // Otherwise if biomass goes very large or very small, set to 
                 // min or max.
								    MIN_THRESHOLD(new_B, B_BaseRef[sp] * EPSILON);
                    MAX_THRESHOLD(new_B, B_BaseRef[sp] * BIGNUM);                    
                    
 										state_BB[sp] = new_B;                                
 						        
                 // sum up Catch at every time step (catch is cumulative)
                     if (fabs(state_BB[sp]/old_B-1.0) > EPSILON){
 								       out_CC(dd, sp) = ((FishingLoss[sp] * DELTA_T) / old_B) *
 									                       (state_BB[sp] - old_B) /
 										 	    							 log(state_BB[sp] / old_B);
 									  }
 									  else {out_CC(dd, sp) = (FishingLoss[sp] * DELTA_T);
 								    }		  
  									         
 								 } // End of species loop
		    		       
 					//}  // End of days loop
           
          // Monthly Stanza (split pool) update
 					   update_stanzas(mod, y, m + 1);
             SplitSetPred(mod);
           
				 // As above for non-split groups, check for out-of-bounds, numerical
				 // errors, or burn-in model stoppages for juvs and then adults.   
            for (i = 1; i <= juv_N; i++){
                sp = JuvNum[i];
                new_B = state_BB[sp];
						    if (std::isnan(new_B) || std::isinf(new_B)) {
                   CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;
                }                  
                if (y < BURN_YEARS){
                    if  ( (new_B < B_BaseRef[sp] * LO_DISCARD) ||
                          (new_B > B_BaseRef[sp] * HI_DISCARD) ){
                         CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;                           
                    }                          
                }                
                sp    = AduNum[i];
                new_B = state_BB[sp];
						    if (std::isnan(new_B) || std::isinf(new_B)) {
                    CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;
                }     
                if (y < BURN_YEARS){
                   if  ( (new_B < B_BaseRef[sp] * LO_DISCARD) ||
                         (new_B > B_BaseRef[sp] * HI_DISCARD) ){
                       CRASH_YEAR = y; y = EndYear; m = STEPS_PER_YEAR;                           
                   }                          
                }                
            }
                           
 		    }  // End of main months loop
     
   }// End of years loop
       
   // If you get to the end and haven't crashed, save biomass in the
   // final time slot before exiting.
      if (CRASH_YEAR < 0){
         dd = EndYear * STEPS_PER_YEAR;
         for (sp = 1; sp <= NUM_LIVING + NUM_DEAD; sp++){ 
              out_BB(dd, sp)  = state_BB[sp];
              out_SSB(dd, sp) = state_BB[sp];
              out_rec(dd, sp) = state_BB[sp];
         }
         for (i = 1; i <= juv_N; i++){
              out_SSB(dd, JuvNum[i]) = 0.0;
						  out_SSB(dd, AduNum[i]) = SpawnBio[i];
						  out_rec(dd, AduNum[i]) = NageS(firstMoAdu[i], i) * WageS(firstMoAdu[i], i);
         }
      }  
      
return(0);
} 
