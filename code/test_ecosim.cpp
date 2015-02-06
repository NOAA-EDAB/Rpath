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

// Deriv master calculates the biomass dynamics derivative
// [[Rcpp::export]]
int deriv_master(List mod, int y, int m, int d){
 //if (!mod.inherits("Rpath.sim")) stop("Input must be a Rpath model");

  // Functional response vars     
  int sp, links, prey, pred, i;
  double caught, Q;
  unsigned int LL;
  // Fishing vars
  int gr, dest;
  // double Master_Density, totQ;
  
  // Parse out List mod
  int NUM_GROUPS                 = as<int>(mod["NUM_GROUPS"]);
  int NUM_LIVING                 = as<int>(mod["NUM_LIVING"]);
  int NUM_DEAD                   = as<int>(mod["NUM_DEAD"]);
  //int NUM_GEARS                  = as<int>(mod["NUM_GEARS"]);
  int NumPredPreyLinks           = as<int>(mod["NumPredPreyLinks"]);
  int NumFishingLinks            = as<int>(mod["NumFishingLinks"]);
  int NumDetLinks                = as<int>(mod["NumDetLinks"]);
  int juv_N                      = as<int>(mod["juv_N"]);
  int COUPLED                    = as<int>(mod["COUPLED"]);
  
  CharacterVector spname         = as<CharacterVector>(mod["spname"]);
  
  NumericVector spnum            = as<NumericVector>(mod["spnum"]);
  NumericVector B_BaseRef        = as<NumericVector>(mod["B_BaseRef"]);
  NumericVector MzeroMort        = as<NumericVector>(mod["MzeroMort"]);
  NumericVector UnassimRespFrac  = as<NumericVector>(mod["UnassimRespFrac"]);
  NumericVector ActiveRespFrac   = as<NumericVector>(mod["ActiveRespFrac"]);
  NumericVector FtimeAdj         = as<NumericVector>(mod["FtimeAdj"]);
  NumericVector FtimeQBOpt       = as<NumericVector>(mod["FtimeQBOpt"]);
  NumericVector PBopt            = as<NumericVector>(mod["PBopt"]);
  NumericVector fish_Effort      = as<NumericVector>(mod["fish_Effort"]);
  NumericVector NoIntegrate      = as<NumericVector>(mod["NoIntegrate"]);
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
  NumericVector PredTotWeight    = as<NumericVector>(mod["PredTotWeight"]);
  NumericVector PreyTotWeight    = as<NumericVector>(mod["PreyTotWeight"]);
  NumericVector FishFrom         = as<NumericVector>(mod["FishFrom"]);
  NumericVector FishThrough      = as<NumericVector>(mod["FishThrough"]);
  NumericVector FishQ            = as<NumericVector>(mod["FishQ"]);
  NumericVector FishTo           = as<NumericVector>(mod["FishTo"]);
  NumericVector DetFrac          = as<NumericVector>(mod["DetFrac"]);
  NumericVector DetFrom          = as<NumericVector>(mod["DetFrom"]);
  NumericVector DetTo            = as<NumericVector>(mod["DetTo"]);
  NumericVector state_BB         = as<NumericVector>(mod["state_BB"]);
  NumericVector state_Ftime      = as<NumericVector>(mod["state_Ftime"]);
  NumericVector WageS            = as<NumericVector>(mod["WageS"]);
  NumericVector WWa              = as<NumericVector>(mod["WWa"]);
  NumericVector NageS            = as<NumericVector>(mod["NageS"]);
  NumericVector SplitAlpha       = as<NumericVector>(mod["SplitAlpha"]);
  NumericVector state_NN         = as<NumericVector>(mod["state_NN"]);
  NumericVector stanzaPred       = as<NumericVector>(mod["stanzaPred"]);
  NumericVector stanzaBasePred   = as<NumericVector>(mod["stanzaBasePred"]);
  NumericVector SpawnBio         = as<NumericVector>(mod["SpawnBio"]);
  NumericVector EggsStanza       = as<NumericVector>(mod["EggsStanza"]);
  NumericVector stanzaGGJuv      = as<NumericVector>(mod["stanzaGGJuv"]);
  NumericVector stanzaGGAdu      = as<NumericVector>(mod["stanzaGGAdu"]);
  NumericVector SpawnEnergy      = as<NumericVector>(mod["SpawnEnergy"]);
  NumericVector SpawnX           = as<NumericVector>(mod["SpawnX"]);
  NumericVector SpawnAllocR      = as<NumericVector>(mod["SpawnAllocR"]);
  NumericVector SpawnAllocG      = as<NumericVector>(mod["SpawnAllocG"]);
  NumericVector recruits         = as<NumericVector>(mod["recruits"]);
  NumericVector RzeroS           = as<NumericVector>(mod["RzeroS"]);
  NumericVector baseEggsStanza   = as<NumericVector>(mod["baseEggsStanza"]);
  NumericVector baseSpawnBio     = as<NumericVector>(mod["baseSpawnBio"]);
  NumericVector Rbase            = as<NumericVector>(mod["Rbase"]);
  NumericVector RscaleSplit      = as<NumericVector>(mod["RscaleSplit"]);
  NumericVector JuvNum           = as<NumericVector>(mod["JuvNum"]);
  NumericVector AduNum           = as<NumericVector>(mod["AduNum"]);
  NumericVector RecMonth         = as<NumericVector>(mod["RecMonth"]);
  NumericVector VonBD            = as<NumericVector>(mod["VonBD"]);
  NumericVector aduEqAgeZ        = as<NumericVector>(mod["aduEqAgeZ"]);
  NumericVector juvEqAgeZ        = as<NumericVector>(mod["juvEqAgeZ"]);
  NumericVector RecPower         = as<NumericVector>(mod["RecPower"]);
  NumericVector Wmat001          = as<NumericVector>(mod["Wmat001"]);
  NumericVector Wmat50           = as<NumericVector>(mod["Wmat50"]);
  NumericVector Amat001          = as<NumericVector>(mod["Amat001"]);
  NumericVector Amat50           = as<NumericVector>(mod["Amat50"]);
  NumericVector WmatSpread       = as<NumericVector>(mod["WmatSpread"]);
  NumericVector AmatSpread       = as<NumericVector>(mod["AmatSpread"]);
  NumericVector vBM              = as<NumericVector>(mod["vBM"]);
  NumericVector firstMoJuv       = as<NumericVector>(mod["firstMoJuv"]);
  NumericVector lastMoJuv        = as<NumericVector>(mod["lastMoJuv"]);
  NumericVector firstMoAdu       = as<NumericVector>(mod["firstMoAdu"]);
  NumericVector lastMoAdu        = as<NumericVector>(mod["lastMoAdu"]);
  NumericVector YEARS            = as<NumericVector>(mod["YEARS"]);
  NumericVector BURN_YEARS       = as<NumericVector>(mod["BURN_YEARS"]);
  NumericVector CRASH_YEAR       = as<NumericVector>(mod["CRASH_YEAR"]);
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
  NumericVector TotDetOut        = as<NumericVector>(mod["TotDetOut"]);
  NumericVector preyYY           = as<NumericVector>(mod["preyYY"]);
  NumericVector predYY           = as<NumericVector>(mod["predYY"]);
  NumericVector TerminalF        = as<NumericVector>(mod["TerminalF"]);
  NumericVector TARGET_BIO       = as<NumericVector>(mod["TARGET_BIO"]);
  NumericVector TARGET_F         = as<NumericVector>(mod["TARGET_F"]);
  NumericVector ALPHA            = as<NumericVector>(mod["ALPHA"]);
  
  DataFrame force_byprey     = as<DataFrame>(mod["force_byprey"]);
  DataFrame force_bymort     = as<DataFrame>(mod["force_bymort"]);
  DataFrame force_byrecs     = as<DataFrame>(mod["force_byrecs"]);
  DataFrame force_bysearch   = as<DataFrame>(mod["force_bysearch"]);
  DataFrame FORCED_FRATE     = as<DataFrame>(mod["FORCED_FRATE"]);
  DataFrame FORCED_CATCH     = as<DataFrame>(mod["FORCED_CATCH"]);
  DataFrame out_BB           = as<DataFrame>(mod["out_BB"]);
  DataFrame out_CC           = as<DataFrame>(mod["out_CC"]);
  DataFrame out_SSB          = as<DataFrame>(mod["out_SSB"]);
  DataFrame out_rec          = as<DataFrame>(mod["out_rec"]);
    
  
  // Some derivative parts need to be set to zero
//  LL =  (NUM_GROUPS + 1) * sizeof(double);
//  memset(FoodLoss,         0, LL);
//  memset(FoodGain,         0, LL);
//  memset(UnAssimLoss,      0, LL);
//  memset(ActiveRespLoss,   0, LL);   
//  memset(DetritalGain,     0, LL);
//  memset(FishingGain,      0, LL);
//  memset(MzeroLoss,        0, LL);
//  memset(FishingLoss,      0, LL);
//  memset(DetritalLoss,     0, LL);
//  memset(FishingThru,      0, LL);
//  memset(PredSuite,        0, LL);
//  memset(HandleSuite,      0, LL);
// 				     			  
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
            NumericVector sp_force_bysearch = force_bysearch[sp];
            predYY[sp] *= sp_force_bysearch[y * STEPS_PER_YEAR + m]; 
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
 				  NumericVector prey_force_byprey = force_byprey[prey];
           Q *= prey_force_byprey[y * STEPS_PER_YEAR + m]; 

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
             NumericVector sp_FORCED_CATCH = FORCED_CATCH[sp];
             NumericVector sp_FORCED_FRATE = FORCED_FRATE[sp];
             caught = sp_FORCED_CATCH[y] + sp_FORCED_FRATE[y] * state_BB[sp];
             // KYA Aug 2011 removed terminal effort option to allow negative fishing pressure 
                // if (caught <= -EPSILON) {caught = TerminalF[sp] * state_BB[sp];}
             if (caught>=state_BB[sp]){caught=(1.0-EPSILON)*(state_BB[sp]);}
             FishingLoss[sp] += caught;
             FishingThru[0]  += caught;
             FishingGain[0]  += caught;
             TerminalF[sp] = caught/state_BB[sp];
        }
    
    // KINKED CONTROL RULE - NEEDS INPUT of TARGET BIOMASS and TARGET CATCH
        double RefBio, maxcaught;
        for (sp=1; sp<=NUM_LIVING+NUM_DEAD; sp++){
            if (TARGET_BIO[sp] > EPSILON){
               RefBio    = state_BB[sp]/TARGET_BIO[sp];
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
        NumericVector i_force_bymort = force_bymort[i];
        FoodLoss[i]  *= i_force_bymort[y * STEPS_PER_YEAR + m];
        MzeroLoss[i] *= i_force_bymort[y * STEPS_PER_YEAR + m];
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
double propRepro;

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
  DataFrame force_byrecs         = as<DataFrame>(mod["force_byrecs"]);
   
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
 			   NumericVector jv_force_byrecs = force_byrecs[JuvNum[i]];
          EggsStanza[i] *= jv_force_byrecs[yr * STEPS_PER_YEAR + mon];		

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
