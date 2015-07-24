// Note: this is compiled in C++ to interface with R declarations (extern C) but
// the code is pretty much straight c.

#include <math.h>
#include <stdlib.h>
#include <string.h>
using namespace std;             

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

// newsim is giant structure of all variables to be brought over from R.
// It is (mostly) pointers to set pointing to R vectors

 struct newsim{  
  // Dim Variables
     int YEARS            ;
     int NUM_GROUPS       ;//   1   -none- numeric  
     int NUM_LIVING       ;//   1   -none- numeric  
     int NUM_DEAD         ;//   1   -none- numeric  
     int NUM_GEARS        ;//   1   -none- numeric  

     int BURN_YEARS       ;
     int COUPLED          ;
     int CRASH_YEAR       ;
  // flags
     int init_run;
//  double *spname           ;//  35   -none- character
//  double *spnum            ;//  35   -none- numeric  
  double *B_BaseRef        ;//  35   -none- numeric  
  double *MzeroMort        ;//  35   -none- numeric  
  double *UnassimRespFrac  ;//  35   -none- numeric  
  double *ActiveRespFrac   ;//  35   -none- numeric  
  double *FtimeAdj         ;//  35   -none- numeric  
  double *FtimeQBOpt       ;//  35   -none- numeric  
  double *PBopt            ;//  35   -none- numeric  
  double *NoIntegrate      ;//  35   -none- numeric  
  double *HandleSelf       ;//  35   -none- numeric  
  double *ScrambleSelf     ;//  35   -none- numeric  
  double *PredTotWeight    ;//  35   -none- numeric  
  double *PreyTotWeight    ;//  35   -none- numeric  

  int NumPredPreyLinks ;//   1   -none- numeric 
  int *PreyFrom         ;// 175   -none- numeric  
  int *PreyTo           ;// 175   -none- numeric  
  double *QQ               ;// 175   -none- numeric  
  double *DD               ;// 175   -none- numeric  
  double *VV               ;// 175   -none- numeric  
  double *HandleSwitch     ;// 175   -none- numeric  
  double *PredPredWeight   ;// 175   -none- numeric  
  double *PreyPreyWeight   ;// 175   -none- numeric  

  int NumFishingLinks  ;//   1   -none- numeric   
  int *FishingFrom         ;//  89   -none- numeric  
  int *FishingThrough      ;//  89   -none- numeric  
  int *FishingTo           ;//  89   -none- numeric 
  double *FishingQ            ;//  89   -none- numeric  
 

  int NumDetLinks      ;//   1   -none- numeric  
  int *DetFrom          ;//  53   -none- numeric  
  int *DetTo            ;//  53   -none- numeric  
  double *DetFrac          ;//  53   -none- numeric  

  int juv_N            ;//   1   -none- numeric    
  int *JuvNum           ;//   5   -none- numeric  
  int *AduNum           ;//   5   -none- numeric 
  int *firstMoJuv       ;//   5   -none- numeric  
  int *lastMoJuv        ;//   5   -none- numeric  
  int *firstMoAdu       ;//   5   -none- numeric  
  int *lastMoAdu        ;//   5   -none- numeric  
  double *SpawnBio         ;//   5   -none- numeric  
  double *EggsStanza       ;//   5   -none- numeric  
  double *stanzaGGJuv      ;//   5   -none- numeric  
  double *stanzaGGAdu      ;//   5   -none- numeric  
  double *SpawnEnergy      ;//   5   -none- numeric  
  double *SpawnX           ;//   5   -none- numeric  
  double *SpawnAllocR      ;//   5   -none- numeric  
  double *SpawnAllocG      ;//   5   -none- numeric  
  double *recruits         ;//   5   -none- numeric  
  double *RzeroS           ;//   5   -none- numeric  
  double *baseEggsStanza   ;//   5   -none- numeric  
  double *baseSpawnBio     ;//   5   -none- numeric  
  double *Rbase            ;//   5   -none- numeric  
  double *RecMonth         ;//   5   -none- numeric  
  double *VonBD            ;//   5   -none- numeric  
  double *aduEqAgeZ        ;//   5   -none- numeric  
  double *juvEqAgeZ        ;//   5   -none- numeric  
  double *RecPower         ;//   5   -none- numeric  
  double *Wmat001          ;//   5   -none- numeric  
  double *Wmat50           ;//   5   -none- numeric  
  double *Amat001          ;//   5   -none- numeric  
  double *Amat50           ;//   5   -none- numeric  
  double *WmatSpread       ;//   5   -none- numeric  
  double *AmatSpread       ;//   5   -none- numeric  
  double *vBM              ;//   5   -none- numeric  
  double *RscaleSplit      ;//   5   -none- numeric 

  double *state_BB         ;//  35   -none- numeric  
  double *state_Ftime      ;//  35   -none- numeric  
  double *state_NN         ;//  35   -none- numeric  
  double *stanzaPred       ;//  35   -none- numeric  
  double *stanzaBasePred   ;//  35   -none- numeric  
  double *fish_Effort; //[NUM_GROUPS+1];
  
  double **WageS            ;//2005   -none- numeric  
  double **WWa              ;//2005   -none- numeric  
  double **NageS            ;//2005   -none- numeric  
  double **SplitAlpha       ;//2005   -none- numeric 

// Variables prefaced by deriv_ are derivative parts (e.g. current consumption)
   double *deriv_TotGain;
   double *deriv_TotLoss;
   double *deriv_LossPropToB;
   double *deriv_LossPropToQ;

   double *deriv_DerivT;
   double *deriv_dyt;
   double *deriv_biomeq;

   double *deriv_FoodGain;
   double *deriv_DetritalGain;
   double *deriv_FoodLoss;
   double *deriv_UnAssimLoss;
   double *deriv_ActiveRespLoss;
   //float deriv_PassiveRespLoss[NUM_GROUPS + 1];
   double *deriv_MzeroLoss;
   double *deriv_DetritalLoss;
   double *deriv_TotDetOut;
   
   double *deriv_FishingLoss;
   double *deriv_FishingGain;
   double *deriv_FishingThru;

   double *deriv_preyYY;
   double *deriv_predYY;
   double *deriv_PredSuite;
 	 double *deriv_HandleSuite;
   double *deriv_TerminalF; //[NUM_GROUPS+1];

	    double *TARGET_BIO;	 
	    double *TARGET_FRATE;
	    double *ALPHA;
	    
// TIME SERIES FORCING FUNCTIONS
      double **FORCED_FRATE;   // float FORCED_FRATE[NUM_GROUPS+1][MAX_YEARS+1];
      double **FORCED_CATCH;   // float FORCED_CATCH[NUM_GROUPS+1][MAX_YEARS+1];
     //RFISH      int *FORCED_TARGET;   // int FORCED_TARGET[NUM_GROUPS+1];
	   //RFISH int *FORCED_FTARGET; //int FORCED_FTARGET[NUM_GROUPS+1];
	   //RFISH double **FORCED_EFFORT;  // float FORCED_EFFORT[NUM_GROUPS+1][MAX_YEARS+1]; 

   // Force prey to be consumed more (including prey 0 for PrimProd))
      double **force_byprey; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      double **force_bymort; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      double **force_byrecs; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];
      double **force_bysearch; //[NUM_GROUPS + 1][MAX_YEARS*STEPS_PER_YEAR+1];

     double **out_BB;
     double **out_CC;
     double **out_SSB;
     double **out_rec;
};  // END OF GIANT STRUCTURE

// List of functions in main code
   int deriv_master(struct newsim *v, int y, int m, int d);
   int Adams_Basforth(struct newsim *v, int StartYear, int EndYear);
   int SplitSetPred(struct newsim *v);
   int update_stanzas(struct newsim *v, int yr, int mon);

// This Function transposes matrixes from R input order to get correct vector of vectors 
   double **rmat_trans(int cols, int rows, double *rvec){
   double **lm;
   int r,c;
   lm = (double **) malloc((size_t)((cols)*sizeof(double*)));
   for (c=0; c<cols; c++){lm[c] = rvec + (rows)*c;}
   return(lm);
} 
