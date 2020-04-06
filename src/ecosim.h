#include <Rcpp.h>
using namespace Rcpp;
 
// Various Constants
#define MAX_MONTHS_STANZA  400      // Max number age pools within stanza calcs (months)
//#define LO_DISCARD         1e-4     // Discard Limit for sense - KYA 2020 made this an rsim.sense input
//#define HI_DISCARD         1e+4     // Discard Limit for sense - KYA 2020 made this an rsim.sense input
#define STEPS_PER_YEAR     12                   // Steps per year between stanzas ("months")
#define DELTA_T            0.083333333333333333 // 1/12 for steps per month
#define SORWT              0.5                  // Fast Equilibrium smoother
// KYA April 2020 - changed below from (1e-8, 1e8) to (1e-15, 1e15) for greater precision
// this may change recovery times from "0" biomass - sufficiently consistent from old?
#define EPSILON            1.0E-15                 // Test threshold for "too close to zero"
#define BIGNUM             1.0E+15                 // Test threshold for "too big"

#define MIN_THRESHOLD(x,y) if ((x)<(y))(x)=(y)  // Set x to min y if too small
#define MAX_THRESHOLD(x,y) if ((x)>(y))(x)=(y)  // Set x to max y if too big

// Functions in ecosim.cpp 
   List deriv_vector(List params, List state, List forcing, List fishing, List stanzas, 
                     int y, int m, double tt);
   //List Adams_run(List params, List instate, List forcing, List fishing, int StartYear, 
  //                int EndYear);
   List rk4_run(List params, List instate, List forcing, List fishing, int StartYear, 
                int EndYear);

   int deriv_old(List mod, int y, int m, int d);
   int SplitSetPred(List stanzas, List state);
   int SplitUpdate(List stanzas, List state, List forcing, List deriv, int yr, int mon);
   int Adams_Basforth_old (List mod, int StartYear, int EndYear);

// Power function that works with NumericVectors
struct pow_wrapper{
   public: double operator()(double a, double b){
      return ::pow(a, b);
   }
};

NumericVector vpow(const NumericVector base, const NumericVector exp) {
   NumericVector out(base.size());
   std::transform(base.cbegin(), base.cend(), exp.cbegin(), out.begin(), pow_wrapper());
   return out;
}
         