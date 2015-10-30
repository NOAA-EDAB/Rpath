#include <Rcpp.h>
using namespace Rcpp;
 
// Various Constants
#define MAX_MONTHS_STANZA  400      // Max number age pools within stanza calcs (months)
#define LO_DISCARD         1e-4     // Discard Limit for sense
#define HI_DISCARD         1e+4     // Discard Limit for sense
#define STEPS_PER_YEAR     12                   // Steps per year between stanzas ("months")
#define DELTA_T            0.083333333333333333 // 1/12 for steps per month
#define SORWT              0.5                  // Fast Equilibrium smoother
#define EPSILON            1E-8                 // Test threshold for "too close to zero"
#define BIGNUM             1E+8                 // Test threshold for "too big"

#define MIN_THRESHOLD(x,y) if ((x)<(y))(x)=(y)  // Set x to min y if too small
#define MAX_THRESHOLD(x,y) if ((x)>(y))(x)=(y)  // Set x to max y if too big

// Functions in ecosim.cpp 
   List deriv_vector(List params, List state, List forcing, List fishing, List stanzas, 
                     int y, int m, double tt);
   List Adams_run(List params, List instate, List forcing, List fishing, int StartYear, 
                  int EndYear);
   List rk4_run(List params, List instate, List forcing, List fishing, int StartYear, 
                int EndYear);

   int deriv_old(List mod, int y, int m, int d);
   int SplitSetPred(List stanzas);
   int update_stanzas(List mod, int yr, int mon);
   int Adams_Basforth_old (List mod, int StartYear, int EndYear);

// Power function that works with NumericVectors
   NumericVector vpow(const NumericVector base, const NumericVector exp) {
   NumericVector out(base.size());
   std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow); return out; }
         