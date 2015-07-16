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

List deriv_test(List par, List force, int y, int m, int d);
int deriv_master(List mod, int y, int m, int d);
int SplitSetPred(List mod);
int update_stanzas(List mod, int yr, int mon);
int Adams_Basforth (List mod, int StartYear, int EndYear);

// [[Rcpp::export]]
NumericVector vpow(const NumericVector base, const NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow); return out; }
         