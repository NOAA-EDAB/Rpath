// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rk4_run
List rk4_run(List params, List instate, List forcing, List fishing, List stanzas, int StartYear, int EndYear);
RcppExport SEXP _Rpath_rk4_run(SEXP paramsSEXP, SEXP instateSEXP, SEXP forcingSEXP, SEXP fishingSEXP, SEXP stanzasSEXP, SEXP StartYearSEXP, SEXP EndYearSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type instate(instateSEXP);
    Rcpp::traits::input_parameter< List >::type forcing(forcingSEXP);
    Rcpp::traits::input_parameter< List >::type fishing(fishingSEXP);
    Rcpp::traits::input_parameter< List >::type stanzas(stanzasSEXP);
    Rcpp::traits::input_parameter< int >::type StartYear(StartYearSEXP);
    Rcpp::traits::input_parameter< int >::type EndYear(EndYearSEXP);
    rcpp_result_gen = Rcpp::wrap(rk4_run(params, instate, forcing, fishing, stanzas, StartYear, EndYear));
    return rcpp_result_gen;
END_RCPP
}
// Adams_run
List Adams_run(List params, List instate, List forcing, List fishing, List stanzas, int StartYear, int EndYear, List InitDeriv);
RcppExport SEXP _Rpath_Adams_run(SEXP paramsSEXP, SEXP instateSEXP, SEXP forcingSEXP, SEXP fishingSEXP, SEXP stanzasSEXP, SEXP StartYearSEXP, SEXP EndYearSEXP, SEXP InitDerivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type instate(instateSEXP);
    Rcpp::traits::input_parameter< List >::type forcing(forcingSEXP);
    Rcpp::traits::input_parameter< List >::type fishing(fishingSEXP);
    Rcpp::traits::input_parameter< List >::type stanzas(stanzasSEXP);
    Rcpp::traits::input_parameter< int >::type StartYear(StartYearSEXP);
    Rcpp::traits::input_parameter< int >::type EndYear(EndYearSEXP);
    Rcpp::traits::input_parameter< List >::type InitDeriv(InitDerivSEXP);
    rcpp_result_gen = Rcpp::wrap(Adams_run(params, instate, forcing, fishing, stanzas, StartYear, EndYear, InitDeriv));
    return rcpp_result_gen;
END_RCPP
}
// deriv_vector
List deriv_vector(List params, List state, List forcing, List fishing, List stanzas, int inyear, int m, double tt);
RcppExport SEXP _Rpath_deriv_vector(SEXP paramsSEXP, SEXP stateSEXP, SEXP forcingSEXP, SEXP fishingSEXP, SEXP stanzasSEXP, SEXP inyearSEXP, SEXP mSEXP, SEXP ttSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< List >::type forcing(forcingSEXP);
    Rcpp::traits::input_parameter< List >::type fishing(fishingSEXP);
    Rcpp::traits::input_parameter< List >::type stanzas(stanzasSEXP);
    Rcpp::traits::input_parameter< int >::type inyear(inyearSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type tt(ttSEXP);
    rcpp_result_gen = Rcpp::wrap(deriv_vector(params, state, forcing, fishing, stanzas, inyear, m, tt));
    return rcpp_result_gen;
END_RCPP
}
// SplitSetPred
int SplitSetPred(List stanzas, List state);
RcppExport SEXP _Rpath_SplitSetPred(SEXP stanzasSEXP, SEXP stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type stanzas(stanzasSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    rcpp_result_gen = Rcpp::wrap(SplitSetPred(stanzas, state));
    return rcpp_result_gen;
END_RCPP
}
// SplitUpdate
int SplitUpdate(List stanzas, List state, List forcing, List deriv, int yr, int mon);
RcppExport SEXP _Rpath_SplitUpdate(SEXP stanzasSEXP, SEXP stateSEXP, SEXP forcingSEXP, SEXP derivSEXP, SEXP yrSEXP, SEXP monSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type stanzas(stanzasSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    Rcpp::traits::input_parameter< List >::type forcing(forcingSEXP);
    Rcpp::traits::input_parameter< List >::type deriv(derivSEXP);
    Rcpp::traits::input_parameter< int >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< int >::type mon(monSEXP);
    rcpp_result_gen = Rcpp::wrap(SplitUpdate(stanzas, state, forcing, deriv, yr, mon));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rpath_rk4_run", (DL_FUNC) &_Rpath_rk4_run, 7},
    {"_Rpath_Adams_run", (DL_FUNC) &_Rpath_Adams_run, 8},
    {"_Rpath_deriv_vector", (DL_FUNC) &_Rpath_deriv_vector, 8},
    {"_Rpath_SplitSetPred", (DL_FUNC) &_Rpath_SplitSetPred, 2},
    {"_Rpath_SplitUpdate", (DL_FUNC) &_Rpath_SplitUpdate, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rpath(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
