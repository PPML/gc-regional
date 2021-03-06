// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// mv_mult
arma::vec mv_mult(arma::mat& lhs, arma::vec& rhs);
RcppExport SEXP _gcRegional_mv_mult(SEXP lhsSEXP, SEXP rhsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rhs(rhsSEXP);
    rcpp_result_gen = Rcpp::wrap(mv_mult(lhs, rhs));
    return rcpp_result_gen;
END_RCPP
}
// mat_mult
arma::mat mat_mult(arma::mat& a, arma::mat& b);
RcppExport SEXP _gcRegional_mat_mult(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_mult(a, b));
    return rcpp_result_gen;
END_RCPP
}
// gcSim
Rcpp::List gcSim(List x, double dt, List cm, int nYrs, NumericVector initPop, NumericVector n_sa, arma::mat births, arma::mat births_sa, arma::mat births_nsa, arma::mat aging, arma::mat aging_nsa, bool debug, bool quarterly, bool increase_risk_all_pops);
RcppExport SEXP _gcRegional_gcSim(SEXP xSEXP, SEXP dtSEXP, SEXP cmSEXP, SEXP nYrsSEXP, SEXP initPopSEXP, SEXP n_saSEXP, SEXP birthsSEXP, SEXP births_saSEXP, SEXP births_nsaSEXP, SEXP agingSEXP, SEXP aging_nsaSEXP, SEXP debugSEXP, SEXP quarterlySEXP, SEXP increase_risk_all_popsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< List >::type cm(cmSEXP);
    Rcpp::traits::input_parameter< int >::type nYrs(nYrsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initPop(initPopSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_sa(n_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births(birthsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_sa(births_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_nsa(births_nsaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging(agingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging_nsa(aging_nsaSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< bool >::type quarterly(quarterlySEXP);
    Rcpp::traits::input_parameter< bool >::type increase_risk_all_pops(increase_risk_all_popsSEXP);
    rcpp_result_gen = Rcpp::wrap(gcSim(x, dt, cm, nYrs, initPop, n_sa, births, births_sa, births_nsa, aging, aging_nsa, debug, quarterly, increase_risk_all_pops));
    return rcpp_result_gen;
END_RCPP
}
// gcSim2
Rcpp::List gcSim2(List x, double dt, List cm, int nYrs, NumericVector initPop, NumericVector n_sa, arma::mat births, arma::mat births_sa, arma::mat births_nsa, arma::mat aging, arma::mat aging_nsa, bool debug, bool quarterly, bool increase_risk_all_pops);
RcppExport SEXP _gcRegional_gcSim2(SEXP xSEXP, SEXP dtSEXP, SEXP cmSEXP, SEXP nYrsSEXP, SEXP initPopSEXP, SEXP n_saSEXP, SEXP birthsSEXP, SEXP births_saSEXP, SEXP births_nsaSEXP, SEXP agingSEXP, SEXP aging_nsaSEXP, SEXP debugSEXP, SEXP quarterlySEXP, SEXP increase_risk_all_popsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< List >::type cm(cmSEXP);
    Rcpp::traits::input_parameter< int >::type nYrs(nYrsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initPop(initPopSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_sa(n_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births(birthsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_sa(births_saSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type births_nsa(births_nsaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging(agingSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aging_nsa(aging_nsaSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< bool >::type quarterly(quarterlySEXP);
    Rcpp::traits::input_parameter< bool >::type increase_risk_all_pops(increase_risk_all_popsSEXP);
    rcpp_result_gen = Rcpp::wrap(gcSim2(x, dt, cm, nYrs, initPop, n_sa, births, births_sa, births_nsa, aging, aging_nsa, debug, quarterly, increase_risk_all_pops));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gcRegional_mv_mult", (DL_FUNC) &_gcRegional_mv_mult, 2},
    {"_gcRegional_mat_mult", (DL_FUNC) &_gcRegional_mat_mult, 2},
    {"_gcRegional_gcSim", (DL_FUNC) &_gcRegional_gcSim, 14},
    {"_gcRegional_gcSim2", (DL_FUNC) &_gcRegional_gcSim2, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_gcRegional(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
