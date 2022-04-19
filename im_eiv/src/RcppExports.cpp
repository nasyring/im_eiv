#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


// randsetsMCMC
Rcpp::List plauscontour(NumericVector par, NumericVector stat, NumericVector del, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz);
RcppExport SEXP imeiv_plauscontour(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP, SEXP truebzSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd(propsdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    __result = Rcpp::wrap(plauscontour(par,stat,del,n,propsd,truebx,truebz));
    return __result;
END_RCPP
}

