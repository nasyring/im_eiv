#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

Rcpp::List plauscontourMCMC(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector truebz, NumericVector bzseq, NumericVector randsettype);
RcppExport SEXP imeiv_plauscontourMCMC(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP , SEXP bxseqSEXP, SEXP truebzSEXP, SEXP bzseqSEXP, SEXP randsettypeSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd(propsdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bzseq(bzseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type randsettype(randsettypeSEXP);   
    __result = Rcpp::wrap(plauscontourMCMC(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq));
    return __result;
END_RCPP
}

Rcpp::List plauscontourMC(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq,NumericVector truebz, NumericVector bzseq);
RcppExport SEXP imeiv_plauscontourMC(SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP bxseqSEXP, SEXP sampsizeSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type sampsize(sampsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bzseq(bzseqSEXP);
    __result = Rcpp::wrap(plauscontourMC(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq));
    return __result;
END_RCPP
}    
    

Rcpp::NumericMatrix sortmat(NumericMatrix x, unsigned int col);
RcppExport SEXP imeiv_sortmat(SEXP xSEXP, SEXP colSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type col(colSEXP);
    __result = Rcpp::wrap(sortmat(x,col));
    return __result;
END_RCPP    
}





