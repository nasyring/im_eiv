#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



Rcpp::List plauscontourGF(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq);
RcppExport SEXP imeiv_plauscontourGF(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP, SEXP bzseqSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bzseq(bzseqSEXP);
    __result = Rcpp::wrap(plauscontourGF(par,stat,del,type,n,propsd,truebx,truebz, bxseq, bzseq));
    return __result;
END_RCPP
}

Rcpp::List plauscontourGFu(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq);
RcppExport SEXP imeiv_plauscontourGFu(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP, SEXP bzseqSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bzseq(bzseqSEXP);
    __result = Rcpp::wrap(plauscontourGFu(par,stat,del,type,n,propsd,truebx,truebz, bxseq, bzseq));
    return __result;
END_RCPP
}

Rcpp::List plauscontourGFv(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector randsettype);
RcppExport SEXP imeiv_plauscontourGFv(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP, SEXP bxseqSEXP, SEXP randsettypeSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type randsettype(randsettypeSEXP);
    __result = Rcpp::wrap(plauscontourGFv(par,stat,del,type,n,propsd,truebx,bxseq, randsettype));
    return __result;
END_RCPP
}

Rcpp::List plauscontourGFa(NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq);
RcppExport SEXP imeiv_plauscontourGFa(SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP bxseqSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    __result = Rcpp::wrap(plauscontourGFa(stat,del,type,n,truebx,bxseq));
    return __result;
END_RCPP
}    
    
    
Rcpp::List plauscontourGF2(NumericVector par, NumericVector stat, NumericVector del, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq);
RcppExport SEXP imeiv_plauscontourGF2(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP, SEXP bzseqSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd(propsdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bzseq(bzseqSEXP);
    __result = Rcpp::wrap(plauscontourGF2(par,stat,del,n,propsd,truebx,truebz, bxseq, bzseq));
    return __result;
END_RCPP
}

Rcpp::List plauscontourIM(NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector sxseq, NumericVector seseq);
RcppExport SEXP imeiv_plauscontourIM(SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP, SEXP sxseqSEXP, SEXP seseqSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sxseq(sxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seseq(seseqSEXP);
    __result = Rcpp::wrap(plauscontourIM(stat,del,type,n,truebx,truebz,bxseq,sxseq,seseq));
    return __result;
END_RCPP
}

Rcpp::List plauscontourIMmarg(NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz, NumericVector bxseq);
RcppExport SEXP imeiv_plauscontourIMmarg(SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    __result = Rcpp::wrap(plauscontourIMmarg(stat,del,type,n,truebx,truebz,bxseq));
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





