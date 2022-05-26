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
    __result = Rcpp::wrap(plauscontourMCMC(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettypeSEXP));
    return __result;
END_RCPP
}

Rcpp::List plauscontourMCMCcond(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector propsd, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz,NumericVector bxseq, NumericVector sxseq, NumericVector lenseq, NumericVector plbxseq, NumericVector plbzseq, NumericVector lenplseq, NumericVector se2);
RcppExport SEXP imeiv_plauscontourMCMCcond(SEXP sampsizeSEXP, SEXP statSEXP, SEXP delSEXP, SEXP propsdSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP truebzSEXP, SEXP bxseqSEXP, SEXP sxseqSEXP, SEXP lenseqSEXP, SEXP plbxseqSEXP, SEXP plbzseqSEXP, SEXP lenplseqSEXP, SEXP se2SEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type sampsize(sampsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd(propsdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type type(typeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebx(truebxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type truebz(truebzSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bxseq(bxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sxseq(sxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lenseq(lenseqSEXP);  
    Rcpp::traits::input_parameter< NumericVector >::type plbxseq(plbxseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type plbzseq(plbzseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lenplseq(lenplseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se2(se2SEXP);  
    __result = Rcpp::wrap(plauscontourMCMCcond( sampsize ,  stat ,   del ,  propsd , type , n , truebx , truebz , bxseq , sxseq , lenseq , plbxseq , plbzseq , lenplseq , se2 ));
    return __result;
END_RCPP
}

Rcpp::List plauscontourSIR(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector n, NumericVector mode, NumericVector local_pt, NumericVector se2, NumericVector cond_par) {


Rcpp::List plauscontourSIR(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector n, NumericVector mode, NumericVector local_pt, NumericVector se2, NumericVector cond_par);
RcppExport SEXP imeiv_plauscontourSIR(SEXP sampsizeSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP, SEXP modeSEXP, SEXP local_ptSEXP, SEXP se2SEXP, SEXP cond_parSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type sampsize(sampsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type local_pt(local_ptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se2(se2SEXP);  
    Rcpp::traits::input_parameter< NumericVector >::type cond_par(cond_parSEXP);
    __result = Rcpp::wrap(plauscontourSIR( sampsize ,  stat ,   del  , n , mode , local_pt, se2, cond_par ));
    return __result;
END_RCPP
}

Rcpp::List plauscontourMCMC2(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector truebz, NumericVector bzseq, NumericVector randsettype);
RcppExport SEXP imeiv_plauscontourMCMC2(SEXP parSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP propsdSEXP, SEXP truebxSEXP , SEXP bxseqSEXP, SEXP truebzSEXP, SEXP bzseqSEXP, SEXP randsettypeSEXP){
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
    __result = Rcpp::wrap(plauscontourMCMC2(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettypeSEXP));
    return __result;
END_RCPP
}

Rcpp::List plauscontourMC(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq,NumericVector truebz, NumericVector bzseq, NumericVector randsetpred);
RcppExport SEXP imeiv_plauscontourMC(SEXP sampsizeSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP bxseqSEXP, SEXP truebzSEXP, SEXP bzseqSEXP, SEXP randsetpredSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type randsetpred(randsetpredSEXP);
    __result = Rcpp::wrap(plauscontourMC(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred));
    return __result;
END_RCPP
}    
    

Rcpp::List plauscontourMC2(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq,NumericVector truebz, NumericVector bzseq, NumericVector randsetpred);
RcppExport SEXP imeiv_plauscontourMC2(SEXP sampsizeSEXP, SEXP statSEXP, SEXP delSEXP, SEXP typeSEXP, SEXP nSEXP, SEXP truebxSEXP, SEXP bxseqSEXP, SEXP truebzSEXP, SEXP bzseqSEXP, SEXP randsetpredSEXP){
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
    Rcpp::traits::input_parameter< NumericVector >::type randsetpred(randsetpredSEXP);
    __result = Rcpp::wrap(plauscontourMC2(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred));
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

Rcpp::NumericVector grow(NumericVector x);
RcppExport SEXP imeiv_grow(SEXP xSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    __result = Rcpp::wrap(grow(x));
    return __result;
END_RCPP    
}




