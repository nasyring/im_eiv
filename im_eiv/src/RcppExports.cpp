#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;




Rcpp::List plausMC(NumericVector theta, NumericVector intcpt, NumericMatrix grid, NumericVector stat, NumericVector del, NumericVector df, int m_samps, bool intercept);
RcppExport SEXP imeiv_plausMC(SEXP thetaSEXP, SEXP intcptSEXP, SEXP gridSEXP, SEXP statSEXP, SEXP delSEXP, SEXP dfSEXP, SEXP m_sampsSEXP, SEXP interceptSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intcpt(intcptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type m_samps(m_sampsSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    __result = Rcpp::wrap(plausMC(theta, intcpt, grid, stat, del, df, m_samps, intercept));
    return __result;
END_RCPP
}    
    
Rcpp::List plausMCratio(NumericVector theta, NumericVector intcpt, NumericMatrix grid, NumericVector stat, NumericVector del, NumericVector df, int m_samps, bool intercept);
RcppExport SEXP imeiv_plausMCratio(SEXP thetaSEXP, SEXP intcptSEXP, SEXP gridSEXP, SEXP statSEXP, SEXP delSEXP, SEXP dfSEXP, SEXP m_sampsSEXP, SEXP interceptSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intcpt(intcptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type m_samps(m_sampsSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    __result = Rcpp::wrap(plausMCratio(theta, intcpt, grid, stat, del, df, m_samps, intercept));
    return __result;
END_RCPP
}   

Rcpp::List plausMCvar(NumericVector theta, NumericVector intcpt, NumericMatrix grid, NumericVector stat, NumericVector del, NumericVector df, int m_samps, bool intercept);
RcppExport SEXP imeiv_plausMCvar(SEXP thetaSEXP, SEXP intcptSEXP, SEXP gridSEXP, SEXP statSEXP, SEXP delSEXP, SEXP dfSEXP, SEXP m_sampsSEXP, SEXP interceptSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intcpt(intcptSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type df(dfSEXP);
    Rcpp::traits::input_parameter< int >::type m_samps(m_sampsSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    __result = Rcpp::wrap(plausMCvar(theta, intcpt, grid, stat, del, df, m_samps, intercept));
    return __result;
END_RCPP
}   


    


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
    __result = Rcpp::wrap(plauscontourMCMC(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype));
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

Rcpp::List plauscontourSIR(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector n, NumericVector mode, NumericVector dens, NumericVector se2, NumericVector cond_par);
RcppExport SEXP imeiv_plauscontourSIR(SEXP sampsizeSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP, SEXP modeSEXP, SEXP densSEXP, SEXP se2SEXP, SEXP cond_parSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type sampsize(sampsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dens(densSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se2(se2SEXP);  
    Rcpp::traits::input_parameter< NumericVector >::type cond_par(cond_parSEXP);
    __result = Rcpp::wrap(plauscontourSIR( sampsize ,  stat ,   del  , n , mode , dens, se2, cond_par ));
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


Rcpp::NumericVector loglik(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n);
RcppExport SEXP imeiv_loglik(SEXP thetaSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    __result = Rcpp::wrap(loglik(theta, stat, del, n));
    return __result;
END_RCPP    
}

Rcpp::NumericVector maxloglik(NumericMatrix thetas, NumericVector stat, NumericVector del, NumericVector n);
RcppExport SEXP imeiv_maxloglik(SEXP thetasSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    __result = Rcpp::wrap(maxloglik(thetas, stat, del, n));
    return __result;
END_RCPP    
}


Rcpp::List genIMplaus(NumericMatrix thetas, NumericVector stat, NumericVector del, NumericVector n, NumericVector M);
RcppExport SEXP imeiv_genIMplaus(SEXP thetasSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP, SEXP MSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M(MSEXP);
    __result = Rcpp::wrap(genIMplaus(thetas, stat, del, n, M));
    return __result;
END_RCPP    
}


Rcpp::List optimrcpp(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n);
RcppExport SEXP imeiv_optimrcpp(SEXP thetaSEXP, SEXP statSEXP, SEXP delSEXP, SEXP nSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stat(statSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type del(delSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    __result = Rcpp::wrap(optimrcpp(theta, stat, del, n));
    return __result;
END_RCPP    
}
