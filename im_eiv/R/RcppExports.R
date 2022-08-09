plausMC <- function(theta, intcpt, grid, stat, del, df, m_samps, intercept){
    .Call(`imeiv_plausMC`, theta, intcpt, grid, stat, del, df, m_samps, intercept)
}


plauscontourMCMC <- function(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype){
    .Call(`imeiv_plauscontourMCMC`, par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype)
}

plauscontourMCMCcond <- function(sampsize ,  stat ,   del ,  propsd , type , n , truebx , truebz , bxseq , sxseq , lenseq , plbxseq , plbzseq , lenplseq , se2){
    .Call(`imeiv_plauscontourMCMCcond`, sampsize ,  stat ,   del ,  propsd , type , n , truebx , truebz , bxseq , sxseq , lenseq , plbxseq , plbzseq , lenplseq , se2)
}
 
plauscontourSIR <- function(sampsize ,  stat ,   del  , n , mode , dens, se2, cond_par ){
    .Call(`imeiv_plauscontourSIR`, sampsize ,  stat ,   del  , n , mode , dens, se2, cond_par )
}

plauscontourMCMC2 <- function(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype){
    .Call(`imeiv_plauscontourMCMC2`, par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype)
}

plauscontourMC <- function(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred){
    .Call(`imeiv_plauscontourMC`, sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred)
}

plauscontourMC2 <- function(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred){
    .Call(`imeiv_plauscontourMC2`, sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq,randsetpred)
}

sortmat <- function(x,col){
    .Call(`imeiv_sortmat`, x,col)
}

grow <- function(x){
    .Call(`imeiv_grow`, x)
}

loglik <- function(theta, stat, del, n){
    .Call(`imeiv_loglik`, theta, stat, del, n)
}

optimrcpp <- function(theta, stat, del, n){
    .Call(`imeiv_optimrcpp`, theta, stat, del, n)
}


maxloglik <- function(thetas, stat, del, n){
    .Call(`imeiv_maxloglik`, thetas, stat, del, n)
}


genIMplaus <- function(thetas, stat, del, n, M){
    .Call(`imeiv_genIMplaus`, thetas, stat, del, n, M)
}
