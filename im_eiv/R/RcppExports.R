plauscontourMCMC <- function(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype){
    .Call(`imeiv_plauscontourMCMC`, par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype)
}

plauscontourMCMCcond <- function(sampsize ,  stat ,   del ,  propsd , type , n , truebx , truebz , bxseq , sxseq , lenseq , plbxseq , plbzseq , lenplseq , se2){
    .Call(`imeiv_plauscontourMCMCcond`, sampsize ,  stat ,   del ,  propsd , type , n , truebx , truebz , bxseq , sxseq , lenseq , plbxseq , plbzseq , lenplseq , se2)
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
