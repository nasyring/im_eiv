plauscontourMCMC <- function(par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype){
    .Call(`imeiv_plauscontourMCMC`, par,stat,del,type,n,propsd,truebx,bxseq,truebz,bzseq,randsettype)
}

plauscontourMC <- function(sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq){
    .Call(`imeiv_plauscontourMC`, sampsize,stat,del,type,n,truebx,bxseq,truebz,bzseq)
}

sortmat <- function(x,col){
    .Call(`imeiv_sortmat`, x,col)
}
