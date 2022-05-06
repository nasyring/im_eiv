plauscontourGF <- function(par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq){
    .Call(`imeiv_plauscontourGF`, par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq)
}

plauscontourGFu <- function(par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq){
    .Call(`imeiv_plauscontourGFu`, par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq)
}

plauscontourGFv <- function(par,stat,del,type,n,propsd,truebx,bxseq,randsettype){
    .Call(`imeiv_plauscontourGFv`, par,stat,del,type,n,propsd,truebx,bxseq,randsettype)
}

plauscontourGFa <- function(stat,del,type,n,truebx,bxseq){
    .Call(`imeiv_plauscontourGFa`, stat,del,type,n,truebx,bxseq)
}

plauscontourGF2 <- function(par,stat,del,n,propsd,truebx,truebz,bxseq,bzseq){
    .Call(`imeiv_plauscontourGF2`, par,stat,del,n,propsd,truebx,truebz,bxseq,bzseq)
}

plauscontourIM <- function(stat,del,type,n,truebx,truebz,bxseq,sxseq,seseq){
    .Call(`imeiv_plauscontourIM`, stat,del,type,n,truebx,truebz,bxseq,sxseq,seseq)
}

plauscontourIMmarg <- function(stat,del,type,n,truebx,truebz,bxseq){
    .Call(`imeiv_plauscontourIMmarg`, stat,del,type,n,truebx,truebz,bxseq)
}

sortmat <- function(x,col){
    .Call(`imeiv_sortmat`, x,col)
}
