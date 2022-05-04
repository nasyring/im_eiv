plauscontourGF <- function(par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq){
    .Call(`imeiv_plauscontourGF`, par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq)
}

plauscontourGFu <- function(par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq){
    .Call(`imeiv_plauscontourGFu`, par,stat,del,type,n,propsd,truebx,truebz,bxseq,bzseq)
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
