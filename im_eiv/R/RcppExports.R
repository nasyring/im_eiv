plauscontourGF <- function(par,stat,del,n,propsd,truebx,truebz){
    .Call(`imeiv_plauscontourGF`, par,stat,del,n,propsd,truebx,truebz)
}

plauscontourIM <- function(stat,del,type,n,truebx,truebz,bxseq,sxseq,seseq){
    .Call(`imeiv_plauscontourIM`, stat,del,type,n,truebx,truebz,bxseq,sxseq,seseq)
}

sortmat <- function(x,col){
    .Call(`imeiv_sortmat`, x,col)
}
