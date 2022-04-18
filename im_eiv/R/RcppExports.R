
randsetsMCMC <- function(H,A,rL,dimH,M_samp) {
    .Call(`impred_randsetsMCMC`, H,A,rL,dimH,M_samp)
}

randsetspred <- function(S,dimS,nsize,n_i,dimn_i,k,U,Ybar) {
    .Call(`impred_randsetspred`, S,dimS,nsize,n_i,dimn_i,k,U,Ybar)
}


sigmaSolve <- function(Sampsj, SL, aL, aM, lambdaL) {
    .Call(`impred_sigmaSolve`, Sampsj, SL, aL, aM, lambdaL)
}

zeroin <- function(ax, bx, u, v, y, z, f, tol) {
    .Call(`impred_zeroin`, ax, bx, u, v, y, z, f, tol)
}

root_function <- function(x, Sampsj, SL, aL, lambdaL) {
    .Call(`impred_root_function`, x, Sampsj, SL, aL, lambdaL)
}
