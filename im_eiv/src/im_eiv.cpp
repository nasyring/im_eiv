#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include <Rcpp.h>
#include <math.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>

Rcpp::List plauscontour(NumericVector par, NumericVector stat, NumericVector del, NumericVector n, NumericVector M, NumericVector propsd) {

	List result;
	NumericVector uu(1,0.0);
	NumericVector ct(5,0.0);
	NumericVector bx(1,0.0); bx[0] = par[0];
	NumericVector bz(1,0.0); bz[0] = par[1];	
	NumericVector mux(1,0.0); mux[0] = par[2];
	NumericVector sx(1,0.0); sx[0] = par[3];
	NumericVector se(1,0.0); sz[0] = par[4];

	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
	NumericVector L11(1,0.0); L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
	NumericVector L12(1,0.0); L12[0] = sx[0]*bx[0]/L11[0];
	NumericVector L22(1,0.0); 
	if((sx[0]/del[0]) > (L12[0]*L12[0])){
		L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
	}

	NumericVector v1(1, 0.0); v1[0] = s11[0]/L11[0];
	NumericVector v2(1, 0.0); v2[0] = (s12[0] - v1[0]*L12[0])/L22[0];
	NumericVector v3(1, 0.0); v3[0] = s22[0]/L22[0];
	NumericVector z1(1, 0.0); z1[0] = (ybar[0] - bz[0] - bx[0]*mux[0])/L11[0];
	NumericVector z2(1, 0.0); z2[0] = (wbar[0] - mux[0] - L12[0]*z1[0])/L22[0];
	
	NumericVector dL11dbx(1, 0.0); dL11dbx[0] = bx[0]*sx[0]/L11[0];	
	NumericVector dL11dsx(1, 0.0); dL11dsx[0] = 0.5*bx[0]*bx[0]/L11[0];	
	NumericVector dL11dse(1, 0.0); dL11dse[0] = 0.5/L11[0];	
	
	NumericVector dL12dbx(1, 0.0); dL12dbx[0] = (sx[0]*L11[0]-bx[0]*sx[0]*dL11dbx[0])/(L11[0]*L11[0]);	
	NumericVector dL12dsx(1, 0.0); dL12dsx[0] = (bx[0]*L11[0]-bx[0]*sx[0]*dL11dsx[0])/(L11[0]*L11[0]);
	NumericVector dL12dse(1, 0.0); dL12dse[0] = (-bx[0]*sx[0]*dL11dse[0])/(L11[0]*L11[0]);
	
	NumericVector dL22dbx(1, 0.0); dL22dbx[0] = -L12[0]*dL12dbx[0]/L22[0];
	NumericVector dL22dsx(1, 0.0); dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
	NumericVector dL22dse(1, 0.0); dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	
	NumericVector dv1dbx(1, 0.0); dv1dbx[0] = -s11[0]*dL11dbx[0]/(L11[0]*L11[0]);
	NumericVector dv1dbz(1, 0.0); 
	NumericVector dv1dmux(1, 0.0); 
	NumericVector dv1dsx(1, 0.0); dv1dsx[0] = -s11[0]*dL11dsx[0]/(L11[0]*L11[0]);
	NumericVector dv1dse(1, 0.0); dv1dse[0] = -s11[0]*dL11dse[0]/(L11[0]*L11[0]);
	
	NumericVector dv2dbx(1, 0.0); dv2dbx[0] = ((-L12[0]*dv1dbx[0]-v1[0]*dL12dbx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dbx[0])/(L22[0]*L22[0]);
	NumericVector dv2dbz(1, 0.0); 
	NumericVector dv2dmux(1, 0.0); 
	NumericVector dv2dsx(1, 0.0); dv2dsx[0] = ((-L12[0]*dv1dsx[0]-v1[0]*dL12dsx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dsx[0])/(L22[0]*L22[0]);
	NumericVector dv2dse(1, 0.0); dv2dse[0] = ((-L12[0]*dv1dse[0]-v1[0]*dL12dse[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dse[0])/(L22[0]*L22[0]);
	
	NumericVector dv3dbx(1, 0.0); dv3dbx[0] = -s22[0]*dL22dbx[0]/(L22[0]*L22[0]);
	NumericVector dv3dbz(1, 0.0); 
	NumericVector dv3dmux(1, 0.0); 
	NumericVector dv3dsx(1, 0.0); dv3dsx[0] = -s22[0]*dL22dsx[0]/(L22[0]*L22[0]);
	NumericVector dv3dse(1, 0.0); dv3dse[0] = -s22[0]*dL22dse[0]/(L22[0]*L22[0]);
	
	NumericVector dz1dbx(1, 0.0); dz1dbx[0] = (-mux[0]-z1[0]*dL11dbx[0])/L11[0];
	NumericVector dz1dbz(1, 0.0); dz1dbz[0] = -1/L11[0];
	NumericVector dz1dmux(1, 0.0); dz1dmux[0] = -bx[0]/L11[0];
	NumericVector dz1dsx(1, 0.0); dz1dsx[0] = -z1[0]*dL11dsx[0]/L11[0];
	NumericVector dz1dse(1, 0.0); dz1dse[0] = -z1[0]*dL11dse[0]/L11[0];
	
	NumericVector dz2dbx(1, 0.0); dz2dbx[0] = (-z1[0]*dL12dbx[0]-L12[0]*dz1dbx[0])/L22[0] - z2[0]*dL22dbx[0]/L22[0];
	NumericVector dz2dbz(1, 0.0); dz2dbz[0] = -L12[0]*dz1dbz[0]/L22[0];
	NumericVector dz2dmux(1, 0.0); dz2dmux[0] = (-1-L12[0]*dz1dmux[0])/L22[0];
	NumericVector dz2dsx(1, 0.0); dz2dsx[0] = (-z1[0]*dL12dsx[0]-L12[0]*dz1dsx[0])/L22[0] - z2[0]*dL22dsx[0]/L22[0];
	NumericVector dz2dse(1, 0.0); dz2dse[0] = (-z1[0]*dL12dse[0]-L12[0]*dz1dse[0])/L22[0] - z2[0]*dL22dse[0]/L22[0];
	
	mat J(5,5,fill::zeros);
	
	J = { { dv1dbx[0], dv1dbz[0], dv1dmux[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dbz[0], dv2dmux[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dbz[0], dv3dmux[0], dv3dsx[0], dv3dse[0] },
            { dz1dbx[0], dz1dbz[0], dz1dmux[0], dz1dsx[0], dz1dse[0] },
            { dz2dbx[0], dz2dbz[0], dz2dmux[0], dz2dsx[0], dz2dse[0] } };
	
	NumericVector det_J(1,0.0); det_J[0] = log(std::abs(arma::det(J)));
	
	NumericVector log_dens(1, 0.0);
	
	log_dens[0] = det_J[0] + (n[0]-2.0)*log(v1[0])-0.5*(v1[0]*v1[0]) + (n[0]-3.0)*log(v3[0])-0.5*(v3[0]*v3[0]) - 0.5*n[0]*(z1[0]*z1[0] + z2[0]*z2[0]) -  0.5*v2[0]*v2[0];



	/*  Begin MCMC  */
	
	NumericMatrix samples = NumericMatrix(10000, 5, zeroes.begin());
	NumericVector sampdens(10000,0.0);
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = log_dens[0];
	
	for(int j=0; j<100000; j++) {
		if(j>0){
			for(int i=0; i<4; i++){
				propsamp[i] = samples[j-1,i];
				currsamp[i] = samples[j-1,i];
			}
		}
		for(int i=0; i<4; i++){
			if(i>2){
				currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			}else {
				currsamp[i] = Rcpp::rgamma( 1, propsamp[i]/propsd[i], propsd[i] )
			}
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
			L12[0] = sx[0]*bx[0]/L11[0];
			if((sx[0]/del[0]) < (L12[0]*L12[0])){
				currdens[0] = -99.0;
			}else{ 	
				L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
				v1[0] = s11[0]/L11[0];
				v2[0] = (s12[0] - v1[0]*L12[0])/L22[0];
				v3[0] = s22[0]/L22[0];
				z1[0] = (ybar[0] - bz[0] - bx[0]*mux[0])/L11[0];
				z2[0] = (wbar[0] - mux[0] - L12[0]*z1[0])/L22[0];
	
				dL11dbx[0] = bx[0]*sx[0]/L11[0];	
				dL11dsx[0] = 0.5*bx[0]*bx[0]/L11[0];	
				dL11dse[0] = 0.5/L11[0];	
	
				dL12dbx[0] = (sx[0]*L11[0]-bx[0]*sx[0]*dL11dbx[0])/(L11[0]*L11[0]);	
				dL12dsx[0] = (bx[0]*L11[0]-bx[0]*sx[0]*dL11dsx[0])/(L11[0]*L11[0]);
				dL12dse[0] = (-bx[0]*sx[0]*dL11dse[0])/(L11[0]*L11[0]);
	
				dL22dbx[0] = -L12[0]*dL12dbx[0]/L22[0];
				dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
				dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	
				dv1dbx[0] = -s11[0]*dL11dbx[0]/(L11[0]*L11[0]);
				dv1dsx[0] = -s11[0]*dL11dsx[0]/(L11[0]*L11[0]);
				dv1dse[0] = -s11[0]*dL11dse[0]/(L11[0]*L11[0]);
	
				dv2dbx[0] = ((-L12[0]*dv1dbx[0]-v1[0]*dL12dbx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dbx[0])/(L22[0]*L22[0]);
				dv2dsx[0] = ((-L12[0]*dv1dsx[0]-v1[0]*dL12dsx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dsx[0])/(L22[0]*L22[0]);
				dv2dse[0] = ((-L12[0]*dv1dse[0]-v1[0]*dL12dse[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dse[0])/(L22[0]*L22[0]);
	
				dv3dbx[0] = -s22[0]*dL22dbx[0]/(L22[0]*L22[0]);
				dv3dsx[0] = -s22[0]*dL22dsx[0]/(L22[0]*L22[0]);
				dv3dse[0] = -s22[0]*dL22dse[0]/(L22[0]*L22[0]);
	
				dz1dbx[0] = (-mux[0]-z1[0]*dL11dbx[0])/L11[0];
				dz1dbz[0] = -1/L11[0];
				dz1dmux[0] = -bx[0]/L11[0];
				dz1dsx[0] = -z1[0]*dL11dsx[0]/L11[0];
				dz1dse[0] = -z1[0]*dL11dse[0]/L11[0];
	
				dz2dbx[0] = (-z1[0]*dL12dbx[0]-L12[0]*dz1dbx[0])/L22[0] - z2[0]*dL22dbx[0]/L22[0];
				dz2dbz[0] = -L12[0]*dz1dbz[0]/L22[0];
				dz2dmux[0] = (-1-L12[0]*dz1dmux[0])/L22[0];
				dz2dsx[0] = (-z1[0]*dL12dsx[0]-L12[0]*dz1dsx[0])/L22[0] - z2[0]*dL22dsx[0]/L22[0];
				dz2dse[0] = (-z1[0]*dL12dse[0]-L12[0]*dz1dse[0])/L22[0] - z2[0]*dL22dse[0]/L22[0];
	
				J = { { dv1dbx[0], dv1dbz[0], dv1dmux[0], dv1dsx[0], dv1dse[0] },
			            { dv2dbx[0], dv2dbz[0], dv2dmux[0], dv2dsx[0], dv2dse[0] },
  			 	         { dv3dbx[0], dv3dbz[0], dv3dmux[0], dv3dsx[0], dv3dse[0] },
  				          { dz1dbx[0], dz1dbz[0], dz1dmux[0], dz1dsx[0], dz1dse[0] },
     				       { dz2dbx[0], dz2dbz[0], dz2dmux[0], dz2dsx[0], dz2dse[0] } };
	
				det_J[0] = log(std::abs(arma::det(J)));

				currdens[0] = det_J[0] + (n[0]-2.0)*log(v1[0])-0.5*(v1[0]*v1[0]) + (n[0]-3.0)*log(v3[0])-0.5*(v3[0]*v3[0]) - 0.5*n[0]*(z1[0]*z1[0] + z2[0]*z2[0]) -  0.5*v2[0]*v2[0];
			}
			if(currdens[0] == -99.0){
				uu[0] = 1.0;
			}else {
				uu[0] = R::runif(0.0,1.0);
			}
			if(i>2){
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);	
			}else {
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0] + Rcpp::dgamma(currsamp[i], propsamp[i]/propsd[i], propsd[i], true )  - Rcpp::dgamma(propsamp[i], currsamp[i]/propsd[i], propsd[i], true )), 1.0);
			}	
			if(uu[0] < densdiff[0]){
				propsamp[i] = currsamp[i];
				propdens[0] = currdens[0];
				ct[i] = ct[i]+1.0;
			}
		}
		if( (j % 10) == 0 ){
			for(int i=0; i<4; i++){
				sample[j, i] = propsamp[i];
				sampdens[j] = propdens[0];
			}
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



Rcpp::List randsetsMCMC(NumericMatrix H, NumericMatrix A, NumericVector rL, NumericVector dimH, NumericVector M_samp) {
	
	List result;
	int M = int(M_samp[0]);
	int H2 = int(dimH[0]+2);
	int H1 = H2-2;
	NumericVector propsd(1,0.0);
	NumericVector u(2,1.0);
	NumericVector uprop(2,0.0);
	NumericVector logjointold(1,0.0);
	NumericVector logjointnew(1,0.0);
	NumericVector logjointdiff(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector zeroes = NumericVector(H2*1, 0.0); 
        NumericMatrix U = NumericMatrix(H2, 1, zeroes.begin());
	arma::mat Aa = as<arma::mat>(A);
	arma::mat arg;
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	arma::mat Ua = as<arma::mat>(U);
	
	for(int h = 0; h<H1; h++){
		U(h,0) = H(h,0);	
	}

for(int j=0; j<M; j++) {
		for(int i=0; i<2; i++){
			if( i==0 ){
				propsd[0] = 10.0;	
			}else {
				propsd[0] = 0.5;	
			}
			uprop(0) = u(0);uprop(1) = u(1);
			U(H1,0) = uprop[0];
			U(H1+1,0) = uprop[1];
			Ua = as<arma::mat>(U);
			arg = Aa*Ua;
			for(int h = 0; h<H2; h++){
				logjointold[0] = logjointold[0] + 0.5*rL[h]*arg(h,0)-0.5*exp(arg(h,0));	
			}
			uprop[i] = R::rnorm(u[i], propsd[0]);
			U(H1,0) = uprop[0];
			U(H1+1,0) = uprop[1];
			Ua = as<arma::mat>(U);
			arg = Aa*Ua;
			for(int h = 0; h<H2; h++){
				logjointnew[0] = logjointnew[0] + 0.5*rL[h]*arg(h,0)-0.5*exp(arg(h,0));	
			}
			logjointdiff[0] = logjointnew[0] - logjointold[0];
			logjointdiff[0] = fmin(std::exp(logjointdiff[0]), 1.0);
			uu[0] = R::runif(0.0,1.0);
			if(uu[0] <= logjointdiff[0]) {
				if(i==0){
					postsamples0[j] = uprop[0];	
				}else {
					postsamples1[j] = uprop[1];
				}
				u(0)=uprop(0);u(1)=uprop(1);
			}else {
				if(i==0){
					postsamples0[j] = u[0];	
				}else {
					postsamples1[j] = u[1];
				}				
			}
			logjointold[0] = 0.0; logjointnew[0] = 0.0;
		}
	}
result = Rcpp::List::create(Rcpp::Named("samples1") = postsamples0,Rcpp::Named("samples2") = postsamples1);

	return result;
	
	
}



Rcpp::List randsetspred(NumericMatrix S, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector U, NumericVector Ybar) {
	
	List result;
	int M = int(dimS[0]);
	int n = int(nsize[0]);
	int index = int(1);
	int dn_i = int(dimn_i[0]);
	NumericVector sumn_i2(1,0.0);
	NumericVector Qsampsw(10000,0.0);
	NumericVector Qsampsn(10000,0.0);
	NumericVector QsampsT(10000,0.0);
	NumericVector Qwl(1,0.0);
	NumericVector Qwu(1,0.0);
	NumericVector Qnl(1,0.0);
	NumericVector Qnu(1,0.0);
	NumericVector QTl(1,0.0);
	NumericVector QTu(1,0.0);
	NumericVector Ul(1,0.0);
	NumericVector Uu(1,0.0);
	NumericVector zeroes = NumericVector(10000*6, 0.0); 
        NumericMatrix randsetpred = NumericMatrix(10000, 6, zeroes.begin());
	NumericVector Z(1,0.0);
	
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
		
	for(int j=0; j < 10000; j++){
		Z[0] = R::rnorm(0.0,1.0);
		index = rand() % M;
		Qsampsw[j] = Z[0]*std::sqrt(S(index,0)*(1-(2*n_i[dn_i-1]/n)+(1/(n*n))*sumn_i2[0])+S(index,1)*((1/n)+1/(k[0])));
		Qsampsn[j] = Z[0]*std::sqrt(S(index,0)*(1+(1/(n*n))*sumn_i2[0])+S(index,1)*((1/n)+1/(k[0])));
		QsampsT[j] = Z[0]*std::sqrt(S(index,0)*(1+(1/(n*n))*sumn_i2[0])+S(index,1)*(1/n));
	}


	std::sort(Qsampsw.begin(), Qsampsw.end());
	std::sort(Qsampsn.begin(), Qsampsn.end());
	std::sort(QsampsT.begin(), QsampsT.end());
	
	for(int j=0; j < 10000; j++){
		Ul[0] = 0.5-std::fabs(U[j]-0.5);
		Uu[0] = 1.0-Ul[0];
		Qwl[0] = Qsampsw[std::round(10000*Ul[0])];
		Qwu[0] = Qsampsw[std::round(10000*Uu[0])];
		Qnl[0] = Qsampsn[std::round(10000*Ul[0])];
		Qnu[0] = Qsampsn[std::round(10000*Uu[0])];
		QTl[0] = QsampsT[std::round(10000*Ul[0])];
		QTu[0] = QsampsT[std::round(10000*Uu[0])];
		randsetpred(j,0) = Ybar[0]+Qwl[0];
		randsetpred(j,1) = Ybar[0]+Qwu[0];
		randsetpred(j,2) = Ybar[0]+Qnl[0];
		randsetpred(j,3) = Ybar[0]+Qnu[0];
		randsetpred(j,4) = Ybar[0]+QTl[0];
		randsetpred(j,5) = Ybar[0]+QTu[0];

	}
	

result = Rcpp::List::create(Rcpp::Named("randsetpred") = randsetpred);

	return result;
	
	
}
//NumericVector(*f)(NumericVector x, NumericVector uu, NumericVector vv, NumericVector yy, NumericVector zz)
Rcpp::NumericVector zeroin(NumericVector ax, NumericVector bx, NumericVector u, NumericVector v, NumericVector y, NumericVector z, Function f , NumericVector tol) {
    // code here
/*
double zeroin(ax,bx,f,tol)		An estimate to the root	
Xdouble ax;				Left border | of the range	
Xdouble bx;  				Right border| the root is seeked
Xdouble (*f)(double x);			Function under investigation	
Xdouble tol;				Acceptable tolerance		
	*/

NumericVector a(1,0.0); NumericVector b(1,0.0); NumericVector c(1,0.0); NumericVector fa(1,0.0); NumericVector fb(1,0.0); NumericVector fc(1,0.0);


  a[0] = ax[0];  b[0] = bx[0];  fa = f(a, u, v, y, z);  fb = f(b, u, v, y, z);
  c[0] = a[0];   fc[0] = fa[0];

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b[0]-a[0];		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
    double dbl_eps = std::numeric_limits<double>::epsilon();
	  
    if( fabs(fc[0]) < fabs(fb[0]) )
    {                         		/* Swap data for b to be the 	*/
	a[0] = b[0];  b[0] = c[0];  c[0] = a[0];          /* best approximation		*/
	fa[0]=fb[0];  fb[0]=fc[0];  fc[0]=fa[0];
    }
    tol_act = 2*dbl_eps*fabs(b[0]) + tol[0]/2;
    new_step = (c[0]-b[0])/2;

    if( fabs(new_step) <= tol_act || fb[0] == (double)0 )
    {
      return b;				/* Acceptable approx. is found	*/
    }

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	&& fabs(fa[0]) > fabs(fb[0]) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	register double t1,cb,t2;
	cb = c[0]-b[0];
	if( a[0]==c[0] )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb[0]/fa[0];			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa[0]/fc[0];  t1 = fb[0]/fc[0];  t2 = fb[0]/fa[0];
	  p = t2 * ( cb*q*(q-t1) - (b[0]-a[0])*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )
	{/* p was calculated with the op-*/
	  q = -q;	
	}/* posite sign; make p positive	*/
	else	
	{/* and assign possible minus to	*/
	  p = -p;
	}/* q				*/

	if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )
	{/* and isn't too large	*/
	  new_step = p/q;
	}/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act )
    {/* Adjust the step to be not less*/
      if( new_step > (double)0 )
      {/* than tolerance		*/
	new_step = tol_act;
      }
      else
      {
	new_step = -tol_act;
      }
    }

    a[0] = b[0];  fa[0] = fb[0];			/* Save the previous approx.	*/
    b[0] += new_step;  fb = f(b, u, v, y, z);	/* Do step to a new approxim.	*/
    if( (fb[0] > 0 && fc[0] > 0) || (fb[0] < 0 && fc[0] < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c[0] = a[0];  fc[0] = fa[0];                  /* opposite to that of b	*/
    }
  }

}

Rcpp::NumericVector root_function(NumericVector x, NumericVector Sampsj, NumericVector SL, NumericVector aL, NumericVector lambdaL) {
	
	int L = int(aL[0]);
	NumericVector sigsolnsj(1,0.0);
	NumericVector f(1,0.0);
	
	sigsolnsj[0] = std::exp(std::log(SL[L-1])-Sampsj[1]);	
	for(int k = 0; k < (L-1); k++){
		f[0] += (std::log(lambdaL[k]*x[0]+sigsolnsj[0])-std::log(SL[k]));
	}
	f[0] += Sampsj[0];
	
	return(f);
}



Rcpp::List sigmaSolve(NumericMatrix Samps, NumericVector SL, NumericVector aL, NumericVector aM, NumericVector lambdaL) {

	Rcpp::Function zeroin("zeroin");
	Rcpp::Function root_function("root_function");
	List result;
	int L = int(aL[0]);
	int M = int(aM[0]);
	NumericVector tol(1,0.0001);
	NumericVector sigsolnsj(1,0.0);
	NumericVector solnj(1,0.0);
	NumericVector soln(1,99.0);
	NumericVector zeroes = NumericVector(M*2, 0.0); 
        NumericMatrix solution = NumericMatrix(M, 2, zeroes.begin());
	NumericVector u(1,0.00001);
	NumericVector l(1,20000.0);
	NumericVector fu(1,0.0);
	NumericVector fl(1,0.0);
	NumericVector Sampsj(2,0.0);
	
	for(int k = 0; k < (M-1); k++){
		Sampsj[0] = Samps(k,0);Sampsj[1] = Samps(k,1);
		sigsolnsj[0] = std::exp(std::log(SL[L-1])-Sampsj[1]);	
		fl = root_function(l, Sampsj, SL, aL, lambdaL);
		fu = root_function(u, Sampsj, SL, aL, lambdaL);
		if(fl[0]*fu[0] < 0.0) 
		{
			soln = zeroin(l, u, Sampsj, SL, aL, lambdaL, root_function, tol);
		}
		solution(k,0) = soln[0]; solution(k,1) = sigsolnsj[0];
	}
	
	result = Rcpp::List::create(Rcpp::Named("solution") = solution);

	return result;
	
	
}
	
	


