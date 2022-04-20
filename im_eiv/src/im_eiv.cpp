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

Rcpp::List plauscontour(NumericVector par, NumericVector stat, NumericVector del, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz) {

	List result;
	NumericVector tempdens(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector ct(5,0.0);
	NumericVector bx(1,0.0); bx[0] = par[0];
	NumericVector bz(1,0.0); bz[0] = par[1];	
	NumericVector mux(1,0.0); mux[0] = par[2];
	NumericVector sx(1,0.0); sx[0] = par[3];
	NumericVector se(1,0.0); se[0] = par[4];
	NumericVector sd(1,0.0);  sd[0] = std::sqrt(1.0/n[0]);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
	NumericVector L11(1,0.0); L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
	NumericVector L12(1,0.0); L12[0] = sx[0]*bx[0]/L11[0];
	NumericVector L22(1,1.0); 
	if((sx[0]/del[0]) > (L12[0]*L12[0])){
		L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
	}

	NumericVector v1(1, 0.0); v1[0] = s11[0]/L11[0];
	NumericVector v2(1, 0.0); v2[0] = (s12[0] - v1[0]*L12[0])/L22[0];
	NumericVector v3(1, 0.0); v3[0] = s22[0]/L22[0];
	NumericVector z1(1, 0.0); z1[0] = (ybar[0] - bz[0] - bx[0]*mux[0])/L11[0];
	NumericVector z2(1, 0.0); z2[0] = (wbar[0] - mux[0] - L12[0]*z1[0])/L22[0];
	
	if(((sx[0]/del[0]) > (L12[0]*L12[0])) & (v1[0] > 0) & (v3[0] > 0)){
	
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
	
	NumericVector detJ(1,0.0); 
	
	detJ[0] = log(std::abs(arma::det(J)));
	
	NumericVector logdens(1, 0.0);
	
	//logdens[0] = detJ[0] + (n[0]-2.0)*log(v1[0])-0.5*(v1[0]*v1[0]) + (n[0]-3.0)*log(v3[0])-0.5*(v3[0]*v3[0]) - 0.5*n[0]*(z1[0]*z1[0] + z2[0]*z2[0]) -  0.5*v2[0]*v2[0];
	logdens[0] = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
	tempdens[0] = logdens[0]; 
	//  Begin MCMC  
	
	NumericVector zeroes(50000,0.0);
	NumericMatrix samples = NumericMatrix(10000, 5, zeroes.begin());
	samples(0,0) = bx[0];samples(0,1) = bz[0];samples(0,2) = mux[0];samples(0,3) = sx[0];samples(0,4) = se[0];
	NumericVector sampdens(10000,0.0);
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = logdens[0];
	
	
	for(int j=0; j<100000; j++) {
		for(int i=0; i<5; i++){
			for(int k=0; k<5; k++){
				currsamp[k] = propsamp[k];
			}
			if(i<3){
				currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			}else {
				currsamp[i] = R::rgamma(propsd[i] + propsamp[i]/propsd[i], propsd[i] );
			}
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
			L12[0] = sx[0]*bx[0]/L11[0];
			if(  ((sx[0]/del[0]) > (L12[0]*L12[0]))  ){
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
	
				detJ[0] = log(std::abs(arma::det(J)));
				if((v1[0]>0) & (v3[0]>0)){
					currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
				}else {
					currdens[0]  = -1000000;
				}
			}else {
				currdens[0]  = -1000000;
			}
			uu[0] = R::runif(0.0,1.0);
			if(i<3){
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);	
			}else {
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0] + R::dgamma(currsamp[i], propsd[i]+propsamp[i]/propsd[i], propsd[i], true )  - R::dgamma(propsamp[i], propsd[i]+currsamp[i]/propsd[i], propsd[i], true )), 1.0);
			}	
			if(uu[0] < densdiff[0]){
				propsamp[i] = currsamp[i];
				propdens[0] = currdens[0];
				ct[i] = ct[i]+1.0;
			}
		}
		if( (j % 10) == 0 ){
			for(int i=0; i<5; i++){
				samples(j/10, i) = propsamp[i];
				sampdens[j/10] = propdens[0];
			}
		}
	}
	for(int i=0; i<5; i++){
		ct[i]=ct[i]/100000.0;
	}
	
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("samples") = samples, Rcpp::Named("logdens") = sampdens, Rcpp::Named("tempdens") = tempdens);
	}	
	return result;
	
	
	/*
	
	NumericVector unifs(1,0.0);NumericVector unifs_hi(10000,0.0);NumericVector unifs_lo(10000,0.0);
	NumericVector bxs(10000,0.0);NumericVector bzs(10000,0.0);
	for(int i=0; i<9999; i++){
		unifs[0] = R::runif(0.0,1.0);
		unifs_hi[i] = 0.5 + fabs(unifs[0] - 0.5); 
		unifs_lo[i] = 1.0-unifs_hi[i];
		bxs[i] = samples(i, 0);
		bzs[i] = samples(i, 1);
	}
	std::sort(bxs.begin(), bxs.end());
	std::sort(bzs.begin(), bzs.end());
	NumericVector bxseq(501,0.0);NumericVector bzseq(501,0.0);
	for(int i=0; i<499; i++){
		bxseq[i] = bxs[i*20];bzseq[i] = bzs[i*20];
	}
	bxseq[500] = bx[9999]; bxseq[500] = bx[9999];
	NumericVector plausesx(501,0.0);NumericVector plausesz(501,0.0);
	for(int i=0; i<500; i++){
		for(int j=0; j<9999; j++){
			if(   (bxseq[i] < bxs[floor(9999*unifs_lo[j])]) & (bxseq[i] > bxs[ceil(9999*unifs_hi[j])])   ){
				plausesx[i] = plausesx[i]+0.0001;
			}			
			if(   (bzseq[i] < bzs[floor(9999*unifs_lo[j])]) & (bzseq[i] > bzs[ceil(9999*unifs_hi[j])])   ){
				plausesz[i] = plausesz[i]+0.0001;
			}
		}
	}
	NumericVector plausestrux(1,0.0);NumericVector plausestruz(1,0.0);
	for(int j=0; j<9999; j++){
		if(   (truebx[0] < bxs[floor(9999*unifs_lo[j])]) & (truebx[0] > bxs[ceil(9999*unifs_hi[j])])   ){
			plausestrux[0] = plausestrux[0]+0.0001;
		}			
		if(   (truebz[0] < bzs[floor(9999*unifs_lo[j])]) & (truebz[0] > bzs[ceil(9999*unifs_hi[j])])   ){
			plausestruz[0] = plausestruz[0]+0.0001;
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("beta_x_seq") = bxseq, Rcpp::Named("beta_z_seq") = bzseq);

	return result;
	*/
}
	
	
	
	
