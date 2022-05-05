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

Rcpp::List plauscontourGF(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq) {

	List result;
	Rcpp::Function sortmat("sortmat");
	
	NumericVector check(1, 0.0);
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
	
	if(type[0] == 1.0){
		check[0] = (sx[0]/del[0]) - (L12[0]*L12[0]);
	}else if(type[0] == 2.0){
		check[0] = (sx[0] + se[0]/del[0]) - (L12[0]*L12[0]);	
	}else if(type[0] == 3.0){
		check[0] = (sx[0] + del[0]) - (L12[0]*L12[0]);	
	}else {
		check[0] = -1;	
	}	
	
	if(check[0] > 0){
		
	if(type[0] == 1.0){
		L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
	}else if(type[0] == 2.0){
		L22[0] = std::sqrt((sx[0] + se[0]/del[0]) - (L12[0]*L12[0]));
	}else {
		L22[0] = std::sqrt((sx[0] + del[0]) - (L12[0]*L12[0]));
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
	NumericVector dL22dsx(1, 0.0); 
	NumericVector dL22dse(1, 0.0); 
		
	if(type[0] == 1.0){
		dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	}else if(type[0] == 2.0){
		dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = (0.5/L22[0])*(1.0/del[0] - 2.0*L12[0]*dL12dse[0]);
	}else {
		dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	}	
	
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
	
	
	mat J(5,5,fill::zeros); mat J2(3,3,fill::zeros);
	
	J = { { dv1dbx[0], dv1dbz[0], dv1dmux[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dbz[0], dv2dmux[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dbz[0], dv3dmux[0], dv3dsx[0], dv3dse[0] },
            { dz1dbx[0], dz1dbz[0], dz1dmux[0], dz1dsx[0], dz1dse[0] },
            { dz2dbx[0], dz2dbz[0], dz2dmux[0], dz2dsx[0], dz2dse[0] } };
		
	J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
	NumericVector detJ(1,0.0); NumericVector detJ2(1,0.0); 
	
	detJ[0] = log(std::abs(arma::det(J)));
	detJ2[0] = log(std::abs(arma::det(J2)));
	
	NumericVector logdens(1, 0.0); NumericVector logdens2(1, 0.0);
	
	logdens[0] = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
	logdens2[0] = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
	
	tempdens[0] = logdens[0]; 
	//  Begin MCMC  
	
	NumericVector zeroes(70000,0.0);
	NumericMatrix samples = NumericMatrix(10000, 7, zeroes.begin());
	NumericMatrix samples2 = NumericMatrix(10000, 7, zeroes.begin());
	samples(0,0) = bx[0];samples(0,1) = bz[0];samples(0,2) = mux[0];samples(0,3) = sx[0];samples(0,4) = se[0];samples(0,5) = logdens[0];samples(0,6) = logdens2[0];
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = logdens[0];
	NumericVector currdens2(1,0.0);NumericVector propdens2(1,0.0); propdens2[0] = logdens2[0];
	
	
	for(int j=0; j<200000; j++) {
		for(int i=0; i<5; i++){
			for(int k=0; k<5; k++){
				currsamp[k] = propsamp[k];
			}
			currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			if(!((i > 2) & (currsamp[i]<0))){
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
			L12[0] = sx[0]*bx[0]/L11[0];
			if(type[0] == 1.0){
				check[0] = (sx[0]/del[0]) - (L12[0]*L12[0]);
			}else if(type[0] == 2.0){
				check[0] = (sx[0] + se[0]/del[0]) - (L12[0]*L12[0]);	
			}else if(type[0] == 3.0){
				check[0] = (sx[0] + del[0]) - (L12[0]*L12[0]);	
			}else {
				check[0] = -1;	
			}	
			if(check[0] > 0){	
				if(type[0] == 1.0){
					L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
				}else if(type[0] == 2.0){
					L22[0] = std::sqrt((sx[0] + se[0]/del[0]) - (L12[0]*L12[0]));
				}else {
					L22[0] = std::sqrt((sx[0] + del[0]) - (L12[0]*L12[0]));
				}
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
				if(type[0] == 1.0){
					dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
				}else if(type[0] == 2.0){
					dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = (0.5/L22[0])*(1.0/del[0] - 2.0*L12[0]*dL12dse[0]);
				}else {
					dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
				}
	
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
				
				J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
			            { dv2dbx[0],dv2dsx[0], dv2dse[0] },
  			 	         { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
				detJ[0] = log(std::abs(arma::det(J))); detJ2[0] = log(std::abs(arma::det(J2)));
				if((v1[0]>0) & (v3[0]>0)){
					currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
					currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				}else {
					currdens[0]  = -1000000;
					currdens2[0]  = -1000000;
				}
			}else {
				currdens[0]  = -1000000;
				currdens2[0]  = -1000000;
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
				propdens2[0] = currdens2[0];
				ct[i] = ct[i]+1.0;
			}
			}
		}
		if( ((j % 10) == 0) & (j < 100000) ){
			for(int i=0; i<5; i++){
				samples(j/10, i) = propsamp[i];
			}
			samples(j/10, 5) = propdens[0];
			samples(j/10, 6) = propdens2[0];
		}
		if( ((j % 10) == 0) & (j >= 100000) ){
			for(int i=0; i<5; i++){
				samples2((j-100000)/10, i) = propsamp[i];

			}
			samples2((j-100000)/10, 5) = propdens[0];
			samples2((j-100000)/10, 6) = propdens2[0];
		}
	}
	for(int i=0; i<5; i++){
		ct[i]=ct[i]/100000.0;
	}
		
		
	// plausibility of betax
		
	NumericVector plausestrux(1, 0.0);
	samples = sortmat(samples,6);
	NumericVector plausbxseq(500,0.0);
	NumericVector randsetlo(10000,0.0);NumericVector randsethi(10000,0.0);
	int ind = 0;	
		
	for(int i = 0; i < 10000; i++){
		ind = 0;
		if(samples2(i,6) > samples(0,6)){
			if(samples2(i,6) < samples(9999,6)){
				while(samples2(i,6) > samples(ind,6)){
					ind = ind+1;
				}
				NumericVector subset(10000-ind,0.0);
				for(int j = ind; j < 10000; j++){
					subset[j-ind] = samples(j,0);
				}
				randsetlo[i] = 	Rcpp::min(subset);
				randsethi[i] =  Rcpp::max(subset);
			}else {
				randsetlo[i] = 	0;
				randsethi[i] =  0;
			}
		}else {
			NumericVector subset(10000,0.0);
			for(int j = 0; j < 10000; j++){
				subset[j] = samples(j,0);
			}
			randsetlo[i] = 	Rcpp::min(subset);
			randsethi[i] =  Rcpp::max(subset);
		}
		for(int j = 0; j < 500; j++){
			if((bxseq[j] > randsetlo[i]) & (bxseq[j] < randsethi[i])){
				plausbxseq[j] = plausbxseq[j] + 0.0001;	
			}
		}
		if((truebx[0] > randsetlo[i]) & (truebx[0] < randsethi[i])){
			plausestrux[0] = plausestrux[0] + 0.0001;	
		}
	}
		
	NumericVector plausestruz(1, 0.0);
	samples = sortmat(samples,5);
	NumericVector plausbzseq(500,0.0);
	NumericVector randsetloz(10000,0.0);NumericVector randsethiz(10000,0.0);
	ind = 0;	
		
	for(int i = 0; i < 10000; i++){
		ind = 0;
		if(samples2(i,5) > samples(0,5)){
			if(samples2(i,5) < samples(9999,5)){
				while(samples2(i,5) > samples(ind,5)){
					ind = ind+1;
				}
				NumericVector subset(10000-ind,0.0);
				for(int j = ind; j < 10000; j++){
					subset[j-ind] = samples(j,1);
				}
				randsetloz[i] = Rcpp::min(subset);
				randsethiz[i] =  Rcpp::max(subset);
			}else {
				randsetloz[i] = 0;
				randsethiz[i] =  0;
			}
		}else {
			NumericVector subset(10000,0.0);
			for(int j = 0; j < 10000; j++){
				subset[j] = samples(j,1);
			}
			randsetloz[i] = Rcpp::min(subset);
			randsethiz[i] =  Rcpp::max(subset);
		}
		for(int j = 0; j < 500; j++){
			if((bzseq[j] > randsetloz[i]) & (bzseq[j] < randsethiz[i])){
				plausbzseq[j] = plausbzseq[j] + 0.0001;	
			}
		}
		if((truebz[0] > randsetloz[i]) & (truebz[0] < randsethiz[i])){
			plausestruz[0] = plausestruz[0] + 0.0001;	
		}
	}
		
		
		
		
	/*
	
	NumericVector unifs(1,0.0);NumericVector unifs_hi(10000,0.0);NumericVector unifs_lo(10000,0.0);
	NumericVector bxs(10000,0.0);NumericVector bzs(10000,0.0);
	for(int i=0; i<10000; i++){
		unifs[0] = R::runif(0.0,1.0);
		unifs_hi[i] = 0.5 + fabs(unifs[0] - 0.5); 
		unifs_lo[i] = 1.0-unifs_hi[i];
		bxs[i] = samples(i, 0);
		bzs[i] = samples(i, 1);
	}
	std::sort(bxs.begin(), bxs.end());
	std::sort(bzs.begin(), bzs.end());
	NumericVector bxseq(501,0.0);NumericVector bzseq(501,0.0);
	for(int i=0; i<500; i++){
		bxseq[i] = bxs[i*20];bzseq[i] = bzs[i*20];
	}
	bxseq[500] = bxs[9999]; bzseq[500] = bzs[9999];
	NumericVector plausesx(501,0.0);NumericVector plausesz(501,0.0);
	for(int i=0; i<500; i++){
		for(int j=0; j<10000; j++){
			if(   (bxseq[i] > bxs[floor(9999*unifs_lo[j])]) & (bxseq[i] < bxs[ceil(9999*unifs_hi[j])])   ){
				plausesx[i] = plausesx[i]+0.0001;
			}			
			if(   (bzseq[i] > bzs[floor(9999*unifs_lo[j])]) & (bzseq[i] < bzs[ceil(9999*unifs_hi[j])])   ){
				plausesz[i] = plausesz[i]+0.0001;
			}
		}
	}
	NumericVector plausestrux(1,0.0);NumericVector plausestruz(1,0.0);
	for(int j=0; j<10000; j++){
		if(   (truebx[0] > bxs[floor(9999*unifs_lo[j])]) & (truebx[0] < bxs[ceil(9999*unifs_hi[j])])   ){
			plausestrux[0] = plausestrux[0]+0.0001;
		}			
		if(   (truebz[0] > bzs[floor(9999*unifs_lo[j])]) & (truebz[0] < bzs[ceil(9999*unifs_hi[j])])   ){
			plausestruz[0] = plausestruz[0]+0.0001;
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("beta_x_seq") = bxseq, Rcpp::Named("beta_z_seq") = bzseq);
	*/
		
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_x") = plausbxseq, Rcpp::Named("plauses_beta_z") = plausbzseq, Rcpp::Named("randsetlo") = randsetlo, Rcpp::Named("randsethi") = randsethi, Rcpp::Named("samples") = samples, Rcpp::Named("samples2") = samples2);		
	}
	return result;
}



Rcpp::List plauscontourGFu(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq) {

	List result;
	Rcpp::Function sortmat("sortmat");
	
	NumericVector check(1, 0.0);
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
	
	if(type[0] == 1.0){
		check[0] = (sx[0]/del[0]) - (L12[0]*L12[0]);
	}else if(type[0] == 2.0){
		check[0] = (sx[0] + se[0]/del[0]) - (L12[0]*L12[0]);	
	}else if(type[0] == 3.0){
		check[0] = (sx[0] + del[0]) - (L12[0]*L12[0]);	
	}else {
		check[0] = -1;	
	}	
	
	if(check[0] > 0){
		
	if(type[0] == 1.0){
		L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
	}else if(type[0] == 2.0){
		L22[0] = std::sqrt((sx[0] + se[0]/del[0]) - (L12[0]*L12[0]));
	}else {
		L22[0] = std::sqrt((sx[0] + del[0]) - (L12[0]*L12[0]));
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
	NumericVector dL22dsx(1, 0.0); 
	NumericVector dL22dse(1, 0.0); 
		
	if(type[0] == 1.0){
		dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	}else if(type[0] == 2.0){
		dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = (0.5/L22[0])*(1.0/del[0] - 2.0*L12[0]*dL12dse[0]);
	}else {
		dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
		dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
	}	
	
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
	
	
	mat J(5,5,fill::zeros); mat J2(3,3,fill::zeros);
	
	J = { { dv1dbx[0], dv1dbz[0], dv1dmux[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dbz[0], dv2dmux[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dbz[0], dv3dmux[0], dv3dsx[0], dv3dse[0] },
            { dz1dbx[0], dz1dbz[0], dz1dmux[0], dz1dsx[0], dz1dse[0] },
            { dz2dbx[0], dz2dbz[0], dz2dmux[0], dz2dsx[0], dz2dse[0] } };
		
	J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
	NumericVector detJ(1,0.0); NumericVector detJ2(1,0.0); 
	
	detJ[0] = log(std::abs(arma::det(J)));
	detJ2[0] = log(std::abs(arma::det(J2)));
	
	NumericVector logdens(1, 0.0); NumericVector logdens2(1, 0.0);
	
	logdens[0] = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
	logdens2[0] = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
	
	tempdens[0] = logdens[0]; 
	//  Begin MCMC  
	
	NumericVector zeroes(140000,0.0);NumericVector zeroes2(70000,0.0);
	NumericMatrix samples = NumericMatrix(20000, 7, zeroes.begin());
	samples(0,0) = bx[0];samples(0,1) = bz[0];samples(0,2) = mux[0];samples(0,3) = sx[0];samples(0,4) = se[0];samples(0,5) = logdens[0];samples(0,6) = logdens2[0];
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = logdens[0];
	NumericVector currdens2(1,0.0);NumericVector propdens2(1,0.0); propdens2[0] = logdens2[0];
	
	
	for(int j=0; j<200000; j++) {
		for(int i=0; i<5; i++){
			for(int k=0; k<5; k++){
				currsamp[k] = propsamp[k];
			}
			currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			if(!((i > 2) & (currsamp[i]<0.0))){
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			L11[0] = std::sqrt(se[0]+sx[0]*bx[0]*bx[0]);
			L12[0] = sx[0]*bx[0]/L11[0];
			if(type[0] == 1.0){
				check[0] = (sx[0]/del[0]) - (L12[0]*L12[0]);
			}else if(type[0] == 2.0){
				check[0] = (sx[0] + se[0]/del[0]) - (L12[0]*L12[0]);	
			}else if(type[0] == 3.0){
				check[0] = (sx[0] + del[0]) - (L12[0]*L12[0]);	
			}else {
				check[0] = -1;	
			}	
			if(check[0] > 0){	
				if(type[0] == 1.0){
					L22[0] = std::sqrt(sx[0]/del[0] - L12[0]*L12[0]);
				}else if(type[0] == 2.0){
					L22[0] = std::sqrt((sx[0] + se[0]/del[0]) - (L12[0]*L12[0]));
				}else {
					L22[0] = std::sqrt((sx[0] + del[0]) - (L12[0]*L12[0]));
				}
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
				if(type[0] == 1.0){
					dL22dsx[0] = (0.5/L22[0])*((1/del[0]) - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
				}else if(type[0] == 2.0){
					dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = (0.5/L22[0])*(1.0/del[0] - 2.0*L12[0]*dL12dse[0]);
				}else {
					dL22dsx[0] = (0.5/L22[0])*(1.0 - 2.0*L12[0]*dL12dsx[0]);
					dL22dse[0] = -L12[0]*dL12dse[0]/L22[0];
				}
	
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
				
				J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
			            { dv2dbx[0],dv2dsx[0], dv2dse[0] },
  			 	         { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
				detJ[0] = log(std::abs(arma::det(J))); detJ2[0] = log(std::abs(arma::det(J2)));
				if((v1[0]>0) & (v3[0]>0)){
					currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
					currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				}
				uu[0] = R::runif(0.0,1.0);
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);
			} else {
				uu[0] = 2.0;
				densdiff[0] = 0;
			}		
			if(uu[0] < densdiff[0]){
				propsamp[i] = currsamp[i];
				propdens[0] = currdens[0];
				propdens2[0] = currdens2[0];
				ct[i] = ct[i]+1.0;
			}
			}
		}
		if( ((j % 10) == 0) ){
			for(int i=0; i<5; i++){
				samples(j/10, i) = propsamp[i];
			}
			samples(j/10, 5) = propdens[0];
			samples(j/10, 6) = propdens2[0];
		}
	}
	for(int i=0; i<5; i++){
		ct[i]=ct[i]/200000.0;
	}
	
	//result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("samples") = samples);		
	//}
	//return result;
		
		
	// plausibility

	NumericVector unifs(3,0.0);NumericVector maxunifs(1,0.0);NumericVector unifs_hi(20000,0.0);NumericVector unifs_lo(20000,0.0);
	NumericVector unifsz(3,0.0);NumericVector maxunifsz(1,0.0);NumericVector unifs_hiz(20000,0.0);NumericVector unifs_loz(20000,0.0);
	NumericVector bxs(20000,0.0);NumericVector bzs(20000,0.0);
	for(int i=0; i<20000; i++){
		unifs[0] = R::runif(0.0,1.0); unifs[1] = R::runif(0.0,1.0); unifs[2] = R::runif(0.0,1.0);
		maxunifs[0] = fmax(unifs[0], unifs[1]); maxunifs[0] = fmax(maxunifs[0], unifs[2]); 
		unifs_hi[i] = 0.5 + fabs(maxunifs[0] - 0.5); 
		unifs_lo[i] = 1.0-unifs_hi[i];
		unifsz[0] = R::runif(0.0,1.0); unifsz[1] = R::runif(0.0,1.0);unifsz[2] = R::runif(0.0,1.0);unifsz[3] = R::runif(0.0,1.0);unifsz[4] = R::runif(0.0,1.0);
		maxunifsz[0] = fmax(unifsz[0], unifsz[1]); maxunifsz[0] = fmax(maxunifsz[0], unifsz[2]); maxunifsz[0] = fmax(maxunifsz[0], unifsz[3]);maxunifsz[0] = fmax(maxunifsz[0], unifsz[4]);  
		unifs_hiz[i] = 0.5 + fabs(maxunifsz[0] - 0.5); 
		unifs_loz[i] = 1.0-unifs_hiz[i];	
		bxs[i] = samples(i, 0);
		bzs[i] = samples(i, 1);
	}
	std::sort(bxs.begin(), bxs.end());
	std::sort(bzs.begin(), bzs.end());
	NumericVector plausesx(500,0.0);NumericVector plausesz(500,0.0);
	int intlo = 0;  int inthi = 0;
	for(int i=0; i<500; i++){
		for(int j=0; j<20000; j++){
			intlo = round(floor(19999*unifs_lo[j]));
			if(intlo < 0){
				intlo = 0;	
			}
			inthi = round(ceil(19999*unifs_hi[j]));
			if(inthi > 19999){
				inthi = 19999;	
			}
			if(   (bxseq[i] > bxs[intlo]) & (bxseq[i] < bxs[inthi])   ){
				plausesx[i] = plausesx[i]+0.00005;
			}	
			intlo = round(floor(19999*unifs_loz[j]));
			if(intlo < 0){
				intlo = 0;	
			}
			inthi = round(ceil(19999*unifs_hiz[j]));
			if(inthi > 19999){
				inthi = 19999;	
			}
			if(   (bzseq[i] > bzs[intlo]) & (bzseq[i] < bzs[inthi])   ){
				plausesz[i] = plausesz[i]+0.00005;
			}
		}
	}
	NumericVector plausestrux(1,0.0);NumericVector plausestruz(1,0.0);
	for(int j=0; j<20000; j++){
			intlo = round(floor(19999*unifs_lo[j]));
			if(intlo < 0){
				intlo = 0;	
			}
			inthi = round(ceil(19999*unifs_hi[j]));
			if(inthi > 19999){
				inthi = 19999;	
			}
		if(   (truebx[0] > bxs[intlo]) & (truebx[0] < bxs[inthi])   ){
			plausestrux[0] = plausestrux[0]+0.00005;
		}	
			intlo = round(floor(19999*unifs_loz[j]));
			if(intlo < 0){
				intlo = 0;	
			}
			inthi = round(ceil(19999*unifs_hiz[j]));
			if(inthi > 19999){
				inthi = 19999;	
			}
		if(   (truebz[0] > bzs[intlo]) & (truebz[0] < bzs[inthi])   ){
			plausestruz[0] = plausestruz[0]+0.00005;
		}
	}
		
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("samples") = samples);		
	}
	return result;
	
}




Rcpp::List plauscontourGF2(NumericVector par, NumericVector stat, NumericVector del, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector bzseq) {

	List result;
	Rcpp::Function sortmat("sortmat");
	
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
	
	
	mat J(5,5,fill::zeros); mat J2(3,3,fill::zeros);
	
	J = { { dv1dbx[0], dv1dbz[0], dv1dmux[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dbz[0], dv2dmux[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dbz[0], dv3dmux[0], dv3dsx[0], dv3dse[0] },
            { dz1dbx[0], dz1dbz[0], dz1dmux[0], dz1dsx[0], dz1dse[0] },
            { dz2dbx[0], dz2dbz[0], dz2dmux[0], dz2dsx[0], dz2dse[0] } };
		
	J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
            { dv2dbx[0], dv2dsx[0], dv2dse[0] },
            { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
	NumericVector detJ(1,0.0); NumericVector detJ2(1,0.0); 
	
	detJ[0] = log(std::abs(arma::det(J)));
	detJ2[0] = log(std::abs(arma::det(J2)));
	
	NumericVector logdens(1, 0.0); NumericVector logdens2(1, 0.0);
	
	logdens[0] = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
	logdens2[0] = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
	
	tempdens[0] = logdens[0]; 
	//  Begin MCMC  
	
	NumericVector zeroes(700000,0.0);
	NumericMatrix samples = NumericMatrix(100000, 7, zeroes.begin());
	NumericMatrix samples2 = NumericMatrix(100000, 7, zeroes.begin());
	samples(0,0) = bx[0];samples(0,1) = bz[0];samples(0,2) = mux[0];samples(0,3) = sx[0];samples(0,4) = se[0];samples(0,5) = logdens[0];samples(0,6) = logdens2[0];
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = logdens[0];
	NumericVector currdens2(1,0.0);NumericVector propdens2(1,0.0); propdens2[0] = logdens2[0];
	
	// plauscontourGF2 eliminates thinning of the mcmc samples	
	
	for(int j=0; j<200000; j++) {
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
				L22[0] = std::sqrt((sx[0]/del[0]) - (L12[0]*L12[0]));
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
				
				J2 = { { dv1dbx[0], dv1dsx[0], dv1dse[0] },
			            { dv2dbx[0],dv2dsx[0], dv2dse[0] },
  			 	         { dv3dbx[0], dv3dsx[0], dv3dse[0] }};
	
				detJ[0] = log(std::abs(arma::det(J))); detJ2[0] = log(std::abs(arma::det(J2)));
				if((v1[0]>0) & (v3[0]>0)){
					currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,sd[0],1) + R::dnorm(z2[0],0.0,sd[0],1);
					currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				}else {
					currdens[0]  = -1000000;
					currdens2[0]  = -1000000;
				}
			}else {
				currdens[0]  = -1000000;
				currdens2[0]  = -1000000;
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
				propdens2[0] = currdens2[0];
				ct[i] = ct[i]+1.0;
			}
		}
		if(j < 100000) {
			for(int i=0; i<5; i++){
				samples(j, i) = propsamp[i];
			}
			samples(j, 5) = propdens[0];
			samples(j, 6) = propdens2[0];
		}
		if(j >= 100000){
			for(int i=0; i<5; i++){
				samples2((j-100000), i) = propsamp[i];

			}
			samples2((j-100000), 5) = propdens[0];
			samples2((j-100000), 6) = propdens2[0];
		}
	}
	for(int i=0; i<5; i++){
		ct[i]=ct[i]/100000.0;
	}
		
		
	// plausibility of betax
		
	NumericVector plausestrux(1, 0.0);
	samples = sortmat(samples,6);
	NumericVector plausbxseq(500,0.0);
	NumericVector randsetlo(100000,0.0);NumericVector randsethi(100000,0.0);
	int ind = 0;	
		
	for(int i = 0; i < 100000; i++){
		ind = 0;
		if(samples2(i,6) > samples(0,6)){
			if(samples2(i,6) < samples(99999,6)){
				while(samples2(i,6) > samples(ind,6)){
					ind = ind+1;
				}
				NumericVector subset(100000-ind,0.0);
				for(int j = ind; j < 100000; j++){
					subset[j-ind] = samples(j,0);
				}
				randsetlo[i] = 	Rcpp::min(subset);
				randsethi[i] =  Rcpp::max(subset);
			}else {
				randsetlo[i] = 	0;
				randsethi[i] =  0;
			}
		}else {
			NumericVector subset(100000,0.0);
			for(int j = 0; j < 100000; j++){
				subset[j] = samples(j,0);
			}
			randsetlo[i] = 	Rcpp::min(subset);
			randsethi[i] =  Rcpp::max(subset);
		}
		for(int j = 0; j < 500; j++){
			if((bxseq[j] > randsetlo[i]) & (bxseq[j] < randsethi[i])){
				plausbxseq[j] = plausbxseq[j] + 0.00001;	
			}
		}
		if((truebx[0] > randsetlo[i]) & (truebx[0] < randsethi[i])){
			plausestrux[0] = plausestrux[0] + 0.00001;	
		}
	}
		
	NumericVector plausestruz(1, 0.0);
	samples = sortmat(samples,5);
	NumericVector plausbzseq(500,0.0);
	NumericVector randsetloz(100000,0.0);NumericVector randsethiz(100000,0.0);
	ind = 0;	
		
	for(int i = 0; i < 100000; i++){
		ind = 0;
		if(samples2(i,5) > samples(0,5)){
			if(samples2(i,5) < samples(99999,5)){
				while(samples2(i,5) > samples(ind,5)){
					ind = ind+1;
				}
				NumericVector subset(100000-ind,0.0);
				for(int j = ind; j < 100000; j++){
					subset[j-ind] = samples(j,1);
				}
				randsetloz[i] = Rcpp::min(subset);
				randsethiz[i] =  Rcpp::max(subset);
			}else {
				randsetloz[i] = 0;
				randsethiz[i] =  0;
			}
		}else {
			NumericVector subset(100000,0.0);
			for(int j = 0; j < 100000; j++){
				subset[j] = samples(j,1);
			}
			randsetloz[i] = Rcpp::min(subset);
			randsethiz[i] =  Rcpp::max(subset);
		}
		for(int j = 0; j < 500; j++){
			if((bzseq[j] > randsetloz[i]) & (bzseq[j] < randsethiz[i])){
				plausbzseq[j] = plausbzseq[j] + 0.00001;	
			}
		}
		if((truebz[0] > randsetloz[i]) & (truebz[0] < randsethiz[i])){
			plausestruz[0] = plausestruz[0] + 0.00001;	
		}
	}
		
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_x") = plausbxseq, Rcpp::Named("plauses_beta_z") = plausbzseq, Rcpp::Named("randsetlo") = randsetlo, Rcpp::Named("randsethi") = randsethi, Rcpp::Named("samples") = samples, Rcpp::Named("samples2") = samples2);		
	}
	return result;
}

	


Rcpp::List plauscontourIM(NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz, NumericVector bxseq, NumericVector sxseq, NumericVector seseq) {
	
	Rcpp::Function sortmat("sortmat");
	
	List result;
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	NumericVector check(1,0.0);
	
	NumericVector L11(1,0.0); 
	NumericVector L12(1,0.0); 
	NumericVector L22(1,1.0); 

	NumericVector v1(1, 0.0); 
	NumericVector v2(1, 0.0); 
	NumericVector v3(1, 0.0); 
	NumericVector u(1, 0.0); 
	NumericVector logdensseq(1, 0.0); 

	NumericVector plausbetax(500, 0.0);
	NumericVector plausbetaxseq(1, 0.0);
	
	// Generate MC sample of aux rvs
	
	NumericVector zeroes = NumericVector(10000*4, 0.0); 
        NumericMatrix samps = NumericMatrix(10000, 4, zeroes.begin());
	NumericVector V2(10000,0.0); V2 = Rcpp::rnorm( 10000, 0.0, 1.0 );
	NumericVector V1(10000,0.0); V1 = Rcpp::rchisq( 10000, n[0]-1 );
	NumericVector V3(10000,0.0); V3 = Rcpp::rchisq( 10000, n[0]-2 );
	NumericVector U(10000,0.0);
	NumericVector logdens(10000,0.0);
	for(int i=0; i < 10000; i++){
		V1[i] = std::sqrt(V1[i]);	
		V3[i] = std::sqrt(V3[i]);
		U[i] = V2[i]/V3[i];
		logdens[i] = log(V3[i]) + R::dchisq(V1[i]*V1[i],n[0]-1,true) + R::dchisq(V3[i]*V3[i],n[0]-2,true) + R::dnorm(U[i]*V3[i],0.0,1.0,true);
		samps(i,0) = logdens[i];samps(i,1) = V1[i];samps(i,2) = V3[i];samps(i,3) = U[i];
	}
	samps = sortmat(samps,0);

	// Computing plausibility contour of beta_x using grid of variance components and MC density random set
	
	int ind=0;
	
	for(int i=0; i < 500; i++){
		for(int j=0; j < 500; j++){
			for(int k=0; k < 500; k++){
				L11[0] = std::sqrt(seseq[k]+sxseq[j]*bxseq[i]*bxseq[i]);
				L12[0] = sxseq[j]*bxseq[i]/L11[0];
				if(type[0] == 1.0){
					check[0] = (sxseq[j]/del[0]) - (L12[0]*L12[0]);
				}else if(type[0] == 2.0){
					check[0] = (sxseq[j] + seseq[k]/del[0]) - (L12[0]*L12[0]);	
				}else if(type[0] == 3.0){
					check[0] = (sxseq[j] + del[0]) - (L12[0]*L12[0]);	
				}else {
					check[0] = -1;	
				}	
				if(check[0] > 0){
					if(type[0] == 1.0){
						L22[0] = std::sqrt(sxseq[j]/del[0] - L12[0]*L12[0]);
					}else if(type[0] == 2.0){
						L22[0] = std::sqrt((sxseq[j] + seseq[k]/del[0]) - (L12[0]*L12[0]));	
					}else {
						L22[0] = std::sqrt((sxseq[j] + del[0]) - (L12[0]*L12[0]));
					}
					v1[0] = s11[0]/L11[0];
					v2[0] = (s12[0] - v1[0]*L12[0])/L22[0];
					v3[0] = s22[0]/L22[0];
					u[0] = v2[0]/v3[0];
					logdensseq[0] = log(v3[0]) + R::dchisq(v1[0]*v1[0],n[0]-1,true) + R::dchisq(v3[0]*v3[0],n[0]-2,true) + R::dnorm(u[0]*v3[0],0.0,1.0,true);
					if(logdensseq[0] < samps(0,0)){
						plausbetaxseq[0] = 0.0;	
					}else if(logdensseq[0] > samps(9999,0)){
						plausbetaxseq[0] = 1.0;	
					}else {
						while((ind < 10000) & (logdensseq[0] > samps(ind,0))){
							ind = ind + 1;
							plausbetaxseq[0] = plausbetaxseq[0] + 0.0001;	
						}
						ind = 0;
					}
				}
				if(plausbetaxseq[0] > plausbetax[i]){
					plausbetax[i] = plausbetaxseq[0];
				}
				plausbetaxseq[0] = 0.0;
			}	
		}
	}
	
	
	// Computing plausibility of true beta_x using grid of variance components and MC density random set
	
	NumericVector plaustruebetax(1,0.0);
	
		for(int j=0; j < 500; j++){
			for(int k=0; k < 500; k++){
				L11[0] = std::sqrt(seseq[k]+sxseq[j]*truebx[0]*truebx[0]);
				L12[0] = sxseq[j]*truebx[0]/L11[0];
				if(type[0] == 1.0){
					check[0] = (sxseq[j]/del[0]) - (L12[0]*L12[0]);
				}else if(type[0] == 2.0){
					check[0] = (sxseq[j] + seseq[k]/del[0]) - (L12[0]*L12[0]);	
				}else if(type[0] == 3.0){
					check[0] = (sxseq[j] + del[0]) - (L12[0]*L12[0]);	
				}else {
					check[0] = -1;	
				}	
				if(check[0] > 0){
					if(type[0] == 1.0){
						L22[0] = std::sqrt(sxseq[j]/del[0] - L12[0]*L12[0]);
					}else if(type[0] == 2.0){
						L22[0] = std::sqrt((sxseq[j] + seseq[k]/del[0]) - (L12[0]*L12[0]));	
					}else {
						L22[0] = std::sqrt((sxseq[j] + del[0]) - (L12[0]*L12[0]));
					}
					v1[0] = s11[0]/L11[0];
					v2[0] = (s12[0] - v1[0]*L12[0])/L22[0];
					v3[0] = s22[0]/L22[0];
					u[0] = v2[0]/v3[0];
					logdensseq[0] = log(v3[0]) + R::dchisq(v1[0]*v1[0],n[0]-1,true) + R::dchisq(v3[0]*v3[0],n[0]-2,true) + R::dnorm(u[0]*v3[0],0.0,1.0,true);
					if(logdensseq[0] < samps(0,0)){
						plausbetaxseq[0] = 0.0;	
					}else if(logdensseq[0] > samps(9999,0)){
						plausbetaxseq[0] = 1.0;	
					}else {
						while((ind < 10000) & (logdensseq[0] > samps(ind,0))){
							ind = ind + 1;
							plausbetaxseq[0] = plausbetaxseq[0] + 0.0001;	
						}
						ind = 0;
					}
				}
				if(plausbetaxseq[0] > plaustruebetax[0]){
					plaustruebetax[0] = plausbetaxseq[0];
				}
				plausbetaxseq[0] = 0.0;
			}	
		}
	
	
	
	result = Rcpp::List::create(Rcpp::Named("plaus_beta_x") = plaustruebetax, Rcpp::Named("plauses_beta_x") = plausbetax);
	
	return result;
	
}
	

Rcpp::List plauscontourIMmarg(NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz, NumericVector bxseq) {
	
	List result;
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector check(1,0.0);
	
	NumericVector L11(1,0.0); 
	NumericVector L12(1,0.0); 
	NumericVector L22(1,1.0); 


	NumericVector plausbetax(500, 0.0);
	
	// Generate MC sample of aux rvs
	
        NumericVector sampslo(10000, 0.0);
	NumericVector sampshi(10000, 0.0);
	NumericVector V2(10000,0.0); V2 = Rcpp::rnorm( 10000, 0.0, 1.0 );
	NumericVector V1(10000,0.0); V1 = Rcpp::rchisq( 10000, n[0]-1 );
	//NumericVector V3(10000,0.0); V3 = Rcpp::rchisq( 10000, n[0]-2 );
	NumericVector U(10000,0.0);
	for(int i=0; i < 10000; i++){
		V1[i] = std::sqrt(V1[i]);	
		//V3[i] = std::sqrt(V3[i]);
		U[i] = V2[i]/V1[i];
		sampslo[i] = fmin(U[i]*std::sqrt((1/del[0])-1.0),U[i]*std::sqrt((1/del[0])*0.01-0.0001));
		sampshi[i] = 1.0 + fmax(U[i]*std::sqrt((1/del[0])-1.0), U[i]*std::sqrt((1/del[0])*0.01-0.0001));
	}
	std::sort(sampslo.begin(), sampslo.end());
	std::sort(sampshi.begin(), sampshi.end());
	
	


	// Computing plausibility contour of beta_x using grid of variance components and MC density random set
	
	NumericVector uni(1,0.0);
	NumericVector unilo(1,0.0);
	NumericVector unihi(1,0.0);
	NumericVector betaxlo(1,0.0);
	NumericVector betaxhi(1,0.0);
	int indlo;indlo=0;
	int indhi;indhi=0;
	NumericVector plaustruebetax(1,0.0);
	
	
	for(int j=0; j < 10000; j++){
			uni[0] = R::runif(0.0,1.0); 
			unilo[0] = 0.5 - std::abs(0.5-uni[0]); 
			unihi[0] = 1.0 - unilo[0];
			indlo = floor(unilo[0]*9999);
			indhi = ceil(unihi[0]*9999);
			betaxlo[0] = fmin(sampslo[indlo]/(s12[0]/s11[0]), sampshi[indhi]/(s12[0]/s11[0]));
			betaxhi[0] = fmax(sampslo[indlo]/(s12[0]/s11[0]), sampshi[indhi]/(s12[0]/s11[0]));
		for(int i=0; i < 500; i++){
			if((bxseq[i] > betaxlo[0]) & (bxseq[i] < betaxhi[0])){
				plausbetax[i] = plausbetax[i] + 0.0001;
			}
		}
		if((truebx[0] > betaxlo[0]) & (truebx[0] < betaxhi[0])){
			plaustruebetax[0] = plaustruebetax[0] + 0.0001;
		}
	}
		
	result = Rcpp::List::create(Rcpp::Named("plaus_beta_x") = plaustruebetax, Rcpp::Named("plauses_beta_x") = plausbetax);
	
	return result;
	
	
}
	

Rcpp::NumericMatrix sortmat(NumericMatrix x, unsigned int col){
  arma::mat xx = as<arma::mat>(x);
  arma::uvec id = arma::sort_index(xx.col(col));
  
  for(unsigned int i = 0; i<xx.n_cols; i++){
    arma::vec sub = xx.col(i);
    xx.col(i) = sub.elem(id);
  }
  
  NumericMatrix	y = wrap(xx); 
  return y;
}
