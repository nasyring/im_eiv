#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include <Rcpp.h>
#include <math.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;
#include <cmath>
#include <algorithm>


Rcpp::List plauscontourMCMC(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector truebz, NumericVector bzseq, NumericVector randsettype) {

	List result;
	Rcpp::Function sortmat("sortmat");
	
	NumericVector check(1, 0.0);
	NumericVector tempdens(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector ct(3,0.0);
	NumericVector bx(1,0.0); bx[0] = par[0];
	NumericVector bz(1,0.0); bz[0] = par[1];
	NumericVector mux(1,0.0); mux[0] = par[2];
	NumericVector sx(1,0.0); sx[0] = par[3];
	NumericVector se(1,0.0); se[0] = par[4];
	
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
	NumericVector dv1dsx(1, 0.0); dv1dsx[0] = -s11[0]*dL11dsx[0]/(L11[0]*L11[0]);
	NumericVector dv1dse(1, 0.0); dv1dse[0] = -s11[0]*dL11dse[0]/(L11[0]*L11[0]);
	NumericVector dv1dbz(1, 0.0);
	NumericVector dv1dmux(1, 0.0);
	
	NumericVector dv2dbx(1, 0.0); dv2dbx[0] = ((-L12[0]*dv1dbx[0]-v1[0]*dL12dbx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dbx[0])/(L22[0]*L22[0]);
	NumericVector dv2dsx(1, 0.0); dv2dsx[0] = ((-L12[0]*dv1dsx[0]-v1[0]*dL12dsx[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dsx[0])/(L22[0]*L22[0]);
	NumericVector dv2dse(1, 0.0); dv2dse[0] = ((-L12[0]*dv1dse[0]-v1[0]*dL12dse[0])*L22[0] - (s12[0]-v1[0]*L12[0])*dL22dse[0])/(L22[0]*L22[0]);
	NumericVector dv2dbz(1, 0.0);
	NumericVector dv2dmux(1, 0.0);
		
	NumericVector dv3dbx(1, 0.0); dv3dbx[0] = -s22[0]*dL22dbx[0]/(L22[0]*L22[0]);
	NumericVector dv3dsx(1, 0.0); dv3dsx[0] = -s22[0]*dL22dsx[0]/(L22[0]*L22[0]);
	NumericVector dv3dse(1, 0.0); dv3dse[0] = -s22[0]*dL22dse[0]/(L22[0]*L22[0]);
	NumericVector dv3dbz(1, 0.0);
	NumericVector dv3dmux(1, 0.0);
		
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
	
	
	detJ[0] = log(std::abs(arma::det(J))); detJ2[0] = log(std::abs(arma::det(J2)));
	
	NumericVector logdens(1, 0.0); NumericVector logdens2(1, 0.0);
	
	logdens[0] = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,std::sqrt(1.0/n[0]),1) + R::dnorm(z2[0],0.0,std::sqrt(1.0/n[0]),1);
	logdens2[0] = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
	
	tempdens[0] = logdens[0]; 
	//  Begin MCMC  
	
	NumericVector zeroes(280000,0.0);
	NumericMatrix samples = NumericMatrix(40000, 7, zeroes.begin());
	samples(0,0) = bx[0];samples(0,1) = bz[0];samples(0,2) = mux[0];samples(0,3) = sx[0];samples(0,4) = se[0];samples(0,5) = logdens[0];samples(0,6) = logdens2[0];
	NumericVector densdiff(1,0.0);
	NumericVector propsamp(5,0.0);
	propsamp[0] = bx[0];propsamp[1] = bz[0];propsamp[2] = mux[0];propsamp[3] = sx[0];propsamp[4] = se[0];
	NumericVector currsamp(5,0.0);
	currsamp[0] = bx[0];currsamp[1] = bz[0];currsamp[2] = mux[0];currsamp[3] = sx[0];currsamp[4] = se[0];
	NumericVector currdens(1,0.0);NumericVector propdens(1,0.0); propdens[0] = logdens[0];
	NumericVector currdens2(1,0.0);NumericVector propdens2(1,0.0); propdens2[0] = logdens2[0];

	
	for(int j=0; j<400000; j++) {
		for(int i=0; i<5; i++){
			for(int k=0; k<5; k++){
				currsamp[k] = propsamp[k];
			}
			currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			if(!((i > 0) & (currsamp[i]<0.0))){
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
					currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,std::sqrt(1.0/n[0]),1) + R::dnorm(z2[0],0.0,std::sqrt(1.0/n[0]),1);
					currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				}
				uu[0] = R::runif(0.0,1.0);
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);
			}else {
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
		ct[i]=ct[i]/400000.0;
	}
			
		
	// plausibility
	NumericVector plausesx(500,0.0);
	NumericVector plausestrux(1,0.0);
	NumericVector plausesz(500,0.0);
	NumericVector plausestruz(1,0.0);
	int dim = round(randsettype[0]);
	NumericVector unifs_hi(40000,0.0);NumericVector unifs_lo(40000,0.0);
	NumericVector bxs(40000,0.0);NumericVector bzs(40000,0.0);
	NumericVector unifs(1,0.0);NumericVector maxunifs(1,0.0);
	if(randsettype[0] > 0.0){	
		for(int i=0; i<40000; i++){
			for(int j = 0; j < dim; j++){
				unifs[0] = R::runif(0.0,1.0);
				maxunifs[0] = std::max(maxunifs[0], unifs[0]);	
			}
			unifs_hi[i] = 0.5 + std::abs(maxunifs[0] - 0.5); 
			unifs_lo[i] = 1.0-unifs_hi[i];
			maxunifs[0] = 0.0;
			bxs[i] = samples(i, 0);	bzs[i] = samples(i, 1);
		}
		std::sort(bxs.begin(), bxs.end());std::sort(bzs.begin(), bzs.end());
		int intlo = 0;  int inthi = 0;
		for(int j=0; j<40000; j++){
			intlo = int(floor(39999*unifs_lo[j]));
			inthi = int(ceil(39999*unifs_hi[j]));
			if(   (truebx[0] > bxs[intlo]) & (truebx[0] < bxs[inthi])   ){
				plausestrux[0] = plausestrux[0]+0.000025;
			}	
			if(   (truebz[0] > bzs[intlo]) & (truebz[0] < bzs[inthi])   ){
				plausestruz[0] = plausestruz[0]+0.000025;
			}	
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > bxs[intlo]) & (bxseq[i] < bxs[inthi])   ){
					plausesx[i] = plausesx[i]+0.000025;
				}
				if(   (bzseq[i] > bzs[intlo]) & (bzseq[i] < bzs[inthi])   ){
					plausesz[i] = plausesz[i]+0.000025;
				}
			}
		}
	}else {
		samples = sortmat(samples,6);
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			NumericVector subset(40000-j-1, 0.0);
			for(int i=0; i<(40000-j-1); i++){
				subset[i] = samples(i+j+1,0);	
			}
			randsetslo[0] = Rcpp::min(subset);randsetshi[0] = Rcpp::max(subset);
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > randsetslo[0]) & (bxseq[i] < randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/39999.0);
				}
			}
		}
		samples = sortmat(samples,5);
		for(int j=0; j<39999; j++){
			NumericVector subset(40000-j-1, 0.0);
			for(int i=0; i<(40000-j-1); i++){
				subset[i] = samples(i+j+1,1);	
			}
			randsetslo[0] = Rcpp::min(subset);randsetshi[0] = Rcpp::max(subset);
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bzseq[i] > randsetslo[0]) & (bzseq[i] < randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1/39999.0);
				}
			}
		}
		
	}	
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("samples") = samples,Rcpp::Named("unifs_hi") = unifs_hi,Rcpp::Named("unifs_lo") = unifs_lo,Rcpp::Named("bxs") = bxs,Rcpp::Named("bzs") = bzs, Rcpp::Named("dim") = dim);		
	}
	return result;
	
}


Rcpp::List plauscontourMC(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq,NumericVector truebz, NumericVector bzseq, NumericVector randsettype) {

	List result;
	Rcpp::Function sortmat("sortmat");
	int size = round(sampsize[0]);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
	// Generate MC sample of aux rvs

	NumericVector V2(1,0.0); 
	NumericVector V1(1,0.0); 
	NumericVector V3(1,0.0); 
	NumericVector Z1(1,0.0); 
	NumericVector Z2(1,0.0); 
	NumericVector zeroes(2*size,0.0);
	NumericMatrix bxs(size,2,zeroes.begin()); NumericMatrix bzs(size,2,zeroes.begin());
	NumericVector bx(1,0.0); NumericVector sx(1,0.0); NumericVector se(1,0.0); NumericVector mux(1,0.0); NumericVector bz(1,0.0); 
	int ind = 0;
	while(ind < size){
		V2[0] = R::rnorm( 0.0, 1.0 );V1[0] = R::rchisq( n[0]-1 );V3[0] = R::rchisq( n[0]-2 );Z1[0] = R::rnorm( 0.0, std::sqrt(1.0/n[0]) );Z2[0] = R::rnorm( 0.0, std::sqrt(1.0/n[0]) );
		V1[0] = std::sqrt(V1[0]);	
		V3[0] = std::sqrt(V3[0]);
		if(type[0] == 2.0){
			NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
			L11[0] = s11[0]/V1[0]; L22[0] = s22[0]/V3[0]; L12[0] = (s12[0] - V2[0]*L22[0])/V1[0]; 
			sx[0] = 0.5*(-(L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0]) + std::sqrt(((L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0])*(L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0]))+4*L11[0]*L11[0]*L12[0]*L12[0]/del[0]));
			bx[0] = L11[0]*L12[0]/sx[0];
			se[0] = L11[0]*L11[0]-sx[0]*bx[0]*bx[0];
			mux[0] = wbar[0] - L12[0]*Z1[0]-L22[0]*Z2[0];
			bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[0];
			if((se[0] > 0.0) & (sx[0]>0.0)){
				bxs(ind,0) = bx[0]; 
				bzs(ind,0) = bz[0]; 
				bxs(ind,1) = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);	
				bzs(ind,1) = bxs(ind,1) + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
				ind = ind+1;
			}
		}else if(type[0] == 1.0){
			NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
			L11[0] = s11[0]/V1[0]; L22[0] = s22[0]/V3[0]; L12[0] = (s12[0] - V2[0]*L22[0])/V1[0]; 
			sx[0] = del[0]*(L22[0]*L22[0]+L12[0]*L12[0]);
			bx[0] = L11[0]*L12[0]/sx[0];
			se[0] = L11[0]*L11[0]-sx[0]*bx[0]*bx[0];
			mux[0] = wbar[0] - L12[0]*Z1[0]-L22[0]*Z2[0];
			bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[0];
			if((se[0] > 0.0) & (sx[0]>0.0)){
				bxs(ind,0) = bx[0]; 
				bzs(ind,0) = bz[0]; 
				bxs(ind,1) = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);	
				bzs(ind,1) = bxs(ind,1) + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
				ind = ind+1;
			}		
		}else {
			NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
			L11[0] = s11[0]/V1[0]; L22[0] = s22[0]/V3[0]; L12[0] = (s12[0] - V2[0]*L22[0])/V1[0]; 
			sx[0] = L22[0]*L22[0]+L12[0]*L12[0] - del[0];
			bx[0] = L11[0]*L12[0]/sx[0];
			se[0] = L11[0]*L11[0]-sx[0]*bx[0]*bx[0];
			mux[0] = wbar[0] - L12[0]*Z1[0]-L22[0]*Z2[0];
			bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[0];
			if((se[0] > 0.0) & (sx[0]>0.0)){
				bxs(ind,0) = bx[0]; 
				bzs(ind,0) = bz[0]; 
				bxs(ind,1) = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);	
				bzs(ind,1) = bxs(ind,1) + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
				ind = ind+1;
			}
		}
	}

	NumericVector zeroes2(2*ind,0.0);
	NumericMatrix samples_bx(ind,2,zeroes.begin()); NumericMatrix samples_bz(ind,2,zeroes.begin());
	NumericVector bx_s(size, 0.0); NumericVector bz_s(size, 0.0);
	for(int i=0; i < ind; i++){
		bx_s[i] = bxs(i,0); samples_bx(i,0) = bxs(i,0); samples_bx(i,1) = bxs(i,1);
		bz_s[i] = bzs(i,0); samples_bz(i,0) = bzs(i,0); samples_bz(i,1) = bzs(i,1);
	}
	samples_bx = sortmat(samples_bx,1);
	samples_bz = sortmat(samples_bz,1);
	std::sort(bx_s.begin(), bx_s.end()); std::sort(bz_s.begin(), bz_s.end());
	

	

	// plausibility
	NumericVector plausesx(500,0.0);
	NumericVector plausestrux(1,0.0);
	NumericVector plausesz(500,0.0);
	NumericVector plausestruz(1,0.0);
	
	int inc = round(floor(size / 500));
	NumericVector bx_seq(500, 0.0);NumericVector bz_seq(500, 0.0);
	for(int i = 0; i < 500; i++){
		bx_seq[i] = bx_s[i*inc];bz_seq[i] = bz_s[i*inc];	
	}
	bx_seq[499] = bx_s[(size-1)];bz_seq[499] = bz_s[(size-1)];
	int dim = round(randsettype[0]);
	
	if(randsettype[0] == 0){


		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<(ind-1); j++){
			NumericVector subset(ind-j-1, 0.0);
			for(int i=0; i<(ind-j-1); i++){
				subset[i] = samples_bx(i+j+1,0);	
			}
			if(ind-j-1>1){
			std::sort(subset.begin(), subset.end());
			}
			randsetslo[0] = subset[0]; randsetshi[0] = subset[ind-j-2];
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/(ind-1.0));
			}
			for(int i=0; i<500; i++){
				if(   (bx_seq[i] >= randsetslo[0]) & (bx_seq[i] <= randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/(ind-1.0));
				}
			}
		}



		for(int j=0; j<(ind-1); j++){
			NumericVector subset(ind-j-1, 0.0);
			for(int i=0; i<(ind-j-1); i++){
				subset[i] = samples_bz(i+j+1,0);	
			}
			if(ind-j-1>1){
			std::sort(subset.begin(), subset.end());
			}
			randsetslo[0] = subset[0]; randsetshi[0] = subset[ind-j-2];
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/(ind-1.0));
			}
			for(int i=0; i<500; i++){
				if(   (bz_seq[i] >= randsetslo[0]) & (bz_seq[i] <= randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1.0/(ind-1.0));
				}
			}
		}
	}else {
	NumericVector unifs_hi(size,0.0);NumericVector unifs_lo(size,0.0);
	NumericVector bxs(size,0.0);NumericVector bzs(size,0.0);
	NumericVector unifs(1,0.0);NumericVector maxunifs(1,0.0);
		for(int i=0; i<size; i++){
			for(int j = 0; j < dim; j++){
				unifs[0] = R::runif(0.0,1.0);
				maxunifs[0] = std::max(maxunifs[0], unifs[0]);	
			}
			unifs_hi[i] = 0.5 + std::abs(maxunifs[0] - 0.5); 
			unifs_lo[i] = 1.0-unifs_hi[i];
			maxunifs[0] = 0.0;
		}
		int intlo = 0;  int inthi = 0;
		for(int j=0; j<size; j++){
			intlo = int(floor((size - 1)*unifs_lo[j]));
			inthi = int(ceil((size-1)*unifs_hi[j]));
			if(   (truebx[0] > bx_s[intlo]) & (truebx[0] < bx_s[inthi])   ){
				plausestrux[0] = plausestrux[0]+(1.0/size);
			}	
			if(   (truebz[0] > bz_s[intlo]) & (truebz[0] < bz_s[inthi])   ){
				plausestruz[0] = plausestruz[0]+(1.0/size);
			}	
			for(int i=0; i<500; i++){
				if(   (bx_seq[i] >= bx_s[intlo]) & (bx_seq[i] <= bx_s[inthi])   ){
					plausesx[i] = plausesx[i]+(1.0/size);
				}
				if(   (bz_seq[i] >= bz_s[intlo]) & (bz_seq[i] <= bz_s[inthi])   ){
					plausesz[i] = plausesz[i]+(1.0/size);
				}
			}
		}
		
		
	}
	
	result = Rcpp::List::create(Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plauses_beta_x") = plausesx,  Rcpp::Named("samples_bx") = samples_bx, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_z") = plausesz,  Rcpp::Named("samples_bz") = samples_bz, Rcpp::Named("bx_seq") = bx_seq, Rcpp::Named("bz_seq") = bz_seq);		

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
