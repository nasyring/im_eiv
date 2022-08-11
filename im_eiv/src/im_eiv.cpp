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


Rcpp::List plausMC(NumericVector theta, NumericVector intcpt, NumericMatrix grid, NumericVector stat, NumericVector del, NumericVector df, int m_samps, bool intercept){

	List result;
	int m_the = theta.length();
	int m_grid = grid.nrow();
	int m_int = intcpt.length();
	
	NumericVector L11(1,0.0);
	NumericVector L12(1,0.0);
	NumericVector L22(1,0.0);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	

	NumericVector V1(m_samps, 0.0); V1 = Rcpp::rchisq(m_samps,df[0]);
	NumericVector V2(m_samps, 0.0); V2 = Rcpp::rnorm(m_samps,0.0,1.0);	
	NumericVector V3(m_samps, 0.0); V3 = Rcpp::rchisq(m_samps,df[1]);	
	
	NumericVector plaus_theta_temp(1, 0.0);
	NumericVector plaus_theta(m_the, 0.0);
	
	
	NumericVector dens1(1, 0.0);
	NumericVector dens2(1, 0.0);
	NumericVector dens3(1, 0.0);
	NumericVector dens4(1, 0.0);
	
	bool marginalize = TRUE;
	for(int i = 0; i < m_samps; i++){
		if(std::pow((1/std::sqrt(V1[i]))*(s12[0]-s22[0]*V2[i]/std::sqrt(V3[i])), 2.0)>1.0){
			marginalize = FALSE;	
		}
	}
	
	
	if(marginalize){
		if(intercept){
			NumericVector plaus_intcpt(m_int, 0.0);
			NumericVector Z(m_samps, 0.0); Z = Rcpp::rnorm(m_samps,0.0,1.0/std::sqrt(df[0]+1.0));
			NumericVector temp1(1, 0.0); NumericVector temp2(1, 0.0); NumericVector temp3(1, 0.0);
			NumericVector aux_var(m_samps,0.0);
			NumericVector aux_var2(m_samps,0.0);
			for(int k = 0; k< m_samps; k++){
				temp1[0] = (1/std::sqrt(V1[k]))*(s12[0]-s22[0]*V2[k]/std::sqrt(V3[k]));
				temp2[0] = del[0]*(std::pow(s22[0]/std::sqrt(V3[k]),2.0)+std::pow(temp1[0],2.0));
				aux_var[k] = (s11[0]/std::sqrt(V1[k]))*temp1[0]/temp2[0];
				aux_var2[k] = stat[3] - aux_var[k]*stat[4] - Z[k]*std::sqrt(std::pow(aux_var[k]*temp1[0]+(s11[0]/std::sqrt(V1[k])), 2.0) + std::pow(aux_var[k]*(s22[0]/std::sqrt(V3[k])), 2.0));
			}
			NumericVector thetapts(101,0.0);
			NumericVector thetaplaus(101,0.0);
			std::sort(aux_var.begin(), aux_var.end());
			for(int i = 0; i < 100; i++){
				thetapts[i] = aux_var[i*100];
				thetaplaus[i] = 1.0 - std::abs(2.0 * (i/100) - 1.0);
			}
			thetapts[100] = aux_var[9999];
			thetaplaus[0] = 0.0001; thetaplaus[100] = 0.0001;
			NumericVector intpts(101,0.0);
			NumericVector intplaus(101,0.0);
			std::sort(aux_var2.begin(), aux_var2.end());
			for(int i = 0; i < 100; i++){
				intpts[i] = aux_var2[i*100];
				intplaus[i] = 1.0 - std::abs(2.0 * (i/100) - 1.0);
			}
			intpts[100] = aux_var2[9999];
			intplaus[0] = 0.0001; intplaus[100] = 0.0001;
			result = Rcpp::List::create(Rcpp::Named("plauses.theta") = thetaplaus,Rcpp::Named("thetas") = thetapts, Rcpp::Named("plauses.intercept") = intplaus, Rcpp::Named("intercepts") = intpts, Rcpp::Named("marginalize") = marginalize);					      		      
		}else {
			NumericVector temp1(1, 0.0); NumericVector temp2(1, 0.0); NumericVector temp3(1, 0.0);
			NumericVector aux_var(m_samps,0.0);
			for(int k = 0; k< m_samps; k++){
				temp1[0] = (1/std::sqrt(V1[k]))*(s12[0]-s22[0]*V2[k]/std::sqrt(V3[k]));
				temp2[0] = del[0]*(std::pow(s22[0]/std::sqrt(V3[k]),2.0)+std::pow(temp1[0],2.0));
				aux_var[k] = (s11[0]/std::sqrt(V1[k]))*temp1[0]/temp2[0];
			}
			NumericVector thetapts(101,0.0);
			NumericVector thetaplaus(101,0.0);
			std::sort(aux_var.begin(), aux_var.end());
			for(int i = 0; i < 100; i++){
				thetapts[i] = aux_var[i*100];
				thetaplaus[i] = 1.0 - std::abs(2.0 * (i/100) - 1.0);
			}
			thetapts[100] = aux_var[9999];
			thetaplaus[0] = 0.0001; thetaplaus[100] = 0.0001;
			result = Rcpp::List::create(Rcpp::Named("plauses.theta") = thetaplaus,Rcpp::Named("thetas") = thetapts, Rcpp::Named("marginalize") = marginalize, Rcpp::Named("aux.var") = aux_var);					      		      
		}
	}else {
		if(intercept){
			NumericVector plaus_intcpt_temp(1, 0.0);
			NumericVector plaus_intcpt(m_int, 0.0);
			NumericVector Z(m_samps, 0.0); Z = Rcpp::rnorm(m_samps,0.0,1.0/std::sqrt(df[0]+1.0));
			for(int l = 0; l < m_int; l++){	
				for(int i = 0; i < m_the; i++){
					for(int j = 0; j < m_grid; j++){
						L11[0] = std::sqrt( grid(j, 1) + grid(j, 0)*theta[i]*theta[i] );
						L12[0] = grid(j,0)*theta[i]/L11[0];
						L22[0] = (grid(j,0)/del[0]) - std::pow(L12[0], 2.0);
						if(L22[0] > 0.0){
							L22[0] = std::sqrt(L22[0]);
							plaus_theta_temp[0] = 0.0;
							plaus_intcpt_temp[0] = 0.0;
							for(int k = 0; k< m_samps; k++){
								dens1[0] = R::dchisq(V1[k], df[0], 0)*R::dnorm(V2[k],0.0,1.0,0)*R::dchisq(V3[k],df[1],0);
								dens2[0] = R::dchisq(std::pow(s11[0]/L11[0],2.0), df[0], 0)*R::dnorm((s12[0]-L12[0]*s11[0]/L11[0])/L22[0],0.0,1.0,0)*R::dchisq(std::pow(s22[0]/L22[0],2.0),df[1],0);
								dens3[0] = dens1[0]*R::dnorm(Z[k], 0.0, 1.0/std::sqrt(df[0]+1.0), 0);
								dens4[0] = dens2[0]*R::dnorm((ybar[0] - intcpt[l]-theta[i]*wbar[0])/std::sqrt(std::pow(L12[0]*theta[i]+L11[0], 2)+std::pow(theta[i]*L22[0], 2)), 0.0, 1.0/std::sqrt(df[0]+1.0), 0);
								if( dens1[0] <= dens2[0] ){
									plaus_theta_temp[0] = plaus_theta_temp[0] + 1.0/m_samps;	
								}
								if( dens3[0] <= dens4[0] ){
									plaus_intcpt_temp[0] = plaus_intcpt_temp[0] + 1.0/m_samps;	
								}						
							}
							if(plaus_theta[i] < plaus_theta_temp[0]){
								plaus_theta[i] = plaus_theta_temp[0];
							}
							if(plaus_intcpt[l] < plaus_intcpt_temp[0]){
								plaus_intcpt[l] = plaus_intcpt_temp[0];
							}						
						}
					}
				}	
			}
			result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_theta, Rcpp::Named("plauses.intercept") = plaus_intcpt, Rcpp::Named("marginalize") = marginalize);
		}else {
			for(int i = 0; i < m_the; i++){
				for(int j = 0; j < m_grid; j++){
					L11[0] = std::sqrt( grid(j, 1) + grid(j, 0)*theta[i]*theta[i] );
					L12[0] = grid(j,0)*theta[i]/L11[0];
					L22[0] = (grid(j,0)/del[0]) - std::pow(L12[0], 2.0);
					if(L22[0] > 0.0){
						L22[0] = std::sqrt(L22[0]);
						plaus_theta_temp[0] = 0.0;
						for(int k = 0; k< m_samps; k++){
							dens1[0] = R::dchisq(V1[k], df[0], 0)*R::dnorm(V2[k],0.0,1.0,0)*R::dchisq(V3[k],df[1],0);
							dens2[0] = R::dchisq(std::pow(s11[0]/L11[0],2.0), df[0], 0)*R::dnorm((s12[0]-L12[0]*s11[0]/L11[0])/L22[0],0.0,1.0,0)*R::dchisq(std::pow(s22[0]/L22[0],2.0),df[1],0);
							if( dens1[0] <= dens2[0] ){
								plaus_theta_temp[0] = plaus_theta_temp[0] + 1.0/m_samps;	
							}
						}
						if(plaus_theta[i] < plaus_theta_temp[0]){
							plaus_theta[i] = plaus_theta_temp[0];
						}
					}
				}
			}
			result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_theta, Rcpp::Named("marginalize") = marginalize);
		}
	}	
	
	return result;
	
}



Rcpp::List plauscontourMCMC(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector truebz, NumericVector bzseq, NumericVector randsettype) {

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
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			if((se[0]+sx[0]*bx[0]*bx[0] > 0)){
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
			if((check[0] > 0)){	
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

					
				currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,std::sqrt(1.0/n[0]),1) + R::dnorm(z2[0],0.0,std::sqrt(1.0/n[0]),1);
				currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				uu[0] = R::runif(0.0,1.0);
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);
			}else {
				uu[0] = 2.0;
				densdiff[0] = 0;
			}
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
		bool offcheck = true;
		NumericVector offset(1,0.0);
		int index1 = 0;
		while(offcheck){
			if((samples(39999-index1,3)<0) || (samples(39999-index1,4)<0)){
				offset[0] = samples(39999,6) - samples(39999-index1-1,6);	
			}else {
				offcheck = false;	
			}
			index1 = index1+1;
		}
		if(offset[0] > 0.0){
			
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bxs);randsetshi[0] = Rcpp::min(bxs);
			for(int i=0; i<39999; i++){
				if((samples(i,3)>0) & (samples(i,4)>0) & (samples(i,6)>(samples(j,6)-offset[0]))){
					randsetslo[0] = std::min(randsetslo[0], samples(i,0));
					randsetshi[0] = std::max(randsetshi[0], samples(i,0));
				}
			}
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > randsetslo[0]) & (bxseq[i] < randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/39999.0);
				}
			}
		}
			
		}else {
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bxs);randsetshi[0] = Rcpp::min(bxs);
			for(int i=0; i<(40000-j-1); i++){
				if((samples(i+j+1,3)>0) & (samples(i+j+1,4)>0)){
					randsetslo[0] = std::min(randsetslo[0], samples(i+j+1,0));
					randsetshi[0] = std::max(randsetshi[0], samples(i+j+1,0));
				}
			}
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > randsetslo[0]) & (bxseq[i] < randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/39999.0);
				}
			}
		}
		}
		samples = sortmat(samples,5);
		offcheck = true;
		offset[0] = 0.0;
		index1 = 0;
		while(offcheck){
			if((samples(39999-index1,3)<0) || (samples(39999-index1,4)<0)){
				offset[0] = samples(39999,5) - samples(39999-index1-1,5);	
			}else {
				offcheck = false;	
			}
			index1 = index1+1;
		}
		if(offset[0] > 0.0){
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bzs);randsetshi[0] = Rcpp::min(bzs);
			for(int i=0; i<39999; i++){
				if((samples(i,3)>0) & (samples(i,4)>0) & (samples(i,5)>(samples(j,5)-offset[0]))){
					randsetslo[0] = std::min(randsetslo[0], samples(i,1));
					randsetshi[0] = std::max(randsetshi[0], samples(i,1));
				}
			}
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bzseq[i] > randsetslo[0]) & (bzseq[i] < randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1.0/39999.0);
				}
			}
		}
			
		}else {
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bzs);randsetshi[0] = Rcpp::min(bzs);
			for(int i=0; i<(40000-j-1); i++){
				if((samples(i+j+1,3)>0) & (samples(i+j+1,4)>0)){
					randsetslo[0] = std::min(randsetslo[0], samples(i+j+1,1));
					randsetshi[0] = std::max(randsetshi[0], samples(i+j+1,1));
				}
			}
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bzseq[i] > randsetslo[0]) & (bzseq[i] < randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1.0/39999.0);
				}
			}
		}
	}
	}
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("samples") = samples,Rcpp::Named("unifs_hi") = unifs_hi,Rcpp::Named("unifs_lo") = unifs_lo,Rcpp::Named("bxs") = bxs,Rcpp::Named("bzs") = bzs, Rcpp::Named("dim") = dim);		
	}
	return result;
}
								 
Rcpp::List plauscontourMCMC2(NumericVector par, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector propsd, NumericVector truebx, NumericVector bxseq, NumericVector truebz, NumericVector bzseq, NumericVector randsettype) {

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
		for(int i=0; i<4; i++){
			for(int k=0; k<4; k++){
				currsamp[k] = propsamp[k];
			}
			currsamp[i] = R::rnorm(propsamp[i], propsd[i]);
			bx[0] = currsamp[0];bz[0] = currsamp[1];mux[0] = currsamp[2];sx[0] = currsamp[3];se[0] = currsamp[4];
			if((se[0]+sx[0]*bx[0]*bx[0] > 0)){
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
			if((check[0] > 0)){	
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

					
				currdens[0]  = detJ[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1) + R::dnorm(z1[0],0.0,std::sqrt(1.0/n[0]),1) + R::dnorm(z2[0],0.0,std::sqrt(1.0/n[0]),1);
				currdens2[0]  = detJ2[0] + R::dchisq(v1[0]*v1[0],n[0]-1.0,1) +  R::dchisq(v3[0]*v3[0],n[0]-2.0,1) + R::dnorm(v2[0],0.0,1.0,1);
				uu[0] = R::runif(0.0,1.0);
				densdiff[0] = fmin(std::exp(currdens[0] - propdens[0]), 1.0);
			}else {
				uu[0] = 2.0;
				densdiff[0] = 0;
			}
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
		bool offcheck = true;
		NumericVector offset(1,0.0);
		int index1 = 0;
		while(offcheck){
			if((samples(39999-index1,3)<0) || (samples(39999-index1,4)<0)){
				offset[0] = samples(39999,6) - samples(39999-index1-1,6);	
			}else {
				offcheck = false;	
			}
			index1 = index1+1;
		}
		if(offset[0] > 0.0){
			
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bxs);randsetshi[0] = Rcpp::min(bxs);
			for(int i=0; i<39999; i++){
				if((samples(i,3)>0) & (samples(i,4)>0) & (samples(i,6)>(samples(j,6)-offset[0]))){
					randsetslo[0] = std::min(randsetslo[0], samples(i,0));
					randsetshi[0] = std::max(randsetshi[0], samples(i,0));
				}
			}
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > randsetslo[0]) & (bxseq[i] < randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/39999.0);
				}
			}
		}
			
		}else {
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bxs);randsetshi[0] = Rcpp::min(bxs);
			for(int i=0; i<(40000-j-1); i++){
				if((samples(i+j+1,3)>0) & (samples(i+j+1,4)>0)){
					randsetslo[0] = std::min(randsetslo[0], samples(i+j+1,0));
					randsetshi[0] = std::max(randsetshi[0], samples(i+j+1,0));
				}
			}
			if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
				plausestrux[0] = plausestrux[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bxseq[i] > randsetslo[0]) & (bxseq[i] < randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/39999.0);
				}
			}
		}
		}
		samples = sortmat(samples,5);
		offcheck = true;
		offset[0] = 0.0;
		index1 = 0;
		while(offcheck){
			if((samples(39999-index1,3)<0) || (samples(39999-index1,4)<0)){
				offset[0] = samples(39999,5) - samples(39999-index1-1,5);	
			}else {
				offcheck = false;	
			}
			index1 = index1+1;
		}
		if(offset[0] > 0.0){
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bzs);randsetshi[0] = Rcpp::min(bzs);
			for(int i=0; i<39999; i++){
				if((samples(i,3)>0) & (samples(i,4)>0) & (samples(i,5)>(samples(j,5)-offset[0]))){
					randsetslo[0] = std::min(randsetslo[0], samples(i,1));
					randsetshi[0] = std::max(randsetshi[0], samples(i,1));
				}
			}
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bzseq[i] > randsetslo[0]) & (bzseq[i] < randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1.0/39999.0);
				}
			}
		}
			
		}else {
		
		NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
		for(int j=0; j<39999; j++){
			randsetslo[0] = Rcpp::max(bzs);randsetshi[0] = Rcpp::min(bzs);
			for(int i=0; i<(40000-j-1); i++){
				if((samples(i+j+1,3)>0) & (samples(i+j+1,4)>0)){
					randsetslo[0] = std::min(randsetslo[0], samples(i+j+1,1));
					randsetshi[0] = std::max(randsetshi[0], samples(i+j+1,1));
				}
			}
			if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
				plausestruz[0] = plausestruz[0]+(1.0/39999.0);
			}
			for(int i=0; i<500; i++){
				if(   (bzseq[i] > randsetslo[0]) & (bzseq[i] < randsetshi[0])   ){
					plausesz[i] = plausesz[i]+(1.0/39999.0);
				}
			}
		}
	}
	}
	result = Rcpp::List::create(Rcpp::Named("rate") = ct, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plauses_beta_x") = plausesx, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_z") = plausesz, Rcpp::Named("samples") = samples,Rcpp::Named("unifs_hi") = unifs_hi,Rcpp::Named("unifs_lo") = unifs_lo,Rcpp::Named("bxs") = bxs,Rcpp::Named("bzs") = bzs, Rcpp::Named("dim") = dim);		
	}
	return result;
}	


Rcpp::List plauscontourMCMCcond(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector propsd, NumericVector type, NumericVector n, NumericVector truebx, NumericVector truebz,NumericVector bxseq, NumericVector sxseq, NumericVector lenseq, NumericVector plbxseq, NumericVector plbzseq, NumericVector lenplseq, NumericVector se2) {

	List result;
	Rcpp::Function sortmat("sortmat");
	int size = round(sampsize[0]);
	int L = round(lenseq[0]);  	int pL = round(lenplseq[0]);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
	// Generate MC sample of aux rvs Z1, Z2, V2

	NumericVector V1(1,0.0); 
	NumericVector V2(1,0.0);
	NumericVector V3(1,0.0); 
	NumericVector Z1(size,0.0); Z1 = Rcpp::rnorm(size,0.0,std::sqrt(1.0/n[0]));
	NumericVector Z2(size,0.0); Z2 = Rcpp::rnorm(size,0.0,std::sqrt(1.0/n[0]));
	
	// for grid of (bxseq, sxseq) values generate (V1, omega) r.v.'s by MCMC 
	
	NumericVector denscurr(1,0.0); NumericVector densprop(1,0.0);
	NumericVector c1(1,0.0); NumericVector c2(1,0.0); NumericVector L110(1,0.0); NumericVector L220(1,0.0);NumericVector L120(1,0.0); 
	NumericVector dL110bx(1,0.0); NumericVector dL220bx(1,0.0);NumericVector dL110sx(1,0.0); NumericVector dL220sx(1,0.0);NumericVector dL120bx(1,0.0);NumericVector dL120sx(1,0.0);
	NumericVector dV10bx(1,0.0); NumericVector dV30bx(1,0.0);NumericVector dV10sx(1,0.0); NumericVector dV30sx(1,0.0);NumericVector dV20bx(1,0.0); NumericVector dV20sx(1,0.0);
	NumericVector V10(1,0.0); NumericVector V20(1,0.0); NumericVector V30(1,0.0); NumericVector eta(1,0.0);
	NumericVector sampcurr(2,0.0); NumericVector sampprop(2,0.0);
	
	NumericVector zeroes7(size*10, 0.0);
	
	NumericVector maxplausesx(L, 0.0); NumericVector maxplausesz(pL, 0.0); 
	NumericVector offsetx(1,0.0);NumericVector offsetz(1,0.0);
	
	
	int step = 0;
	int ind = 0; 
	NumericVector unif(1,0.0);
	NumericMatrix samples(size,10, zeroes7.begin());
	
	for(int i = 0; i < L; i++){
		for(int j = 0; j < L; j++){
			L110[0] = std::sqrt(se2[0] + bxseq[i]*bxseq[i]*sxseq[j]);
			L120[0] = bxseq[i]*sxseq[j]/L110[0];
			L220[0] = std::sqrt(sxseq[j]/del[0] - std::pow(sxseq[j]*bxseq[i],2.0)/std::pow(L110[0],2.0));
			V10[0] = s11[0]/L110[0]; V30[0] = s22[0]/L220[0]; V20[0] = (s12[0] - L120[0]*V10[0])/L220[0];
			dL110bx[0] = L120[0]; dL110sx[0] = bxseq[i]*bxseq[i]/(2.0*L110[0]);
			dL120bx[0] = (sxseq[j] - std::pow(L120[0], 2.0))/L110[0]; dL120sx[0] = (bxseq[i]/L110[0])*(1.0 - 0.5*std::pow(L120[0], 2.0)/sxseq[j]);
			dL220bx[0] = -L120[0]*dL120bx[0]/L220[0]; dL220sx[0] = -L120[0]*dL120sx[0]/L220[0];
			dV10bx[0] = (-s11[0]/std::pow(L110[0],2.0)) * dL110bx[0];dV10sx[0] = (-s11[0]/std::pow(L110[0],2.0)) * dL110sx[0];
			dV30bx[0] = (-s22[0]/std::pow(L220[0],2.0)) * dL220bx[0];dV30sx[0] = (-s22[0]/std::pow(L220[0],2.0)) * dL220sx[0];
			dV20bx[0] = (-s12[0]/std::pow(L220[0],2.0)) * dL220bx[0] - (1.0/std::pow(L220[0],2.0))*( (dV10bx[0]*L120[0]+V10[0]*dL120bx[0])*L220[0] - dL220bx[0]*V10[0]*L120[0] );
			dV20sx[0] = (-s12[0]/std::pow(L220[0],2.0)) * dL220sx[0] - (1.0/std::pow(L220[0],2.0))*( (dV10sx[0]*L120[0]+V10[0]*dL120sx[0])*L220[0] - dL220sx[0]*V10[0]*L120[0] );
			c2[0] = (-dV30sx[0]+dV30bx[0]*dV10sx[0]/dV10bx[0])/((-dV20bx[0]/dV10bx[0])+dV20sx[0]);
			c1[0] = (-c2[0]*dV20bx[0]-dV30bx[0])/dV10bx[0];
			eta[0] = c1[0]*V10[0] + c2[0]*V20[0] + V30[0]; 
			denscurr[0] = log(std::abs(1.0/c2[0]))+R::dchisq(sampcurr[0]*sampcurr[0], n[0] - 2.0, 1) + R::dchisq(sampcurr[1]*sampcurr[1], n[0] - 3.0, 1) + R::dnorm((eta[0] -c1[0]*sampcurr[0] - sampcurr[1])/c2[0], 0.0, 1.0, 1);
			ind = 0; step = 0;
			NumericVector densx(1,0.0); NumericVector densz(1,0.0); NumericVector bxs(size,0.0); NumericVector bzs(size,0.0); 
			sampcurr[0] = V10[0]; sampcurr[1] = V30[0];
			while(ind < size){
				if(step > 0){
					densx.push_back(0.0);	densz.push_back(0.0);
				}
				sampprop[0] = R::rnorm(sampcurr[0], propsd[0]); sampprop[1] = R::rnorm(sampcurr[1], propsd[1]);
				if((sampprop[0] > 0.0) & (sampprop[1] > 0.0)){
					densprop[0] = log(std::abs(1.0/c2[0]))+R::dchisq(sampprop[0]*sampprop[0], n[0] - 2.0, 1) + R::dchisq(sampprop[1]*sampprop[1], n[0] - 3.0, 1) + R::dnorm((eta[0] -c1[0]*sampprop[0] - sampprop[1])/c2[0], 0.0, 1.0, 1);
					unif[0] = R::runif(0.0,1.0);
				}else {
					unif[0] = 2.0;	densprop[0] = denscurr[0];
				}
				if(unif[0] < std::exp(densprop[0] - denscurr[0])){
					V1[0] = sampprop[0]; V3[0] = sampprop[1]; V2[0] = (eta[0]-c1[0]*V1[0]-V3[0])/c2[0];
					if(std::pow(s11[0]/V1[0], 2.0)>se2[0]){
						denscurr[0] = densprop[0];
						sampcurr[0] = sampprop[0]; sampcurr[1] = sampprop[1];
						NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
						NumericVector sx(1,0.0);NumericVector bx(1,0.0);NumericVector bz(1,0.0);NumericVector mux(1,0.0);NumericVector se(1,0.0);
						NumericVector t(1,0.0);
						/*t[0] = std::pow(s11[0]/V1[0], 2.0)-se2[0];
						bx[0] = (V1[0]*t[0]/std::sqrt(se2[0] + t[0]))/(s12[0] - s22[0]*V2[0]/V3[0]);
						sx[0] = t[0]/std::pow(bx[0], 2.0);
						L11[0] = std::sqrt(se2[0] + bx[0]*bx[0]*sx[0]);
						L22[0] = std::sqrt(sx[0]/del[0] - std::pow(bx[0]*sx[0]/L11[0], 2.0));
						L12[0] = bx[0]*sx[0]/L11[0];*/
						L11[0] = s11[0]/V1[0];
						sx[0] = std::pow(s22[0]*L11[0]/V3[0], 2.0)/((std::pow(L11[0], 2.0)*(-1.0+1.0/del[0])) +se2[0]);
						bx[0] = std::sqrt((-std::pow(s22[0]/V3[0], 2.0)+sx[0]/del[0])*std::pow(L11[0]/sx[0], 2.0));
						L22[0] = std::sqrt(sx[0]/del[0] - std::pow(bx[0]*sx[0]/L11[0], 2.0));
						if(s12[0] < V2[0]*L22[0]){
							bx[0] = -bx[0];	
						}
						L12[0] = bx[0]*sx[0]/L11[0];
						mux[0] = wbar[0] - L12[0]*Z1[ind]-L22[0]*Z2[ind];
						bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[ind];
						se[0] = (L11[0]*L11[0])-bx[0]*bx[0]*sx[0];
						densx[step] = denscurr[0];
						densz[step] = densx[step] + R::dnorm(Z1[ind], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[ind], 0.0, std::sqrt(1.0/n[0]), 1);
						bxs[ind] = bx[0]; 
						bzs[ind] = bz[0]; 
						samples(ind,0) = bx[0]; samples(ind,1) = bz[0]; samples(ind,2) = mux[0]; samples(ind,3) = sx[0]; samples(ind,4) = se[0]; 
						samples(ind,5) = densx[step];	
						samples(ind,6) = densz[step];
						samples(ind,7) = s11[0] - L11[0]*V1[0];
						samples(ind,8) = s12[0] - L12[0]*V1[0] - L22[0]*V2[0];
						samples(ind,9) = s22[0] - L22[0]*V3[0];
						ind = ind+1;
					}
				}else {
					V1[0] = sampcurr[0]; V3[0] = sampcurr[1]; V2[0] = (eta[0]-c1[0]*V1[0]-V3[0])/c2[0];
					if(std::pow(s11[0]/V1[0], 2.0)>se2[0]){
						NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
						NumericVector sx(1,0.0);NumericVector bx(1,0.0);NumericVector bz(1,0.0);NumericVector mux(1,0.0);NumericVector se(1,0.0);
						NumericVector t(1,0.0);
						/*t[0] = std::pow(s11[0]/V1[0], 2.0)-se2[0];
						bx[0] = (V1[0]*t[0]/std::sqrt(se2[0] + t[0]))/(s12[0] - s22[0]*V2[0]/V3[0]);
						sx[0] = t[0]/std::pow(bx[0], 2.0);
						L11[0] = std::sqrt(se2[0] + bx[0]*bx[0]*sx[0]);
						L22[0] = std::sqrt(sx[0]/del[0] - std::pow(bx[0]*sx[0]/L11[0], 2.0));
						L12[0] = bx[0]*sx[0]/L11[0];*/
						L11[0] = s11[0]/V1[0];
						sx[0] = std::pow(s22[0]*L11[0], 2.0)/((std::pow(L11[0], 2.0)*(-1.0+1.0/del[0])) +se2[0]);
						bx[0] = std::sqrt((-std::pow(s22[0], 2.0)+sx[0]/del[0])*std::pow(L11[0]/sx[0], 2.0));
						L22[0] = std::sqrt(sx[0]/del[0] - std::pow(bx[0]*sx[0]/L11[0], 2.0));
						if(s12[0] < V2[0]*L22[0]){
							bx[0] = -bx[0];	
						}
						L12[0] = bx[0]*sx[0]/L11[0];
						mux[0] = wbar[0] - L12[0]*Z1[ind]-L22[0]*Z2[ind];
						bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[ind];
						se[0] = (L11[0]*L11[0])-bx[0]*bx[0]*sx[0];
						densx[step] = denscurr[0];
						densz[step] = densx[step] + R::dnorm(Z1[ind], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[ind], 0.0, std::sqrt(1.0/n[0]), 1);
						bxs[ind] = bx[0]; 
						bzs[ind] = bz[0]; 
						samples(ind,0) = bx[0]; samples(ind,1) = bz[0]; samples(ind,2) = mux[0]; samples(ind,3) = sx[0]; samples(ind,4) = se[0]; 
						samples(ind,5) = densx[step];	
						samples(ind,6) = densz[step];
						samples(ind,7) = s11[0] - L11[0]*V1[0];
						samples(ind,8) = s12[0] - L12[0]*V1[0] - L22[0]*V2[0];
						samples(ind,9) = s22[0] - L22[0]*V3[0];
						ind = ind+1;
					}	
				}
				step = step + 1;
			}
			
			// Compute plausibility
			
			NumericVector plausesx(L,0.0);
			NumericVector plausesz(pL,0.0);
			std::sort(bxs.begin(), bxs.end());std::sort(bzs.begin(), bzs.end());
			NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
			
			std::sort(densx.begin(), densx.end()); std::sort(densz.begin(), densz.end());
			
			
			samples = sortmat(samples,5);
			offsetx[0] = densx[step - 1] - samples(ind-1,5);
			int unifind =0; int ind2 =0;
			bool comp = true;
			for(int k=0; k<size; k++){
				randsetslo[0] = bxs[(size-1)]; randsetshi[0] = bxs[0];
				ind2 = 0;
				unifind = round(R::runif(0.0,1.0)*(step-1));
				comp = (samples(ind2,5) < (densx[unifind] - offsetx[0]));
				while(comp & (ind2 < (ind - 1))){
					ind2 = ind2 + 1	;
					comp = (samples(ind2,5) < (densx[unifind] - offsetx[0]));
				}
				NumericVector subset(size - ind2, 0.0);
				for(int l = 0; l < (size - ind2); l++){
					subset[l] = samples(ind2+l,0);	
				}
				randsetslo[0] = Rcpp::min(subset);  randsetshi[0] = Rcpp::max(subset); 
				if(   (bxseq[i] >= randsetslo[0]) & (bxseq[i] <= randsetshi[0])   ){
					plausesx[i] = plausesx[i]+(1.0/(size));
				}
			}
			
			maxplausesx[i] = std::max(maxplausesx[i], plausesx[i]);
						

			
			
			samples = sortmat(samples,6);
			offsetz[0] = densz[step - 1] - samples(ind-1,6);
			unifind =0; ind2 =0;
			comp = true;
			for(int k=0; k<size; k++){
				randsetslo[0] = bzs[(size-1)]; randsetshi[0] = bzs[0];
				ind2 = 0;
				unifind = round(R::runif(0.0,1.0)*(step-1));
				comp = (samples(ind2,6) < (densz[unifind] - offsetz[0]));
				while(comp & (ind2 < (ind - 1))){
					ind2 = ind2 + 1	;
					comp = (samples(ind2,6) < (densz[unifind] - offsetz[0]));
				}
				NumericVector subset(size - ind2, 0.0);
				for(int l = 0; l < (size - ind2); l++){
					subset[l] = samples(ind2+l,1);	
				}
				randsetslo[0] = Rcpp::min(subset);  randsetshi[0] = Rcpp::max(subset); 
				for(int l=0; l<pL; l++){
					if(   (plbzseq[l] >= randsetslo[0]) & (plbzseq[l] <= randsetshi[0])   ){
						plausesz[l] = plausesz[l]+(1.0/(size));
					}
				}
			}
			
			for(int l=0; l<pL; l++){
				maxplausesz[l] = std::max(maxplausesz[l], plausesz[l]);
			}
				
		}
	}
			
	result = Rcpp::List::create(Rcpp::Named("plauses_beta_x") = maxplausesx,  Rcpp::Named("plauses_beta_z") = maxplausesz, Rcpp::Named("samples") = samples);		

	return result;
	
}

Rcpp::List plauscontourSIR(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector n, NumericVector mode, NumericVector dens, NumericVector se2, NumericVector cond_par) {

	List result;
	Rcpp::Function sortmat("sortmat");
	int size = round(sampsize[0]);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
		
	// Generate MC sample of aux rvs Z1, Z2, V2

	NumericVector V2(size,0.0); NumericVector V3(size,0.0); NumericVector Z1(size,0.0); NumericVector Z2(size,0.0);
	NumericVector sV2(size,0.0); NumericVector sV3(size,0.0); NumericVector sZ1(size,0.0); NumericVector sZ2(size,0.0);
	NumericVector sumweights(1,0.0); IntegerVector indices(size); NumericVector nacheck(1,0.0); LogicalVector na(1);
	
	
	for(int k = 0; k < size; k++){
		V2[k] = R::rnorm(mode[0], 1.0);V3[k] = R::rchisq(mode[1]);Z1[k] = R::rnorm(0.0, std::sqrt(1.0/n[0]));Z2[k] = R::rnorm(0.0, std::sqrt(1.0/n[0]));
	}
	V2[0] = mode[0];V3[0] = mode[1];
	
	
	NumericVector zeroes7(size*10, 0.0);
	NumericVector offset(1,0.0);
	NumericMatrix samples(size,11, zeroes7.begin());
	

	

	sumweights[0] = 0.0;
	NumericVector weights(size,0.0);
	for(int k = 0; k < size; k++){
		if(V3[k]>0){
			V3[k] = std::sqrt(V3[k]);
			if(((cond_par[0]-V3[k]-cond_par[2]*V2[k])/cond_par[1]) > 0.0){
				weights[k] = log(std::abs(1.0/cond_par[1]))+R::dnorm(V2[k], 0.0, 1.0, 1)-R::dnorm(V2[k], mode[0], 1.0, 1) + R::dchisq(V3[k]*V3[k], n[0]-2.0, 1) - R::dchisq(V3[k]*V3[k], mode[1], 1) + R::dchisq(std::pow((cond_par[0]-V3[k]-cond_par[2]*V2[k])/cond_par[1],2.0), n[0]-2.0, 1);
				nacheck[0]=weights[k];
				na = Rcpp::is_na(nacheck);
				if(na[0]){
					weights[k] = 0.0;	
				}else {
					weights[k] = std::exp(weights[k]);	
				}
			}else {
				weights[k] = 0.0;	
			}
		}
		sumweights[0] = sumweights[0] + weights[k];
	}
	if(sumweights[0]>0.0){
		for(int k = 0; k < size; k++){
			weights[k] = weights[k]/sumweights[0];
		}
		indices = Rcpp::sample(size, size, true, weights, true);
	}
	for(int k = 0; k < size; k++){
		samples(k,0) = V2[k];samples(k,1) = V3[k];samples(k,2) = (cond_par[0]-V3[k]-cond_par[2]*V2[k])/cond_par[1];
		samples(k,4) = log(std::abs(1.0/cond_par[1]))+R::dnorm(V2[k], 0.0, 1.0, 1) + R::dchisq(V3[k]*V3[k], n[0]-2.0, 1) + R::dchisq(std::pow((cond_par[0]-V3[k]-cond_par[2]*V2[k])/cond_par[1],2.0), n[0]-1.0, 1);	
		sV2[k] = V2[indices[k]-1]; sV3[k] = V3[indices[k]-1]; sZ1[k] = Z1[indices[k]-1]; sZ2[k] = Z2[indices[k]-1]; 
		samples(k,5) = log(std::abs(1.0/cond_par[1]))+R::dnorm(sV2[k], 0.0, 1.0, 1) + R::dchisq(sV3[k]*sV3[k], n[0]-2.0, 1) + R::dchisq(std::pow((cond_par[0]-sV3[k]-cond_par[2]*sV2[k])/cond_par[1],2.0), n[0]-1.0, 1);	
		samples(k,6) = samples(k,5) + R::dnorm(sZ1[k], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(sZ2[k], 0.0, std::sqrt(1.0/n[0]), 1);
		samples(k,7) = sV2[k];
		samples(k,8) = (cond_par[0]-sV3[k]-cond_par[2]*sV2[k])/cond_par[1];
		samples(k,9) = sV3[k];
		samples(k,10) = weights[indices[k]-1];
	}	
	
	// Compute plausibility
	NumericVector plaus(1,0.0);
	if(sumweights[0]>0.0){
		
		
		NumericVector randsetdens(1,0.0);
		samples = sortmat(samples,5);

		
		for(int j=0; j<size; j++){
			randsetdens[0] = samples(j,5);
			if(   dens[0]>(randsetdens[0]-offset[0]) ){
				plaus[0] = plaus[0]+(1.0/(size));
			}
		}
		
	}

			
	result = Rcpp::List::create(Rcpp::Named("plaus") = plaus, Rcpp::Named("samples") = samples);		

	return result;
	
}


/*
Rcpp::List plauscontourSIR(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector n, NumericVector mode, NumericVector dens, NumericVector se2, NumericVector cond_par) {

	List result;
	Rcpp::Function sortmat("sortmat");
	int size = round(sampsize[0]);
	
	NumericVector s11(1,0.0); s11[0] = stat[0];
	NumericVector s12(1,0.0); s12[0] = stat[1];	
	NumericVector s22(1,0.0); s22[0] = stat[2];
	NumericVector ybar(1,0.0); ybar[0] = stat[3];
	NumericVector wbar(1,0.0); wbar[0] = stat[4];
	
		
	// Generate MC sample of aux rvs Z1, Z2, V2

	NumericVector V1(size,0.0); NumericVector V2(size,0.0); NumericVector V3(size,0.0); NumericVector Z1(size,0.0); NumericVector Z2(size,0.0);
	NumericVector sV1(size,0.0); NumericVector sV3(size,0.0); NumericVector sZ1(size,0.0); NumericVector sZ2(size,0.0);
	NumericVector sumweights(1,0.0); IntegerVector indices(size); 
	
	
	for(int k = 0; k < size; k++){
		V1[k] = R::rchisq(mode[0]);V3[k] = R::rchisq(mode[1]);Z1[k] = R::rnorm(0.0, std::sqrt(1.0/n[0]));Z2[k] = R::rnorm(0.0, std::sqrt(1.0/n[0]));
		V1[k] = std::sqrt(V1[k]);V3[k] = std::sqrt(V3[k]);
	}
	
	
	
	NumericVector zeroes7(size*10, 0.0);
	NumericVector offset(1,0.0);
	NumericMatrix samples(size,11, zeroes7.begin());
	

	

	sumweights[0] = 0.0;
	NumericVector weights(size,0.0);
	for(int k = 0; k < size; k++){
		weights[k] = std::exp(log(std::abs(1.0/cond_par[2]))+R::dnorm((cond_par[0] -cond_par[1]*V1[k] - V3[k])/cond_par[2], 0.0, 1.0, 1) + R::dchisq(V1[k]*V1[k], n[0]-1.0, 1) + R::dchisq(V3[k]*V3[k], n[0]-2.0, 1) - R::dchisq(V1[k]*V1[k], mode[0], 1) - R::dchisq(V3[k]*V3[k], mode[1], 1));
		sumweights[0] = sumweights[0] + weights[k];
	}
	if(sumweights[0]>0){
	for(int k = 0; k < size; k++){
		weights[k] = weights[k]/sumweights[0];
	}
	indices = Rcpp::sample(size, size, true, weights, true);
	for(int k = 0; k < size; k++){
		sV1[k] = V1[indices[k]-1]; sV3[k] = V3[indices[k]-1]; sZ1[k] = Z1[indices[k]-1]; sZ2[k] = Z2[indices[k]-1]; 
		samples(k,5) = log(std::abs(1.0/cond_par[2]))+R::dnorm((cond_par[0] -cond_par[1]*sV1[k] - sV3[k])/cond_par[2], 0.0, 1.0, 1) + R::dchisq(sV1[k]*sV1[k], n[0]-1.0, 1) + R::dchisq(sV3[k]*sV3[k], n[0]-2.0, 1);	
		samples(k,6) = samples(k,5) + R::dnorm(sZ1[k], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(sZ2[k], 0.0, std::sqrt(1.0/n[0]), 1);
		samples(k,7) = sV1[k];
		samples(k,8) = (cond_par[0] -cond_par[1]*sV1[k] - sV3[k])/cond_par[2];
		samples(k,9) = sV3[k];
		samples(k,10) = weights[indices[k]-1];
	}	
	}



	// Compute plausibility

	NumericVector plaus(1,0.0);
	NumericVector randsetdens(1,0.0);
	samples = sortmat(samples,5);

	

		
		for(int j=0; j<size; j++){
			randsetdens[0] = samples(j,5);
			if(   dens[0]>(randsetdens[0]-offset[0]) ){
				plaus[0] = plaus[0]+(1.0/(size));
			}
		}
		


			
	result = Rcpp::List::create(Rcpp::Named("plaus") = plaus, Rcpp::Named("samples") = samples);		

	return result;
	
}

*/


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

Rcpp::NumericVector grow(NumericVector x) {
  arma::vec xx = as<arma::vec>(x);
  xx.resize( xx.size()+1, 0.0 );
  NumericMatrix	y = wrap(xx); 
  return y;
}

Rcpp::List plauscontourMC2(NumericVector sampsize, NumericVector stat, NumericVector del, NumericVector type, NumericVector n, NumericVector truebx, NumericVector bxseq,NumericVector truebz, NumericVector bzseq, NumericVector randsettype) {

	List result;
	Rcpp::Function sortmat("sortmat");
	Rcpp::Function grow("grow");
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
	NumericVector zeroes(2*size,0.0); ;
	NumericMatrix bxs(size,2,zeroes.begin()); NumericMatrix bzs(size,2,zeroes.begin());
	NumericVector bx(1,0.0); NumericVector sx(1,0.0); NumericVector se(1,0.0); NumericVector mux(1,0.0); NumericVector bz(1,0.0); 
	NumericVector dens_samps_x(1,0.0); NumericVector dens_samps_z(1,0.0);
	NumericVector bx_s(size, 0.0); NumericVector bz_s(size, 0.0);
	NumericVector sxs(size,0.0); 
	int ind = 0; int step = 0;
	while(ind < size){
		if(step > 0){
			dens_samps_x.push_back( 0.0 );
			dens_samps_z.push_back( 0.0 );	
		}
		V2[0] = R::rnorm( 0.0, 1.0 );V1[0] = R::rchisq( n[0]-1 );V3[0] = R::rchisq( n[0]-2 );Z1[0] = R::rnorm( 0.0, std::sqrt(1.0/n[0]) );Z2[0] = R::rnorm( 0.0, std::sqrt(1.0/n[0]) );
		V1[0] = std::sqrt(V1[0]);	
		V3[0] = std::sqrt(V3[0]);
		dens_samps_x[step] = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);
		dens_samps_z[step] = dens_samps_x[ind] + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
		if(type[0] == 2.0){
			NumericVector L11(1,0.0);NumericVector L12(1,0.0);NumericVector L22(1,0.0);
			L11[0] = s11[0]/V1[0]; L22[0] = s22[0]/V3[0]; L12[0] = (s12[0] - V2[0]*L22[0])/V1[0]; 
			sx[0] = 0.5*(-(L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0]) + std::sqrt(((L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0])*(L11[0]*L11[0]/del[0] - L22[0]*L22[0] - L12[0]*L12[0]))+4*L11[0]*L11[0]*L12[0]*L12[0]/del[0]));
			bx[0] = L11[0]*L12[0]/sx[0];
			se[0] = L11[0]*L11[0]-sx[0]*bx[0]*bx[0];
			mux[0] = wbar[0] - L12[0]*Z1[0]-L22[0]*Z2[0];
			bz[0] = ybar[0] - bx[0]*mux[0]-L11[0]*Z1[0];
			if((se[0] > 0.0) & (sx[0]>0.0)){
				bxs(ind,0) = bx[0]; bx_s[ind] = bx[0]; 
				bzs(ind,0) = bz[0]; bz_s[ind] = bz[0]; 
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
				bxs(ind,0) = bx[0]; bx_s[ind] = bx[0]; 
				bzs(ind,0) = bz[0]; bz_s[ind] = bz[0]; 
				bxs(ind,1) = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);	
				bzs(ind,1) = bxs(ind,1) + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
				sxs[ind] = sx[0];
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
				bxs(ind,0) = bx[0]; bx_s[ind] = bx[0]; 
				bzs(ind,0) = bz[0]; bz_s[ind] = bz[0]; 
				bxs(ind,1) = R::dchisq(V1[0]*V1[0], n[0]-1, 1) + R::dchisq(V3[0]*V3[0], n[0]-2, 1)  + R::dnorm(V2[0], 0.0, 1.0, 1);	
				bzs(ind,1) = bxs(ind,1) + R::dnorm(Z1[0], 0.0, std::sqrt(1.0/n[0]), 1) + R::dnorm(Z2[0], 0.0, std::sqrt(1.0/n[0]), 1);
				ind = ind+1;
			}
		}
		step = step + 1;
	}

	bxs = sortmat(bxs,1);
	bzs = sortmat(bzs,1);
	std::sort(bx_s.begin(), bx_s.end()); std::sort(bz_s.begin(), bz_s.end());
	std::sort(dens_samps_x.begin(), dens_samps_x.end()); std::sort(dens_samps_z.begin(), dens_samps_z.end());
	
	// plausibility
	NumericVector plausesx(500,0.0);
	NumericVector plausestrux(1,0.0);
	NumericVector plausesz(500,0.0);
	NumericVector plausestruz(1,0.0);
	
	NumericVector inc(1, 0.0); inc[0] = (bx_s[size-1]-bx_s[0])/500.0;
	NumericVector inc2(1, 0.0); inc2[0] = (bz_s[size-1]-bz_s[0])/500.0;
	NumericVector bx_seq(500, 0.0);NumericVector bz_seq(500, 0.0);
	bx_seq[0] = bx_s[0];bz_seq[0] = bz_s[0];
	for(int i = 1; i < 500; i++){
		bx_seq[i] = bx_seq[i-1]+inc[0]; bz_seq[i] = bz_seq[i-1]+inc2[0];	
	}
	bx_seq[499] = bx_s[(size-1)];bz_seq[499] = bz_s[(size-1)];	
	
	NumericVector offset(1,0.0);
	offset[0] = dens_samps_x[step - 1] - bxs(ind-1,1);
	int ind2 = 0;
	
	NumericVector randsetslo(1,0.0);NumericVector randsetshi(1,0.0);
	int unifind =0;
	bool comp = true;
	for(int j=0; j<size; j++){
		randsetslo[0] = bx_s[(size-1)]; randsetshi[0] = bx_s[0];
		ind2 = 0;
		unifind = round(R::runif(0.0,1.0)*(step-1));
		comp = (bxs(ind2,1) < (dens_samps_x[unifind] - offset[0]));
		while(comp & (ind2 < (ind - 1))){
			ind2 = ind2 + 1	;
			comp = (bxs(ind2,1) < (dens_samps_x[unifind] - offset[0]));
		}
		NumericVector subset(size - ind2, 0.0);
		for(int i = 0; i < (size - ind2); i++){
			subset[i] = bxs(ind2+i,0);	
		}
		randsetslo[0] = Rcpp::min(subset);  randsetshi[0] = Rcpp::max(subset); 
		/*while((bxs(ind-1-ind2,1) >= (dens_samps_x[unifind] - offset[0])) & (ind2 < ind) ){
			randsetslo[0] = std::min(randsetslo[0], bxs(ind-1-ind2,0));
			randsetshi[0] = std::max(randsetshi[0], bxs(ind-1-ind2,0));
			ind2 = ind2+1;
		}*/
		if(   (truebx[0] > randsetslo[0]) & (truebx[0] < randsetshi[0])   ){
			plausestrux[0] = plausestrux[0]+(1.0/size);
		}
		for(int i=0; i<500; i++){
			if(   (bx_seq[i] >= randsetslo[0]) & (bx_seq[i] <= randsetshi[0])   ){
				plausesx[i] = plausesx[i]+(1.0/size);
			}
		}
	}
	
	ind2 = 0;
	offset[0] = dens_samps_z[(step - 1)] - bzs((ind-1),1);
	
	unifind =0;
	for(int j=0; j<size; j++){
		randsetslo[0] = bz_s[(size-1)]; randsetshi[0] = bz_s[0];
		ind2 = 0;
		unifind = round(R::runif(0.0,1.0)*(step-1));
		comp = (bzs(ind2,1) < (dens_samps_z[unifind] - offset[0]));
		while(comp & (ind2 < (ind - 1))){
			ind2 = ind2 + 1	;
			comp = (bzs(ind2,1) < (dens_samps_z[unifind] - offset[0]));
		}
		NumericVector subset(size - ind2, 0.0);
		for(int i = 0; i < (size - ind2); i++){
			subset[i] = bzs(ind2+i,0);	
		}
		randsetslo[0] = Rcpp::min(subset);  randsetshi[0] = Rcpp::max(subset); 
		/*while((bzs(ind-1-ind2,1) >= (dens_samps_z[unifind] - offset[0]) ) & (ind2 < ind)  ){
			randsetslo[0] = std::min(randsetslo[0], bzs(ind-1-ind2,0));
			randsetshi[0] = std::max(randsetshi[0], bzs(ind-1-ind2,0));
			ind2 = ind2+1;
		}*/
		if(   (truebz[0] > randsetslo[0]) & (truebz[0] < randsetshi[0])   ){
			plausestruz[0] = plausestruz[0]+(1.0/size);
		}
		for(int i=0; i<500; i++){
			if(   (bz_seq[i] >= randsetslo[0]) & (bz_seq[i] <= randsetshi[0])   ){
				plausesz[i] = plausesz[i]+(1.0/size);
			}
		}
	}
		
	result = Rcpp::List::create(Rcpp::Named("sxs") = sxs, Rcpp::Named("plaus_beta_x") = plausestrux, Rcpp::Named("plauses_beta_x") = plausesx,  Rcpp::Named("samples_bx") = bxs, Rcpp::Named("plaus_beta_z") = plausestruz, Rcpp::Named("plauses_beta_z") = plausesz,  Rcpp::Named("samples_bz") = bzs, Rcpp::Named("bx_seq") = bx_seq, Rcpp::Named("bz_seq") = bz_seq, Rcpp::Named("step") = step, Rcpp::Named("dens_samps_x") = dens_samps_x, Rcpp::Named("dens_samps_z") = dens_samps_z);		

	return result;
	
}





Rcpp::NumericVector loglik(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n) {

	int nn = round(n[0]);
	
	NumericVector LL(1, 0.0);
	
	NumericVector beta_x(1, 0.0); beta_x[0] = theta[0];
	NumericVector s_x2(1, 0.0); s_x2[0] = theta[1];
	NumericVector s_e2(1, 0.0); s_e2[0] = theta[2];

	NumericVector s11(1, 0.0); s11[0] = stat[0];
	NumericVector s12(1, 0.0); s12[0] = stat[1];
	NumericVector s22(1, 0.0); s22[0] = stat[2];
	
	NumericVector L11(1, 0.0); L11[0] = std::sqrt(beta_x[0]*beta_x[0]*s_x2[0]+s_e2[0]);
	NumericVector L12(1, 0.0); L12[0] = beta_x[0]*s_x2[0]/L11[0];
	NumericVector L22(1, 0.0); L22[0] = s_x2[0]/del[0] - (L12[0]*L12[0]);
	if(L22[0] > 0.0){
		L22[0] = std::sqrt(L22[0]);	
		LL[0] = R::dchisq((s11[0]/L11[0])*(s11[0]/L11[0]), nn-1, 1) + R::dchisq((s22[0]/L22[0])*(s22[0]/L22[0]), nn-2, 1) + R::dnorm((s12[0] - L12[0]*s11[0]/L11[0])/L22[0], 0.0, 1.0, 1);
	}else {
		LL[0] = R_NegInf;
	}
	
	return LL;
	
}


double negloglik(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n) {

	int nn = round(n[0]);
	
	NumericVector LL(1, 0.0);
	
	NumericVector beta_x(1, 0.0); beta_x[0] = theta[0];
	NumericVector s_x2(1, 0.0); s_x2[0] = theta[1];
	NumericVector s_e2(1, 0.0); s_e2[0] = theta[2];

	NumericVector s11(1, 0.0); s11[0] = stat[0];
	NumericVector s12(1, 0.0); s12[0] = stat[1];
	NumericVector s22(1, 0.0); s22[0] = stat[2];
	
	NumericVector L11(1, 0.0); L11[0] = std::sqrt(beta_x[0]*beta_x[0]*s_x2[0]+s_e2[0]);
	NumericVector L12(1, 0.0); L12[0] = beta_x[0]*s_x2[0]/L11[0];
	NumericVector L22(1, 0.0); L22[0] = s_x2[0]/del[0] - (L12[0]*L12[0]);
	if(L22[0] > 0.0){
		L22[0] = std::sqrt(L22[0]);	
		LL[0] = -1.0*(R::dchisq((s11[0]/L11[0])*(s11[0]/L11[0]), nn-1, 1) + R::dchisq((s22[0]/L22[0])*(s22[0]/L22[0]), nn-2, 1) + R::dnorm((s12[0] - L12[0]*s11[0]/L11[0])/L22[0], 0.0, 1.0, 1));
	}else {
		LL[0] = 1000000;
	}
	
	return LL[0];
	
}
/*
Rcpp::NumericVector grloglik(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n) {

	Rcpp::Function loglik("loglik");
	NumericVector thetastar(3,0.0);NumericVector grad(3,0.0);
	thetastar[0] = theta[0]; thetastar[1] = theta[1]; thetastar[2] = theta[2];
	
	thetastar[0] = thetastar[0] + 0.01;
	grad[0] = loglik(thetastar, stat, del, n);
	thetastar[0] = thetastar[0] - 0.02;
	grad[0] = (grad[0] - loglik(thetastar, stat, del, n))/0.02;
	thetastar[0] = thetastar[0] + 0.01;
	
	thetastar[1] = thetastar[1] + 0.01;
	grad[1] = loglik(thetastar, stat, del, n);
	thetastar[1] = thetastar[1] - 0.02;
	grad[1] = (grad[1] - loglik(thetastar, stat, del, n))/0.02;
	thetastar[1] = thetastar[1] + 0.01;
	
	thetastar[2] = thetastar[2] + 0.01;
	grad[2] = loglik(thetastar, stat, del, n);
	thetastar[2] = thetastar[2] - 0.02;
	grad[2] = (grad[2] - loglik(thetastar, stat, del, n))/0.02;
	thetastar[2] = thetastar[2] + 0.01;
	
	return grad;	
}
*/

List optimrcpp(NumericVector theta, NumericVector stat, NumericVector del, NumericVector n){

  // Extract R's optim function
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];

  NumericVector lower(3,0.0); lower[0] = -1000; lower[1] = 0.01; lower[2] = 0.01;

  // Call the optim function from R in C++ 
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = theta,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&negloglik),
                                 Rcpp::_["method"] = "L-BFGS-B",
				 Rcpp::_["lower"] = lower,
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["stat"] = stat,
                                 Rcpp::_["del"] = del,
                                 Rcpp::_["n"] = n);

  // Return results
  return opt_results;
}

Rcpp::NumericVector maxloglik(NumericMatrix thetas, NumericVector stat, NumericVector del, NumericVector n) {
	
	int N = round(thetas.nrow());
	
	int nn = round(n[0]);
	
	NumericVector beta_x(1, 0.0); 
	NumericVector s_x2(1, 0.0); 
	NumericVector s_e2(1, 0.0); 

	NumericVector s11(1, 0.0); s11[0] = stat[0];
	NumericVector s12(1, 0.0); s12[0] = stat[1];
	NumericVector s22(1, 0.0); s22[0] = stat[2];
	
	NumericVector L11(1, 0.0); 
	NumericVector L12(1, 0.0); 
	NumericVector L22(1, 0.0);
	
	NumericVector LL(1, R_NegInf);
	
	for(int i = 0; i < N; i++){
	
		beta_x[0] = thetas(i,0);
		s_x2[0] = thetas(i,1);
		s_e2[0] = thetas(i,2);
		
		L11[0] = std::sqrt(beta_x[0]*beta_x[0]*s_x2[0]+s_e2[0]);
		L12[0] = beta_x[0]*s_x2[0]/L11[0];
		L22[0] = s_x2[0]/del[0] - (L12[0]*L12[0]);
		
		if(L22[0] > 0.0){
			L22[0] = std::sqrt(L22[0]);	
			LL[0] = std::max(LL[0], R::dchisq((s11[0]/L11[0])*(s11[0]/L11[0]), nn-1, 1) + R::dchisq((s22[0]/L22[0])*(s22[0]/L22[0]), nn-2, 1) + R::dnorm((s12[0] - L12[0]*s11[0]/L11[0])/L22[0], 0.0, 1.0, 1));
		}
				
	}

	return LL;
	
}
	
	
Rcpp::List genIMplaus(NumericMatrix thetas, NumericVector stat, NumericVector del, NumericVector n, NumericVector M) {
	
	Rcpp::Function loglik("loglik");
	Rcpp::Function maxloglik("maxloglik");
	
	int nn = round(n[0]);
	int N = round(thetas.nrow());
	int m = round(M[0]);
	
	
	
	NumericVector theta(3, 0.0);
	NumericVector hdata(N,0.0);
	NumericVector num(1,0.0);
	NumericVector denom(1,0.0);
	LogicalVector check(1);
	denom = maxloglik(thetas, stat, del, n);
	for(int i = 0; i < N; i++){
		theta[0] = thetas(i,0); theta[1] = thetas(i,1); theta[2] = thetas(i,2);
		num = loglik(theta, stat, del, n);
		check = is_finite(num);
		if(check[0]){
			hdata[i] = std::exp(num[0] - denom[0]);
		}
	}
	
	
	NumericVector hsims(m, 0.0);
	NumericVector V1(1, 0.0);
	NumericVector V2(1, 0.0);
	NumericVector V3(1, 0.0);
	NumericVector L11(1, 0.0);
	NumericVector L12(1, 0.0);
	NumericVector L22(1, 0.0);
	NumericVector statMC(3, 0.0);
	NumericVector plauses(N,0.0);
		
	for(int i = 0; i < N; i++){
		theta[0] = thetas(i,0); theta[1] = thetas(i,1); theta[2] = thetas(i,2);
		L11[0] = std::sqrt(theta[0]*theta[0]*theta[1]+theta[2]);
		L12[0] = theta[0]*theta[1]/L11[0];
		L22[0] = std::sqrt((theta[1]/del[0]) - (L12[0]*L12[0]));
		for(int j = 0; j < m; j++){
			V1[0] = std::sqrt(R::rchisq(nn-1));
			V3[0] = std::sqrt(R::rchisq(nn-2));
			V2[0] = R::rnorm(0.0, 1.0);
			statMC[0] = L11[0]*V1[0];
			statMC[1] = L12[0]*V1[0]+L22[0]*V2[0];
			statMC[2] = L22[0]*V3[0];
			denom = maxloglik(thetas, statMC, del, n);
			num = loglik(theta, statMC, del, n);
			check = is_finite(num);
			if(check[0]){
				hsims[j] = std::exp(num[0] - denom[0]);
			}else {
				hsims[j] = 0.0;
			}
			if(hsims[j]<=hdata[i]){
				plauses[i] = plauses[i] + (1.0/M[0]);
			}
		}
		
	}
	
	Rcpp::List result;
	
	result = Rcpp::List::create(Rcpp::Named("plauses") = plauses,Rcpp::Named("hsims") = hsims,Rcpp::Named("hdata") = hdata, Rcpp::Named("V1") = V1, Rcpp::Named("V2") = V2, Rcpp::Named("V3") = V3, Rcpp::Named("statMC") = statMC, Rcpp::Named("num") = num, Rcpp::Named("denom") = denom);		

	return result;
	
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

