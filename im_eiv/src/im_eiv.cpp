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
	
	

