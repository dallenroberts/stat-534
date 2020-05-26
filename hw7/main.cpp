/*
 This program uses GSL to create a sequence of 1000
 draws from Uniform(0,1)
*/

#include <stdio.h>
#include <gsl/gsl_rng.h>

// Function that inputs a pointer to a pxp GSL matrix covX and a pointer to a 
// nxp data matrix X and calculates the pxp sample covariance matrix sigma
// Outputs nothing, but updates the gsl_matrix stored at location covX
void makeCovariance(gsl_matrix* covX, gsl_matrix* X) {

	return(1)

}


int main() {

	int i;
  	int n = 1000;

	// Initialize random number generator
  	const gsl_rng_type* T;
  	gsl_rng* r;
	
  	gsl_rng_env_setup();
	
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);
	
	
  	// for(i=0;i<n;i++)
  	// {
  	//   double u = gsl_rng_uniform(r);
  	//   printf("%.5lf\n",u);
  	// }
	
  	gsl_rng_free(r);
	
  	return(1);
}
