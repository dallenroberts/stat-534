/*
 This program uses GSL to create a sequence of 1000
 draws from Uniform(0,1)
*/

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrices.h"

// Function that inputs a pointer to a pxp GSL matrix covX and a pointer to a 
// nxp data matrix X and calculates the pxp sample covariance matrix sigma
// Outputs nothing, but updates the gsl_matrix stored at location covX
void makeCovariance(gsl_matrix* covX, gsl_matrix* X) {

	printf("\n Inside makeCovariance");
	return(1);

}


int main() {

	int i;
  	int n = 158;
  	int p = 51;

	// Initialize random number generator
  	const gsl_rng_type* T;
  	gsl_rng* r;
	
  	gsl_rng_env_setup();
	
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);

  	// Load erdata
  	gsl_matrix * m = gsl_matrix_alloc(n, p);
	FILE * f = fopen("erdata.txt", "r");
	gsl_matrix_fscanf(f, m);
	fclose(f);

	printmatrix("datamat.mat",m);
	

  	// for(i=0;i<n;i++)
  	// {
  	//   double u = gsl_rng_uniform(r);
  	//   printf("%.5lf\n",u);
  	// }

	// Free memory
	gsl_matrix_free(m);
  	gsl_rng_free(r);
	
  	return(1);
}
