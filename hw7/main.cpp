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

	printf("\n Inside makeCovariance \n");

	int i,j;
	double cov;
	gsl_vector_view a, b;

	// Calculate covariance for each combination of columns
	for(i=0;i<X->size2;i++) {

		a = gsl_matrix_column(X, i);

		for(j=i;j<X-> size2;j++) {
          	b = gsl_matrix_column(X, j);
          	cov = gsl_stats_covariance(a.vector.data, a.vector.stride, b.vector.data, b.vector.stride, a.vector.size);
          	gsl_matrix_set(covX, i, j, cov);
		}
	}

	return;

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
  	gsl_matrix * X = gsl_matrix_alloc(n, p);
	FILE * f = fopen("erdata.txt", "r");
	gsl_matrix_fscanf(f, X);
	fclose(f);

	printmatrix("datamat.mat",X);

	// Calculate covariance matrix
	gsl_matrix * covX = gsl_matrix_alloc(p, p);
	makeCovariance(covX,X);

	printmatrix("covmat.mat",X);
	

  	// for(i=0;i<n;i++)
  	// {
  	//   double u = gsl_rng_uniform(r);
  	//   printf("%.5lf\n",u);
  	// }

	// Free memory
	gsl_matrix_free(X);
	gsl_matrix_free(covX);
  	gsl_rng_free(r);
	
  	return(1);
}
