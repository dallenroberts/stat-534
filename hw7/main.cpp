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

// Inputs gsl_matrix* K and outputs Cholesky decomposition as gsl_matrix *
gsl_matrix* makeCholesky(gsl_matrix* K) {

	gsl_matrix* copyK = gsl_matrix_alloc(K->size1,K->size1);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(copyK,K)) {
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}

	if(GSL_SUCCESS!=gsl_linalg_cholesky_decomp(copyK); {
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}

	return(copyK);
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

	// Problem 1
	//Calculate covariance matrix
	gsl_matrix* covX = gsl_matrix_alloc(p, p);
	makeCovariance(covX,X);

	printmatrix("covmat.mat",X);
	
	// Problem 2
	// Calculate Cholesky decomposition
	gsl_matrix* C = makeCholesky(covX);
	printmatrix("cholesky.mat",C);

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
