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
          	gsl_matrix_set(covX, j, i, cov);
		}
	}

	return;

}

// Inputs gsl_matrix* K and outputs Cholesky decomposition as gsl_matrix *
gsl_matrix* makeCholesky(gsl_matrix* K) {

	gsl_matrix* cholesky = gsl_matrix_alloc(K->size1,K->size2);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(cholesky,K)) {
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}

	if(GSL_SUCCESS!=gsl_linalg_cholesky_decomp(cholesky)) {
		printf("GSL failed cholesky decomposition.\n");
		exit(1);
	}

	return(cholesky);
}

// Samples from the multivariate normal distribution using the Cholesky decomposition
// Inputs random number state mystream, gsl_matrix * samples to hold samples, and 
// gsl_matrix * sigma is the covariance matrix of the multivariate normal
void randomMVN(gsl_rng* mystream, gsl_matrix* samples,gsl_matrix* sigma) {

	printf("\n Inside randomMVN \n");

	int i,j;
	gsl_matrix* z = gsl_matrix_alloc(sigma->size2, 1);

	// Calculate Cholesky decomposition
	gsl_matrix* psi = makeCholesky(sigma);

	// Draw samples
	gsl_matrix* X = gsl_matrix_alloc(sigma->size2, 1); // px1 matrix output by dgemm
	gsl_vector* s = gsl_vector_alloc(sigma->size2); // vector to hold samples

	for(i = 0; i < samples->size1; i++) {

		// Generate p independent N(0,1) random numbers
		for(j = 0; j < sigma->size2; j++) {

			gsl_matrix_set(z, i, 0, gsl_ran_ugaussian(mystream));
		}

		// Calculate matrix product X = psi*Z. Note that this is a px1 matrix
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, psi, z, 0.0, X);

		// Store in samples matrix
		gsl_matrix_get_col(s, X, 0);
		gsl_matrix_set_row(samples, i, s);

	}

	// Free memory
	gsl_matrix_free(psi);
	gsl_matrix_free(z);
	gsl_matrix_free(X);
	gsl_vector_free(s);

}

int main() {

	int i;
  	int n = 158;
  	int p = 51;
  	int nsamples = 100;

	// Initialize random number generator
  	const gsl_rng_type* T;
  	gsl_rng* r;
	
  	gsl_rng_env_setup();
	
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);

  	// Load erdata.txt
  	gsl_matrix* X = gsl_matrix_alloc(n, p);
	FILE * f = fopen("erdata.txt", "r");
	gsl_matrix_fscanf(f, X);
	fclose(f);

	//Calculate covariance matrix
	gsl_matrix* covX = gsl_matrix_alloc(p, p);
	makeCovariance(covX,X);
	
	// Sample from the multivariate normal using the Cholesky decomposition
	gsl_matrix* samples = gsl_matrix_alloc(nsamples, p);
	randomMVN(r, samples, covX);

	// Store samples
	printmatrix("samples.txt", samples);

	// Free memory
	gsl_matrix_free(X);
	gsl_matrix_free(covX);
	gsl_matrix_free(samples);
  	gsl_rng_free(r);
	
  	return(1);
}
