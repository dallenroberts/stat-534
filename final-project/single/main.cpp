// Allen Roberts
// June 1, 2020
// Stat 534
// Final Project

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrices.h"

// Function that inputs a pointer to a pxp GSL matrix covX and a pointer to a 
// nxp data matrix X and calculates the pxp sample covariance matrix sigma
// Outputs nothing, but updates the gsl_matrix stored at location covX
void makeCovariance(gsl_matrix* covX, gsl_matrix* X) {

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
// Note that the final matrix returned is lower triangular
gsl_matrix* makeCholesky(gsl_matrix* K) {

	int i, j;

	gsl_matrix* cholesky = gsl_matrix_alloc(K->size1,K->size2);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(cholesky,K)) {
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}

	if(GSL_SUCCESS!=gsl_linalg_cholesky_decomp(cholesky)) {
		printf("GSL failed cholesky decomposition.\n");
		exit(1);
	}

	// Retain only the lower triangle
	for(i=0;i<cholesky->size1;i++) {
		for(j=(i+1);j<cholesky->size2;j++) {
			gsl_matrix_set(cholesky, i, j, 0.0);
		}
	}

	return(cholesky);
}

// Samples from the multivariate normal distribution using the Cholesky decomposition
// Inputs random number state mystream, gsl_matrix * samples to hold samples, and 
// gsl_matrix * sigma is the covariance matrix of the multivariate normal
void randomMVN(gsl_rng* mystream, gsl_matrix* samples,gsl_matrix* sigma) {

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

			gsl_matrix_set(z, j, 0, gsl_ran_ugaussian(mystream));
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

double logisticLogLikStar(int y, int x, gsl_matrix* beta) {

	printf("\n Inside logisticLogLikStar\n");
	return(1);
}

// Loads 534finalprojectdata.txt
int main() {

	int i;
  	int n = 148;
  	int p = 61;
  	int response = 60; // Index of the response column
	
	int index = 0;

  	// Loads 534finalprojectdata.txt. This file has 148 rows (samples) and 61 columns (variables). 
  	// The first 60 columns are associated with 60 explanatory variables X, 
  	// while column 61 (the last column) corresponds with the response binary variable Y
  	gsl_matrix* data = gsl_matrix_alloc(n, p);
	FILE * f = fopen("534finalprojectdata.txt", "r");
	gsl_matrix_fscanf(f, data);
	fclose(f);

	// Calculate log likelihood l*
	double l_star;

	// Initialize beta matrix
	gsl_matrix* beta = gsl_matrix_alloc(2, 1);
	gsl_matrix_set_zero(beta);

	l_star = logisticLogLikStar(response, index, beta);



	// Free memory
	gsl_matrix_free(beta);
	gsl_matrix_free(data);
	
  	return(1);
}
