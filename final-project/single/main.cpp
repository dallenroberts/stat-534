// Allen Roberts
// June 1, 2020
// Stat 534
// Final Project

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrices.h"



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
// Inputs random number state mystream, gsl_matrix * samples to hold samples, 
// gsl_matrix * sigma is the covariance matrix of the multivariate normal, and gsl_matrix*
// means is the px1 mean vector.
void randomMVN(gsl_rng* mystream, gsl_matrix* samples, gsl_matrix* sigma, gsl_matrix* means) {

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

		// Add mean vector
		printf("\nAbout to add matrices\n");
		gsl_matrix_add(X, means);

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

// Inverse logit function
double inverseLogit(double x) {

	return(exp(x)/(1.0+exp(x)));
}

// function for the computation of the Hessian
double inverseLogit2(double x) {

	return(exp(x)/pow(1.0+exp(x), 2.0));
}


// Computes pi_i = P(y_i = 1 | x_i)
gsl_matrix* getPi(int n, gsl_matrix* x, gsl_matrix* beta) {

	gsl_matrix* x0 = gsl_matrix_alloc(n, 2);
	gsl_matrix* out = gsl_matrix_alloc(n, 1);

	int i;

	// Initialize model matrix
	for(i=0;i<n;i++) {
		gsl_matrix_set(x0, i, 0, 1.0); // Intercept column
		gsl_matrix_set(x0, i, 1, gsl_matrix_get(x, i, 0)); // Values of predictor
	}

	// Matrix multiply x0 by beta
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, x0, beta, 0.0, out);

	// Inverse logit transform output
	for(i=0;i<n;i++) {

		gsl_matrix_set(out, i, 0, inverseLogit(gsl_matrix_get(out, i, 0)));
	}

	gsl_matrix_free(x0);

	return(out);

}


// another function for the computation of the Hessian
gsl_matrix* getPi2(int n, gsl_matrix* x, gsl_matrix* beta) {

	gsl_matrix* x0 = gsl_matrix_alloc(n, 2);
	gsl_matrix* out = gsl_matrix_alloc(n, 1);

	int i;

	// Initialize model matrix
	for(i=0;i<n;i++) {
		gsl_matrix_set(x0, i, 0, 1.0); // Intercept column
		gsl_matrix_set(x0, i, 1, gsl_matrix_get(x, i, 0)); // Values of predictor
	}

	// Matrix multiply x0 by beta
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, x0, beta, 0.0, out);

	// Inverse logit transform output
	for(i=0;i<n;i++) {

		gsl_matrix_set(out, i, 0, inverseLogit2(gsl_matrix_get(out, i, 0)));
	}

	gsl_matrix_free(x0);

	return(out);

}

// Logistic log-likelihood star (from Bayesian logistic regression eq 2.5)
double logisticLogLikStar(int n, gsl_matrix* y, gsl_matrix* x, gsl_matrix* beta) {

	double logLik = 0;
	int i,j;
	double yi;
	double pi;
	double lstar;
	double bsum = 0;

	// getPis
	gsl_matrix* Pis = getPi(n, x, beta);
	printmatrix("pis.txt", Pis);

	// Calculate logistic log likelihood
	for(i=0;i<n;i++) {

		yi = gsl_matrix_get(y, i, 0);
		pi = gsl_matrix_get(Pis, i, 0);

		logLik += yi*log(pi) + (1.0-yi)*log(1.0-pi);
	}

	//Calculate sum of squared entries in beta matrix
	for(i=0;i<(beta->size1);i++) {
		for(j=0;j<(beta->size2);j++) {

			bsum += pow(gsl_matrix_get(beta, i, j), 2.0);
		}
	}

	// Calculate lstar
	lstar = -1.0*log(2.0*M_PI) - 0.5*bsum + logLik;

	// printf("\nlogLik = %f", logLik);
	// printf("\nlstar = %f \n", lstar);

	gsl_matrix_free(Pis);

	return(lstar);
}

// Obtain the Hessian for Newton-Raphson
gsl_matrix* getHessian(int n, gsl_matrix* x, gsl_matrix* beta) {

	int i;
	double pi, xi;

	gsl_matrix* hessian = gsl_matrix_alloc(2, 2);

	double h_00 = 0;
	double h_01 = 0;
	double h_10 = 0;
	double h_11 = 0;

	// Get Pi2
	gsl_matrix* Pi2 = getPi2(n, x, beta);

	// Update hessian entries
	for(i=0;i<n;i++) {

		pi = gsl_matrix_get(Pi2, i, 0);
		xi = gsl_matrix_get(x, i, 0);

		h_00 += pi;
		h_01 += pi*xi;
		h_11 += pi*pow(xi, 2.0);

	}

	h_00 += 1.0;
	h_10 = h_01;
	h_11 += 1.0;

	gsl_matrix_set(hessian, 0, 0, -h_00);
	gsl_matrix_set(hessian, 1, 0, -h_10);
	gsl_matrix_set(hessian, 0, 1, -h_01);
	gsl_matrix_set(hessian, 1, 1, -h_11);

	gsl_matrix_free(Pi2);

	return(hessian);

}

// Obtain the gradient for Newton-Raphson
gsl_matrix* getGradient(int n, gsl_matrix* y, gsl_matrix* x, gsl_matrix* beta) {

	int i;
	double pi, xi, yi;

	gsl_matrix* gradient = gsl_matrix_alloc(2, 1);

	double g_00 = 0;
	double g_10 = 0;

	// Get Pis
	gsl_matrix* Pi = getPi(n, x, beta);

	// Update gradient entries
	for(i=0;i<n;i++) {

		pi = gsl_matrix_get(Pi, i, 0);
		xi = gsl_matrix_get(x, i, 0);
		yi = gsl_matrix_get(y, i, 0);

		g_00 += (yi - pi);
		g_10 += (yi - pi)*xi;

	}

	g_00 -= gsl_matrix_get(beta, 0, 0);
	g_10 -= gsl_matrix_get(beta, 1, 0);

	gsl_matrix_set(gradient, 0, 0, g_00);
	gsl_matrix_set(gradient, 1, 0, g_10);

	gsl_matrix_free(Pi);

	return(gradient);
}

// this function implements our own Newton-Raphson procedure
gsl_matrix* getcoefNR(int n, gsl_matrix* y, gsl_matrix* x, int maxIter = 1000) {

	double currentLoglik;
	double newLoglik;
	int iter = 0;
	double tol = 0.00001;

	// Initialize beta matrix
	gsl_matrix* beta = gsl_matrix_alloc(2, 1);
	gsl_matrix_set_zero(beta);

	gsl_matrix* newBeta = gsl_matrix_alloc(2, 1);

	// Calculate log likelihood l*
	currentLoglik = logisticLogLikStar(n, y, x, beta);

	// Matrix to store product of hessian inverse and gradient
	gsl_matrix* hessGrad = gsl_matrix_alloc(2, 1);

	// Infinite loop unless we stop it someplace inside
	while(iter < maxIter) {

		iter += 1;
		// printf("\niter=%d", iter);

		// Get Hessian
		gsl_matrix* hessian = getHessian(n, x, beta);
		// printmatrix("hessian.txt", hessian);
	
		// Get gradient
		gsl_matrix* gradient = getGradient(n, y, x, beta);
		// printmatrix("gradient.txt", gradient);
	
		// Get hessian inverse
		gsl_matrix* hessianInv = inverse(hessian);
		// printmatrix("hessianInv.txt", hessianInv);
	
		// Get product of hessian inverse and gradient
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, hessianInv, gradient, 0.0, hessGrad);
		// printmatrix("hessGrad.txt", hessGrad);
	
		// Update new beta
		gsl_matrix_memcpy(newBeta, beta);
		gsl_matrix_sub(newBeta, hessGrad);
		// printmatrix("newBeta.txt", newBeta);
	
		newLoglik = logisticLogLikStar(n, y, x, newBeta);
		// printf("\n newLoglik=%f", newLoglik);
	
		// At each iteration the log-likelihood must increase
		if(newLoglik < currentLoglik) {
			printf("Coding error! Iteration made logLik worse");
			exit(1);
		}
	
		// Update beta
		gsl_matrix_memcpy(beta, newBeta);
		// printmatrix("beta.txt", beta);
	
		// Stop if the log-likelihood does not improve by too much
		if((newLoglik - currentLoglik) < tol) {

			printf("\n NR algorithm converged after %d iterations.", iter);
			printf(" Log-likelihood is %f \n", currentLoglik);


			// IS THIS MEMORY FREE NECESSARY?
			gsl_matrix_free(newBeta);
			gsl_matrix_free(hessian);
			gsl_matrix_free(gradient);
			gsl_matrix_free(hessianInv);
			gsl_matrix_free(hessGrad);

			return(beta);
		} else {

			currentLoglik = newLoglik;

		}

	}

	// MEMORY FREE HERE?
	if(iter == maxIter) {

		printf("\nNR algorithm reached maximum iterations.");
		printf(" Log-likelihood is %f \n", currentLoglik);
		return(beta);

	}

}

// sampleMH performs Metropolis-Hastings sampling from the posterior distribution
// P(beta0, beta1 | D) of Bayesian univariate logistic regression. It returns
// a matrix with two columns (beta0 and beta1) and niter rows (one for each)
// sample.
gsl_matrix* sampleMH(gsl_rng* mystream, int n,  gsl_matrix* y, gsl_matrix* x, gsl_matrix* betaMode, int niter) {

	int k;

	// Allocate matrix to store samples
	gsl_matrix* samples = gsl_matrix_alloc(niter, 2);
	gsl_matrix_set_zero(samples);

	// Proposal distribution covariance matrix
	gsl_matrix* hessian = getHessian(n, x, betaMode);
	gsl_matrix* covMat = inverse(hessian);
	gsl_matrix_scale(covMat, -1.0);
	printmatrix("covMat.txt", covMat);

	// Initial state
	gsl_matrix* currentBeta = gsl_matrix_alloc(2,1);
	gsl_matrix_memcpy(currentBeta, betaMode);

	gsl_matrix* candidateBeta = gsl_matrix_alloc(2, 1);

	// Start Markov chain
	// for(k=0; k<niter;k++) {
	printf("\nMade it to randomMVN\n");
		// printf("\n MC iter %d", (k+1));
		randomMVN(mystream, candidateBeta, covMat, currentBeta);
		printmatrix("candidateBeta.txt", candidateBeta);


	// }


	// Free memory
	gsl_matrix_free(covMat);
	gsl_matrix_free(currentBeta);
	gsl_matrix_free(candidateBeta);

	return(samples);

}

// Loads 534finalprojectdata.txt
int main() {

	int i;
  	int n = 148;
  	int p = 61;
  	int response = 60; // Index of the response column
	
	int index = 0;

	// Initialize random number generator
  	const gsl_rng_type* T;
  	gsl_rng* r;
	
  	gsl_rng_env_setup();
	
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);

  	// Loads 534finalprojectdata.txt. This file has 148 rows (samples) and 61 columns (variables). 
  	// The first 60 columns are associated with 60 explanatory variables X, 
  	// while column 61 (the last column) corresponds with the response binary variable Y
  	gsl_matrix* data = gsl_matrix_alloc(n, p);
	FILE * f = fopen("534finalprojectdata.txt", "r");
	gsl_matrix_fscanf(f, data);
	fclose(f);

	// Initialize predictor and response columns
	gsl_matrix* x = gsl_matrix_alloc(n, 1);
	gsl_matrix* y = gsl_matrix_alloc(n, 1);
	for(i=0;i<n;i++) {

		gsl_matrix_set(x, i, 0, gsl_matrix_get(data, i, index));
		gsl_matrix_set(y, i, 0, gsl_matrix_get(data, i, response));

	}

	// Calculate beta modes using Newton-Raphson algorithm
	gsl_matrix* betaMode = getcoefNR(n, y, x, 1000);
	printmatrix("betaMode.txt", betaMode);

	gsl_matrix* samples = sampleMH(r, n, y, x, betaMode, 10);

	// Free memory
	gsl_matrix_free(x);
	gsl_matrix_free(y);
	gsl_matrix_free(betaMode);
	gsl_matrix_free(samples);
	gsl_matrix_free(data);
	gsl_rng_free(r);
	
  	return(1);
}
