// Allen Roberts
// June 1, 2020
// Stat 534
// Final Project

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrices.h"



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
		printf("\niter=%d", iter);

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
		printmatrix("newBeta.txt", newBeta);
	
		newLoglik = logisticLogLikStar(n, y, x, newBeta);
		printf("\n newLoglik=%f", newLoglik);
	
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

			gsl_matrix_free(newBeta);
			gsl_matrix_free(hessian);
			gsl_matrix_free(gradient);
			gsl_matrix_free(hessianInv);
			gsl_matrix_free(hessGrad);

			return(beta);
		}
	
		currentLoglik = newLoglik;

	}

	// MEMORY FREE HERE?
	
	return(beta);

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


	// Free memory
	gsl_matrix_free(x);
	gsl_matrix_free(y);
	gsl_matrix_free(betaMode);
	gsl_matrix_free(data);
	
  	return(1);
}
