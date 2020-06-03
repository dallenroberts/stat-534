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
// Inputs random number state mystream, px1 gsl_matrix * samples to hold samples, 
// gsl_matrix * sigma is the pxp covariance matrix of the multivariate normal, and gsl_matrix*
// means is the px1 mean vector.
void randomMVN(gsl_rng* mystream, gsl_matrix* samples, gsl_matrix* sigma, gsl_matrix* means) {

	int i,j;
	gsl_matrix* z = gsl_matrix_alloc(sigma->size2, 1);

	// Calculate Cholesky decomposition
	gsl_matrix* psi = makeCholesky(sigma);

	// Draw samples
	gsl_matrix* X = gsl_matrix_alloc(sigma->size2, 1); // px1 matrix output by dgemm
	gsl_vector* s = gsl_vector_alloc(sigma->size2); // vector to hold samples

	for(i = 0; i < samples->size2; i++) {

		// Generate p independent N(0,1) random numbers
		for(j = 0; j < sigma->size1; j++) {

			gsl_matrix_set(z, j, 0, gsl_ran_ugaussian(mystream));
		}

		// Calculate matrix product X = psi*Z. Note that this is a px1 matrix
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, psi, z, 0.0, X);

		// Add mean vector
		// printf("\nAbout to add matrices\n");
		// printmatrix("X.txt", X);
		// printmatrix("means.txt", means);
		gsl_matrix_add(X, means);

		// Store in samples matrix
		// printf("\n Storing in samples matrix \n");
		// printmatrix("sigma.txt", sigma);
		gsl_matrix_get_col(s, X, 0);
		gsl_matrix_set_col(samples, i, s);
		// printmatrix("samples.txt", samples);

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
void getHessian(int n, gsl_matrix* x, gsl_matrix* beta, gsl_matrix* hessian) {

	int i;
	double pi, xi;

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

}

// Obtain the gradient for Newton-Raphson
void getGradient(int n, gsl_matrix* y, gsl_matrix* x, gsl_matrix* beta, gsl_matrix* gradient) {

	int i;
	double pi, xi, yi;

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

	// Matrix to store hessian
	gsl_matrix* hessian = gsl_matrix_alloc(2, 2);

	// Matrix to store gradient
	gsl_matrix* gradient = gsl_matrix_alloc(2, 1);

	// Matrix to store the hessian inverse
	gsl_matrix* hessianInv = gsl_matrix_alloc(2, 2);

	// Matrix to store product of hessian inverse and gradient
	gsl_matrix* hessGrad = gsl_matrix_alloc(2, 1);

	while(iter < maxIter) {

		iter += 1;
		// printf("\niter=%d", iter);

		// Get Hessian
		getHessian(n, x, beta, hessian);
		// printmatrix("hessian.txt", hessian);
	
		// Get gradient
		getGradient(n, y, x, beta, gradient);
		// printmatrix("gradient.txt", gradient);
	
		// Get hessian inverse
		inverse2(hessian, hessianInv);

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
			printf("\n    l*(beta0, beta1) = %.3f \n", currentLoglik);


			// IS THIS MEMORY FREE NECESSARY?
			gsl_matrix_free(newBeta);

			break;

		} else {

			currentLoglik = newLoglik;

		}

	}

	if(iter == maxIter) {

		printf("\nNR algorithm reached maximum iterations.");
		printf("\n    l*(beta0, beta1) = %.3f \n", currentLoglik);

	}

	// Free memory
	gsl_matrix_free(hessian);
	gsl_matrix_free(gradient);
	gsl_matrix_free(hessianInv);
	gsl_matrix_free(hessGrad);

	return(beta);
}

// sampleMH performs Metropolis-Hastings sampling from the posterior distribution
// P(beta0, beta1 | D) of Bayesian univariate logistic regression. It returns
// a matrix with two columns (beta0 and beta1) and niter rows (one for each)
// sample.
gsl_matrix* sampleMH(gsl_rng* mystream, int n,  gsl_matrix* y, gsl_matrix* x, gsl_matrix* betaMode, int niter) {

	int k;
	double score;
	double u;

	// Allocate matrix to store samples
	gsl_matrix* samples = gsl_matrix_alloc(niter, 2);
	gsl_matrix_set_zero(samples);

	// Proposal distribution covariance matrix
	gsl_matrix* hessian = gsl_matrix_alloc(2, 2);
	getHessian(n, x, betaMode, hessian);

	gsl_matrix* covMat = inverse(hessian);
	gsl_matrix_scale(covMat, -1.0);
	printmatrix("covMat.txt", covMat);

	// Initial state
	gsl_matrix* currentBeta = gsl_matrix_alloc(2,1);
	gsl_matrix_memcpy(currentBeta, betaMode);

	gsl_matrix* candidateBeta = gsl_matrix_alloc(2, 1);

	// Start Markov chain
	printf("\n Starting Metropolis-Hastings Algorithm...");
	for(k=0; k<niter;k++) {
		// printf("\n MC iter %d", (k+1));
		
		// Draw candidate beta from multivariate normal
		randomMVN(mystream, candidateBeta, covMat, currentBeta);
		// printmatrix("candidateBeta.txt", candidateBeta);

		// Accept or reject candidate beta
		score = logisticLogLikStar(n, y, x, candidateBeta) - logisticLogLikStar(n, y, x, currentBeta);
		// printf("\n score = %f \n", score);

		if(score >= 1.0) {

			gsl_matrix_memcpy(currentBeta, candidateBeta);
		
		} else {

			u = gsl_ran_flat(mystream, 0.0, 1.0);

			if(log(u) <= score) {

				gsl_matrix_memcpy(currentBeta, candidateBeta);

			}
		}

		// Update chain
		gsl_matrix_set(samples, k, 0, gsl_matrix_get(currentBeta, 0, 0));
		gsl_matrix_set(samples, k, 1, gsl_matrix_get(currentBeta, 1, 0));

	}

	// Free memory
	gsl_matrix_free(hessian);
	gsl_matrix_free(covMat);
	gsl_matrix_free(currentBeta);
	gsl_matrix_free(candidateBeta);

	return(samples);

}

// Calculates the posterior means of niter samples from the joint distribution
// of the betas given the observed data. Sampling is done via Metropolis-
// Hastings. Returns a 2x1 matrix of beta values.
gsl_matrix* getPosteriorMeans(gsl_rng* mystream, int n,  gsl_matrix* y, gsl_matrix* x, gsl_matrix* betaMode, int niter) {

	int i;
	gsl_vector_view a;

	// Simulate 1000 samples from the posterior distribution using Metropolis-Hastings
	gsl_matrix* samples = sampleMH(mystream, n, y, x, betaMode, niter);
	// printmatrix("MHsamples.txt", samples);

	// Calculate sample means
	printf("Metropolis-Hastings finished.\n");
	gsl_matrix* sampleMeans = gsl_matrix_alloc(2, 1);

	for(i=0; i<(samples->size2); i++) {

		a = gsl_matrix_column(samples, i);
		gsl_matrix_set(sampleMeans, i, 0, 
			gsl_stats_mean(a.vector.data, a.vector.stride, (samples->size1)));

	}

	gsl_matrix_free(samples);

	return(sampleMeans);
}

// getLaplaceApprox uses the Laplace Approximation to calcuate an approximate 
// marginal likelihood for univariate Bayesian logistic regression. Note that 
// this function calculates and returns the log-likehood
double getLaplaceApprox(int n, gsl_matrix* y, gsl_matrix* x, gsl_matrix* betaMode) {

	double loglik_star;
	double lml;

	// Get Hessian and multiply by -1
	gsl_matrix* hessian = gsl_matrix_alloc(2,2);
	getHessian(n, x, betaMode, hessian);
	gsl_matrix_scale(hessian, -1.0);

	// Calculate l*(beta_0, beta_1)
	loglik_star = logisticLogLikStar(n, y, x, betaMode);

	// Calculate log marginal likelihood
	lml = 2.0*log(M_PI) + loglik_star - 0.5*logdet(hessian);

	// Free memory
	gsl_matrix_free(hessian);

	return(lml);

}

// This function evaluates the log marginal likelihood P(D) using 
// Monte Carlo integration. It draws nsamples from the two independent
// normal priors for the regression coefficients beta0 and beta
// and then calculates the average likelihood exp(l*(beta0, beta1))
double getMC(gsl_rng* mystream, int n, gsl_matrix* y, gsl_matrix* x, int nsamples) {

	int j;
	double sum = 0;
	double lml;

	// Means and covariance matrix for prior distributions
	gsl_matrix* priorMeans = gsl_matrix_alloc(2, 1);
	gsl_matrix_set_zero(priorMeans);

	gsl_matrix* priorCovMat = gsl_matrix_alloc(2, 2);
	gsl_matrix_set_identity(priorCovMat);

	// Placeholder for sampled betas
	gsl_matrix* beta_j = gsl_matrix_alloc(2, 1);

	for(j=0;j<nsamples;j++) {

		// Sample from prior distribution
		randomMVN(mystream, beta_j, priorCovMat, priorMeans);

		// Calculate likelihood
		sum += exp(logisticLogLikStar(n,y,x,beta_j));

	}

	lml = sum/(double)nsamples;

	// Free memory
	gsl_matrix_free(priorMeans);
	gsl_matrix_free(priorCovMat);
	gsl_matrix_free(beta_j);

	return(lml);
}

// Loads 534finalprojectdata.txt
int main() {

	int i;
  	int n = 148;
  	int p = 61;
  	int response = 60; // Index of the response column
  	double lml_la;
  	double lml_mc;
	
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

	// Calculate posterior means for betas
	gsl_matrix* sampleMeans = getPosteriorMeans(r, n, y, x, betaMode, 1000);

	printf("Sample means:\n");
	for(i=0;i<(sampleMeans->size1);i++) {

		printf("   beta%i = %.3f\n", i, gsl_matrix_get(sampleMeans, i, 0));
	}
	printmatrix("sampleMeans.txt", sampleMeans);

	// Calculate log marginal likelihood using LaPlace approximation
	lml_la = getLaplaceApprox(n, y, x, betaMode);
	printf("\n Posterior log marginal likelihood P(D) calculated by Laplace approximation = %.3f \n", lml_la);

	// Calculate log marginal likelihood using Monte Carlo integration
	lml_mc = getMC(r, n, y, x, 10000);
	printf("\n Posterior log marginal likelihood P(D) calculated by Monte Carlo integration = %.3f \n", lml_mc);

	
	// Free memory
	gsl_matrix_free(x);
	gsl_matrix_free(y);
	gsl_matrix_free(betaMode);
	gsl_matrix_free(sampleMeans);
	gsl_matrix_free(data);
	gsl_rng_free(r);
	
  	return(1);
}
