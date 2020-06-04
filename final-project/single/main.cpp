// Allen Roberts
// June 1, 2020
// Stat 534
// Final Project

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "matrices.h"
#include "regmodels.h"
#include "bayes.h"

// Adds a regression using predictor index to the LPRegression list regressions
void bayesLogistic(int index, int n, int p, int response, gsl_rng* mystream, LPRegression regressions) {

	int i;
  	double lml_la;
  	double lml_mc;
  	int nMaxRegs = 5; // Maximum number of regressions to keep track of
  	int A[p-1]; //indices of the variables present in the regression
  	int lenA = -1; //number of indices

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
	// printmatrix("betaMode.txt", betaMode);

	// Calculate posterior means for betas
	gsl_matrix* sampleMeans = getPosteriorMeans(mystream, n, y, x, betaMode, 10000);

	printf(" Sample means:\n");
	for(i=0;i<(sampleMeans->size1);i++) {

		printf("    beta%i = %.3f\n", i, gsl_matrix_get(sampleMeans, i, 0));
	}
	// printmatrix("sampleMeans.txt", sampleMeans);

	// Calculate log marginal likelihood using LaPlace approximation
	printf("\n Posterior log marginal likelihood P(D) estimates:\n");
	lml_la = getLaplaceApprox(n, y, x, betaMode);
	printf("    Laplace approximation = %.3f \n", lml_la);

	// Calculate log marginal likelihood using Monte Carlo integration
	lml_mc = log(getMC(mystream, n, y, x, 100000));
	printf("    Monte Carlo integration = %.3f \n", lml_mc);

	// Add to linked list
	lenA = 1;
    A[0] = index+1;
    AddRegression(nMaxRegs, regressions,
      lenA, A, sampleMeans, lml_mc, lml_la);
    
    // Free memory
	gsl_matrix_free(x);
	gsl_matrix_free(y);
	gsl_matrix_free(betaMode);
	gsl_matrix_free(sampleMeans);
	gsl_matrix_free(data);

}

// Loads 534finalprojectdata.txt
int main() {

	int n = 148;
  	int p = 61;
  	int response = 60; // Index of response column
  	int i;

  	char outputfilename[] = "bestregressions.txt";

	// Initialize random number generator
  	const gsl_rng_type* T;
  	gsl_rng* r;
  	gsl_rng_env_setup();
  	T = gsl_rng_default;
  	r = gsl_rng_alloc(T);

  	//create the head of the list of regressions
  	LPRegression regressions = new Regression;
  	regressions->Next = NULL;

  	// Add regression
  	for(i=0;i<response;i++) {

  		bayesLogistic(i, n, p, response, r, regressions);

  	}

    //save the list in a file
  	SaveRegressions(outputfilename,regressions);

  	//delete all regressions
  	DeleteAllRegressions(regressions);
	
	// Free memory
	gsl_rng_free(r);
	delete regressions; regressions = NULL;
	
  	return(1);
}
