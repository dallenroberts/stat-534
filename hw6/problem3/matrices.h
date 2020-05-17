#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_errno.h>

extern "C" {
 void dpotri_(char*,int*,double*,int*,int*);
 void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*,
	 int*, double*, int*, double*, int*, int*);
}

//general functions
double ** allocmatrix(int n,int p);
void freematrix(int n,double** m);
void copymatrix(int n,int p,double** source,double** dest);
void readmatrix(char* filename,int n,int p,double* m[]);
void printmatrix(char* filename,int n,int p,double** m);
double** transposematrix(int n,int p,double** m);
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m);
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m);
void inverse(int p,double** m);
double logdet(int p,double** m);

//functions related to HW2
double** submatrix(int n,int p,double** data,int lenA,int* A);
double marglik(gsl_matrix* data, int lenA, int* A);

//
void subset_gsl_matrix(gsl_matrix* fulldata, gsl_matrix* subdata,
		int n, int* A, int lenA);
gsl_matrix* create_MA(gsl_matrix* D_A);

void subset_gsl_matrix(gsl_matrix* fulldata, gsl_matrix* subdata,
		int n, int* A, int lenA)
