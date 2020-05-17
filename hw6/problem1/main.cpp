#include "matrices.h"

int main()
{
	printf("Hello (GSL matrix) world.\n");
	
	int i,j;
	int n = 10;
	
	//initialise a square matrix of dimension n
	gsl_matrix* m = gsl_matrix_alloc(n,n);
	
	FILE* out = fopen("mymatrixonecolumn.txt","w");
	gsl_matrix_fprintf(out,m,"%.2lf");
	fclose(out);
	
	//this does not look very nice, so
	//we'll do it again with our own function
	printmatrix("mymatrix.txt",m);

	//now make a banded matrix
        gsl_matrix_set(m,0,0,1);
	gsl_matrix_set(m,0,1,0.25);
	for(i=1;i<n-1;i++)
	{
	  gsl_matrix_set(m,i,i-1,0.25);
	  gsl_matrix_set(m,i,i,1);
	  gsl_matrix_set(m,i,i+1,0.25);
	}
	gsl_matrix_set(m,n-1,n-2,0.25);
	gsl_matrix_set(m,n-1,n-1,1);
	
	//print this matrix
	printmatrix("mybandedmatrix.txt",m);
	
	//tranpose the banded matrix
	gsl_matrix* tm = transposematrix(m);
	
	printmatrix("transposedbandedmatrix.txt",tm);
	
	//calculate the dot product of m and its transpose tm
	//first allocate a new matrix with the same dimensions
	gsl_matrix* dm = gsl_matrix_alloc(n,n);
	//copy the matrix m in the copy
	gsl_matrix_memcpy(dm,m);
	//finally, calculate the dot product
	//the result will be in "dm"
	//remark that the initial elements of "dm" are lost
	//after computing the dot product
	gsl_matrix_mul_elements(dm,tm);
	
	printmatrix("dotproductmatrix.txt",dm);

	//calculate the inverse of m
	gsl_matrix* mInverse = inverse(m);
	printmatrix("inversematrix.txt",mInverse);

	//calculate the product of m and its inverse
        gsl_matrix* A = gsl_matrix_alloc(n,n);
	matrixproduct(m,mInverse,A);
	printmatrix("productmatrix.txt",A);

        printf("The log of the determinant of m is: %.4lf\n",
               logdet(m));

	//free the memory
	gsl_matrix_free(m);
	gsl_matrix_free(tm);
	gsl_matrix_free(dm);
	gsl_matrix_free(mInverse);
	gsl_matrix_free(A);

	return(1);
}
