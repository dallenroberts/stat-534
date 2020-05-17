#include "matrices.h"

double getDeterminant(int nrow)
{

	int n;
	double det;
	int i,j,k;

	printf("Hello (matrix) world.\n");

	// Read in banded matrix
	gsl_matrix * m = gsl_matrix_alloc(nrow, nrow);
	FILE * f = fopen("mybandedmatrix.txt", "r");
	gsl_matrix_fscanf(f, m);
	fclose(f);

	n = 3;

	// Edge cases
	if(n == 1) {

		return(gsl_matrix_get(m, 0,0));

	}
	if(n == 2) {

		return(gsl_matrix_get(m,0,0)*gsl_matrix_get(m,1,1) - gsl_matrix_get(m,0,1)*gsl_matrix_get(m,1,0));

	}

	// Define minor matrix by removing the first row and the kth column (where k ranges from 0 to n-1) from an nxn matrix m
	k = 3;
	gsl_matrix * minor = gsl_matrix_alloc(nrow-1, nrow-1);

	for(i=0;i<(nrow-1);i++) {

		printf("\ni = %d", i);


		for(j=0; j<(nrow-1);j++) {

			printf("\nj = %d", j);

			if(j<k){

				gsl_matrix_set(minor, i+1, j, gsl_matrix_get(m, i+1, j));

			} else {

				gsl_matrix_set(minor, i+1, j+1, gsl_matrix_get(m, i+1, j+1));

			}
		}
	}

	printmatrix("minor.mat", minor);

	return(5);

}


int main()
{
	
	// Number of rows of the symmetric matrix for which we want to calculate the determinant
	int n = 10;
	double det;

	det = getDeterminant(n);

	printf("The determinant of the matrix is: %.4lf\n", det);

	return(1);
}
