#include "matrices.h"

double getDeterminant(int nrow)
{

	int n;
	double det;

	printf("Hello (matrix) world.\n");

	// Read in banded matrix
	gsl_matrix * m = gsl_matrix_alloc(nrow, nrow);
	FILE * f = fopen("mybandedmatrix.txt", "r");
	gsl_matrix_fscanf(f, m);
	fclose(f);

	n = 1;

	// Edge cases
	if(n == 1) {

		det = gsl_matrix_get(m, 0,0);

	} else if(n == 2) {

		det = gsl_matrix_get(m,0,0)*gsl_matrix_get(m,1,1) - gsl_matrix_get(m,0,1)*gsl_matrix_get(m,1,0);

	}
	// printmatrix("test.mat", m);

	return(det);

}


int main()
{
	
	int n = 10;
	double det;

	det = getDeterminant(n);

	printf("The determinant of the matrix is: %.4lf\n", det);

	return(1);
}
