#include "matrices.h"

double getDeterminant(int n)
{

	printf("Hello (matrix) world.\n");

	gsl_matrix * m = gsl_matrix_alloc(n, n);
	FILE * f = fopen("mybandedmatrix.txt", "r");

	gsl_matrix_fread(f, m);

	return(1)

}


int main()
{
	
	int n = 10;
	double det;

	det = getDeterminant(n);

	printf("The determinant of the matrix is: %.4lf\n",
               det);

	return(1);
}
