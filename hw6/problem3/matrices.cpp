#include "matrices.h"

//prints the elements of a matrix in a file
void printmatrix(char* filename,gsl_matrix* m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<m->size1;i++)
	{
	        fprintf(out,"%.3lf",gsl_matrix_get(m,i,0));
		for(j=1;j<m->size2;j++)
		{
			fprintf(out,"\t%.3lf",
				gsl_matrix_get(m,i,j));
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//creates the transpose of the matrix m
gsl_matrix* transposematrix(gsl_matrix* m)
{
	int i,j;
	
	gsl_matrix* tm = gsl_matrix_alloc(m->size2,m->size1);
	
	for(i=0;i<tm->size1;i++)
	{
		for(j=0;j<tm->size2;j++)
		{
		  gsl_matrix_set(tm,i,j,gsl_matrix_get(m,j,i));
		}
	}	
	
	return(tm);
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(gsl_matrix* m1,gsl_matrix* m2,gsl_matrix* m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<m->size1;i++)
	{
	  for(k=0;k<m->size2;k++)
	  {
	    s = 0;
	    for(j=0;j<m1->size2;j++)
	    {
	      s += gsl_matrix_get(m1,i,j)*gsl_matrix_get(m2,j,k);
	    }
	    gsl_matrix_set(m,i,k,s);
	  }
	}
	return;
}


//computes the inverse of a positive definite matrix
//the function returns a new matrix which contains the inverse
//the matrix that gets inverted is not modified
gsl_matrix* inverse(gsl_matrix* K)
{
	int j;
	
	gsl_matrix* copyK = gsl_matrix_alloc(K->size1,K->size1);
	if(GSL_SUCCESS!=gsl_matrix_memcpy(copyK,K))
	{
		printf("GSL failed to copy a matrix.\n");
		exit(1);
	}
	
	gsl_matrix* inverse = gsl_matrix_alloc(K->size1,K->size1);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);
	
	if(GSL_SUCCESS!=gsl_linalg_LU_decomp(copyK,myperm,&j))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	if(GSL_SUCCESS!=gsl_linalg_LU_invert(copyK,myperm,inverse))
	{
		printf("GSL failed matrix inversion.\n");
		exit(1);
	}
	gsl_permutation_free(myperm);
	gsl_matrix_free(copyK);
	
	return(inverse);
}

//creates a submatrix of matrix M
//the indices of the rows and columns to be selected are
//specified in the last four arguments of this function
gsl_matrix* MakeSubmatrix(gsl_matrix* M,
			  int* IndRow,int lenIndRow,
			  int* IndColumn,int lenIndColumn)
{
	int i,j;
	gsl_matrix* subM = gsl_matrix_alloc(lenIndRow,lenIndColumn);
	
	for(i=0;i<lenIndRow;i++)
	{
		for(j=0;j<lenIndColumn;j++)
		{
			gsl_matrix_set(subM,i,j,
                                       gsl_matrix_get(M,IndRow[i],IndColumn[j]));
		}
	}
	
	return(subM);
}


//computes the log of the determinant of a symmetric positive definite matrix
double logdet(gsl_matrix* K)
{
        int i;

	gsl_matrix* CopyOfK = gsl_matrix_alloc(K->size1,K->size2);
	gsl_matrix_memcpy(CopyOfK,K);
	gsl_permutation *myperm = gsl_permutation_alloc(K->size1);
	if(GSL_SUCCESS!=gsl_linalg_LU_decomp(CopyOfK,myperm,&i))
	{
		printf("GSL failed LU decomposition.\n");
		exit(1);
	}
	double logdet = gsl_linalg_LU_lndet(CopyOfK);
	gsl_permutation_free(myperm);
	gsl_matrix_free(CopyOfK);
	return(logdet);
}

{
	gsl_vector* tempvec = gsl_vector_alloc(n);
	int i;

	//set the columns of the matrix.
	for(i = 0; i < lenA; i++)
	{
		gsl_matrix_get_col(tempvec, fulldata, (A[i]-1));
		gsl_matrix_set_col(subdata, i, tempvec);
	}
	// clean the memory
	gsl_vector_free(tempvec);
	return;
}

gsl_matrix* create_MA(gsl_matrix* D_A)
{
	int j = 0;
	int len_A = D_A->size2;
	int n = D_A->size1;
	double x;

	gsl_matrix* D_A_transpose = gsl_matrix_alloc(len_A, n);
	gsl_matrix_transpose_memcpy(D_A_transpose, D_A);

	gsl_matrix* M_A = gsl_matrix_alloc(len_A, len_A);

	// Here, compute (X^T X), store in matrix MA. The function gsl_blas_dgemm
	// is a general matrix multiplication algorithm, and is several times faster than
	// the usual nested for-loop routine that you might use naively.

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D_A_transpose, D_A, 0.0, M_A);
	gsl_matrix_free(D_A_transpose);

	for(j = 0; j < len_A;j ++) // Add the identity matrix
	{
		x = gsl_matrix_get(M_A, j, j) + 1.0;
		gsl_matrix_set(M_A, j , j , x);
	}

	return(M_A);
}

double compute_ml_matprod(gsl_vector* response, gsl_matrix* M_A, gsl_matrix* D_A)
{
	double matprod = 0.0;
	int lenA = D_A->size2;
	int n = D_A->size1;

	// allocate temporary vectors for intermediate calculations
	gsl_vector* tempvec = gsl_vector_calloc(lenA);
	gsl_vector* tempvec2 = gsl_vector_calloc(lenA);
	gsl_vector* tempvec3 = gsl_vector_calloc(n);
	gsl_matrix* choldecomp = gsl_matrix_alloc(lenA, lenA);
	gsl_matrix_memcpy(choldecomp, M_A);

	// product of t(D_A) %*% D_1, stored in tempvec
	gsl_blas_dgemv(CblasTrans, 1.0, D_A, response, 0.0, tempvec);

	// next two lines, tempvec2 stores inv(MA) %*% t(D_A) %*% D_1,
	// computed via the cholesky decomposition, since MA is pos-definite.
	gsl_linalg_cholesky_decomp(choldecomp);
	gsl_linalg_cholesky_solve(choldecomp, tempvec, tempvec2);

	// product of D_A %*% inv(MA) %*% t(D_A) %*% D_1, which is a vector of length
	// n, stored in tempvec3
	gsl_blas_dgemv(CblasNoTrans, 1.0, D_A, tempvec2, 0.0, tempvec3);

	// finally, t(D_1) %*% D_A %*% inv(MA) %*% t(D_A) %*% D_1
	gsl_blas_ddot(response, tempvec3, &matprod);

	gsl_matrix_free(choldecomp);
	gsl_vector_free(tempvec);
	gsl_vector_free(tempvec2);
	gsl_vector_free(tempvec3);

	return(matprod);
}

double marglik(gsl_matrix* data, int lenA, int* A)
{
	// declare variables for intermediate calculations
	double gam_coef = 0.0;
	double l_det = 0.0;
	double matprod = 0.0;
	double inner_d = 0.0;
	double total = 0.0;
	int n = data->size1;

	// extract submatrix D_A and create the matrix M_A
	gsl_matrix* D_A = gsl_matrix_alloc(n, lenA);
	subset_gsl_matrix(data, D_A, n, A, lenA);
	gsl_matrix* M_A = create_MA(D_A);

	// extract response vector
	gsl_vector* response = gsl_vector_alloc(n);
	gsl_matrix_get_col(response, data, 0);


	//get the log gamma piece
	gam_coef = lgamma((n + 2.0 + lenA)/2.0) - lgamma((lenA + 2.0)/2.0);
	
	// get the log determinant of MA
	l_det = logdet(M_A);

	// get the inner product of the response, ie compute (y^t y)
	gsl_blas_ddot(response, response, &inner_d);
	
	// compute the long matrix product
	matprod = compute_ml_matprod(response, M_A, D_A);
	
	total = gam_coef - .5*l_det - ((n + lenA + 2.0)/2.0)*log(1.0 + inner_d - matprod);

	// deallocate memory
	gsl_matrix_free(D_A);
	gsl_vector_free(response);
	gsl_matrix_free(M_A);

	return(total);
}

