#include "matrices.h"

//allocates the memory for a matrix with 
//n rows and p columns
double ** allocmatrix(int n,int p)
{
	int i;
	double** m;
	
	m = new double*[n];
	for(i=0;i<n;i++)
	{
		m[i] = new double[p];
		memset(m[i],0,p*sizeof(double));
	}
	return(m);
}

//frees the memory for a matrix with n rows
void freematrix(int n,double** m)
{
	int i;
	
	for(i=0;i<n;i++)
	{
		delete[] m[i]; m[i] = NULL;
	}
	delete[] m; m = NULL;
	return;
}

//creates the copy of a matrix with n rows and p columns
void copymatrix(int n,int p,double** source,double** dest)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			dest[i][j] = source[i][j];
		}
	}
	return;
}

//reads from a file a matrix with n rows and p columns
void readmatrix(char* filename,int n,int p,double* m[])
{
	int i,j;
	double s;
	FILE* in = fopen(filename,"r");
	
	if(NULL==in)
	{
		printf("Cannot open input file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			fscanf(in,"%lf",&s);
			m[i][j] = s;
		}
	}
	fclose(in);
	return;
}

//prints the elements of a matrix in a file
void printmatrix(char* filename,int n,int p,double** m)
{
	int i,j;
	double s;
	FILE* out = fopen(filename,"w");
	
	if(NULL==out)
	{
		printf("Cannot open output file [%s]\n",filename);
		exit(1);
	}
	for(i=0;i<n;i++)
	{
		fprintf(out,"%.3lf",m[i][0]);
		for(j=1;j<p;j++)
		{
			fprintf(out,"\t%.3lf",m[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	return;
}

//creates the transpose of the matrix m
double** transposematrix(int n,int p,double** m)
{
	int i,j;
	
	double** tm = allocmatrix(p,n);
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<n;j++)
		{
			tm[i][j] = m[j][i];
		}
	}	
	
	return(tm);
}

//calculates the dot (element by element) product of two matrices m1 and m2
//with n rows and p columns; the result is saved in m
void dotmatrixproduct(int n,int p,double** m1,double** m2,double** m)
{
	int i,j;
	
	for(i=0;i<n;i++)
	{
		for(j=0;j<p;j++)
		{
			m[i][j] = m1[i][j]*m2[i][j];
		}
	}
	
	return;
}

//calculates the product of a nxp matrix m1 with a pxl matrix m2
//returns a nxl matrix m
void matrixproduct(int n,int p,int l,double** m1,double** m2,double** m)
{
	int i,j,k;
	double s;
	
	for(i=0;i<n;i++)
	{
		for(k=0;k<l;k++)
		{
			s = 0;
			for(j=0;j<p;j++)
			{
				s += m1[i][j]*m2[j][k];
			}
			m[i][k] = s;
		}
	}
	return;
}

void set_mat_identity(int p, double *A)
{
 int i;

 for(i = 0; i < p * p; i++) A[i] = 0;
 for(i = 0; i < p; i++) A[i * p + i] = 1;
 return;
}

//computes the inverse of a symmetric positive definite matrix
void inverse(int p,double** m)
{
  int i,j,k;
  double* m_copy = (double*)malloc((p * p) * sizeof(double));
  double* m_inv = (double*)malloc((p * p) * sizeof(double));

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m_copy[k] = m[i][j];
        k++;
     }
  }

  set_mat_identity(p, m_inv);

  //-----  Use LAPACK  -------
  if(0!=(k=clapack_dposv(CblasRowMajor, CblasUpper, p, p, m_copy, p, m_inv, p)))
  {
    fprintf(stderr,"Something was wrong with clapack_dposv [%d]\n",k);
     exit(1);
  }
  //--------------------------

  k=0;
  for(i=0;i<p;i++)
  {
     for(j=0;j<p;j++)
     {
        m[i][j] = m_inv[k];
        k++;
     }
  }  

  free(m_copy);
  free(m_inv);

  return;
}


//computes the log of the determinant of a symmetric positive definite matrix
double logdet(int p,double** m)
{
        //just take care of the 1x1 case
        if(1==p)
	{
	  return(log(m[0][0]));
	}

	int i,j;
	char jobvl = 'N';
	char jobvr = 'N';
	int lda = p;
	double wr[2*p];
	double wi[2*p];
	double vl[p][p];
	int ldvl = p*p;
	double vr[p][p];
	int ldvr = p*p;
	double work[p*p];
	int lwork = p*p;
	double a[p][p];
	int info = 1;
	
	for(i=0;i<p;i++)
	{
		for(j=0;j<p;j++)
		{
			a[i][j] = m[i][j];
		}
	}
	dgeev_(&jobvl,&jobvr,&p,(double*)a,&lda,(double*)wr,(double*)wi,(double*)vl, 
		  &ldvl,(double*)vr,&ldvr,(double*)work,&lwork,&info);

	if(0!=info)
	{
		printf("Smth wrong in the call of 'dgeev' error is [info = %d]\n",info);
		exit(1);
	}	   
	
	double logdet = 0;
	for(i=0;i<p;i++) logdet+=log(wr[i]);	
	return(logdet);
}

//this returns the columns specified by the set of indices A
//lenA is the number of indices of A
//data is an nxp matrix
//the returned submatrix is n x lenA
double** submatrix(int n,int p,double** data,int lenA,int* A)
{
  double** result = allocmatrix(n,lenA);

  int i,j;
  for(i=0;i<lenA;i++)
  {
    for(j=0;j<n;j++)
    {
      result[j][i] = data[j][A[i]-1];
    }
  }

  return(result);
}

//computes the marginal likelihood
double marglik(int n,int p,double** data,int lenA,int* A)
{
  
  	double** D_A = allocmatrix(n, lenA);
	double** response = allocmatrix(n,1);

	int i,j,k;

	subsetDataMatrix(data, D_A, lenA, n, A);

	// fill the response "matrix"
	for(i = 0; i < n; i++)
	{
		response[i][0] = data[i][0];
	}

	double** D_A_transpose = transposematrix(n, lenA, D_A);
	double** response_transpose = transposematrix(n, 1, response);

	double gam_coef;
	double l_det;
	double matprod;
	double inner_d = 0.0;
	double total;

	//get the funky log gamma piece
	gam_coef = lgamma(((double)n + 2.0 + (double)lenA)/2.0) - lgamma(((double)lenA + 2.0)/2.0);

	double** M_A = allocmatrix(lenA, lenA);
	matrixproduct(lenA, n, lenA, D_A_transpose, D_A, M_A);

	for(j = 0; j < lenA;j ++) // Add the identity matrix
	{
		M_A[j][j] += 1;
	}

	l_det = logdet(lenA, M_A);

	inverse(lenA, M_A); // invert MA. MA will be destroyed, but we don't need it anymore.

	// get the inner product of the response
	for(k = 0; k < n; k++)
	{
		inner_d += response[k][0]*response[k][0];
	}

	// matrices for intermediate calculations
	double** work1 = allocmatrix(1,lenA);
	double** work2 = allocmatrix(1,lenA);
	double** work3 = allocmatrix(1, n);
	double** work4 = allocmatrix(1,1);

	// incrememntally compute t(D_1)%*%(D_A)%*%inv(M_A)%*%t(D_A)%*%D_1
	matrixproduct(1, n, lenA, response_transpose, D_A, work1);
	matrixproduct(1, lenA, lenA, work1, M_A, work2);
	matrixproduct(1, lenA, n, work2, D_A_transpose, work3);
	matrixproduct(1, n, 1, work3, response, work4);

	// technically work4 is a 1x1 matrix, we need to extract its element.
	matprod = work4[0][0];

	total = gam_coef - .5*l_det - ((n + lenA + 2.0)/2.0)*log(1.0 + inner_d - matprod);

	// delete matrices from memory
	freematrix(n, D_A);
	freematrix(n, response);
	freematrix(lenA, D_A_transpose);
	freematrix(1, response_transpose);
	freematrix(lenA, M_A);
	freematrix(1, work1);
	freematrix(1, work2);
	freematrix(1, work3);
	freematrix(1, work4);

	return(total);
}
