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
	int info;
	
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

// Computes the log of the marginal likelihoood for normal linear regression with inverse gamma prior.
double marglik(int n,int p,double** data,int lenA,int* A)
{
	double lml;
	int i,j,k;
	
	// Create matrix just containing response variable
	double** d1 = allocmatrix(n,1);
	
	for(j=0;j<n;j++)
	{
			d1[j][0] = data[j][0];
	}

	printmatrix("d1.mat", n, 1, d1);

	// Create matrix just containing variables specified in A
	double** dA = allocmatrix(n, lenA);

	for(i=0;i<lenA;i++)
	{
		k = A[i] - 1; // C indexes at 0

		for(j=0;j<n;j++)
		{
			dA[j][i] = data[j][k];
		}
	}

	printmatrix("dA.mat", n, lenA, dA);

	// Calculate MA
	double** tdA = transposematrix(n,lenA, dA);
	printmatrix("tdA.mat", lenA, n, tdA);

	double** mA = allocmatrix(lenA, lenA);
	matrixproduct(lenA, n, lenA, tdA, dA, mA);
	
	for(i=0;i<lenA;i++)
	{
		for(j=0;j<lenA;j++)
		{
			if(i==j)
			{
				mA[i][j] += 1;
			}


		}
	}
	printmatrix("mA.mat", lenA, lenA, mA);

	// Calculate log marginal likelihood following the equation given
	double** td1 = transposematrix(n, 1, d1);

	double** td1_d1 = allocmatrix(1, 1);
	matrixproduct(1, n, 1, td1, d1, td1_d1);
	printmatrix("td1_d1.mat", 1, 1, td1_d1);

	double** td1_dA = allocmatrix(1, lenA);
	matrixproduct(1, n, lenA, td1, dA, td1_dA);
	printmatrix("td1_dA.mat", 1, lenA, td1_dA);

	double** mAInverse = allocmatrix(lenA,lenA);
    copymatrix(lenA,lenA,mA,mAInverse);
	inverse(lenA,mAInverse);
	printmatrix("mAInverse.mat", lenA, lenA, mAInverse);

	double** tdA_d1 = allocmatrix(lenA, 1);
	matrixproduct(lenA, n, 1, tdA, d1, tdA_d1);
	printmatrix("tdA_d1.mat", lenA, 1, tdA_d1);

	double** td1_dA_mAInverse = allocmatrix(1,lenA);
	matrixproduct(1, lenA, lenA, td1_dA, mAInverse, td1_dA_mAInverse);
	printmatrix("td1_dA_mAInverse.mat", 1, lenA, td1_dA_mAInverse);

	double** td1_dA_mAInverse_tdA_d1 = allocmatrix(1,1);
	matrixproduct(1, lenA, 1, td1_dA_mAInverse, tdA_d1, td1_dA_mAInverse_tdA_d1);
	printmatrix("td1_dA_mAInverse_tdA_d1.mat", 1, 1, td1_dA_mAInverse_tdA_d1);

	lml = lgamma((n + lenA + 2.0)/2.0) - lgamma((lenA + 2.0)/2.0) - (1.0/2.0)*logdet(lenA, mA) - ((n + lenA + 2.0)/2.0)*log(1.0 + **td1_d1 - **td1_dA_mAInverse_tdA_d1);

	// free(d1);
	// free(dA);
	// free(tdA);
	// free(mA);
	// free(td1);
	// free(td1_d1);
	// free(mAInverse);
	// free(tdA_d1);
	// free(td1_dA_mAInverse);
	// free(td1_dA_mAInverse_tdA_d1);
	freematrix(n,d1);
	freematrix(n,dA);
	freematrix(lenA, tdA);
	freematrix(lenA, mA);
	freematrix(1, td1);
	freematrix(1, td1_d1);
	freematrix(lenA, mAInverse);
	freematrix(lenA, tdA_d1);
	freematrix(1, td1_dA_mAInverse);
	freematrix(1, td1_dA_mAInverse_tdA_d1);

	return(lml);

}


