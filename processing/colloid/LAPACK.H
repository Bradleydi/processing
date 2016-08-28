#ifndef DGEEV_H
#define DGEEV_H

#include "miscellaneous.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//extern "C" {
	void dgeev_(char *JOBVL, char * JOBVR, int *N, double *A,
				int *LDA, double *WR, double *WI, double *VL, int *LDVL,
				double *VR, int *LDVR, double *WORK, int *LWORK, int *INFO);
	
	void dsyev_(char *JOBZ, char *UPLO, int *N, double *A,
				int *LDA, double *W, double *WORK, int *LWORK,
				int *INFO);

	void dgetrf_(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);

	// the inverse of a real symmetric positive definite
	// matrix A using the Cholesky factorization
	void dpotrf_(char *UPLO, int *N, double *A, int *LDA, int *INFO );
	void dpotri_(char *UPLO, int *N, double *A, int *LDA, int *INFO );

	//=======================================================================
	// The multiplication among matrices and vectors should be
	// implimented as fast as possible, since those will be called many times.

	// general matrix-matrix multiplication
	void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
				double *ALPHA, double *A, int *LDA, double *B, int *LDB,
				double *BETA, double *C, int *LDC);

	// symmetric matrix-matrix multiplication
	void dsymm_(char *SIDE, char *UPLO, int *M, int *N, double *ALPHA,
			    double *A, int *LDA, double *B, int *LDB, 
				double *BETA, double *C, int *LDC);

	// general matrix-vector multiplication
	// http://www.netlib.org/blas/dgemv.f
	void dgemv_(char *TRANS, int *M, int *N, double *ALPHA,
				double *A, int *LDA, double *X, int *INCX,
				double *BETA, double *Y, int *INCY);

	// symmetric matrix-vector multiplication
	// http://www.netlib.org/blas/dsymv.f
	void dsymv_(char *UPLO, int *N, double *ALPHA, double *A,
				int *LDA, double *X, int *INCX, 
				double *BETA, double *Y, int *INCY);
	//=======================================================================

	// add of two vectors  y = alpha * x + y
	void daxpy_(int *N, double *ALPHA, double *X, int *INCX, 
			double *Y, int *INCY);

	// A = alpha * xx^T + A
	void dsyr_(char *UPLO, int *N, double *ALPHA, double *X, int *INCX,
				double *A, int *LDA );

	// A = alpha * xy^T + A
	void dger_(int *M, int *N, double *ALPHA, 
				double *X, int *INCX, double *Y, int *INCY, 
				double *A, int *LDA);

	double ddot_(int *N, double *DX, int *INCX, double *DY, int *INCY);

//}

void dgeev(double *A, int N, double *ER, double *EI, double *VR);
//{
//	// ER, EI, V should be allocated memory before calling
//	printf("# Calling dgeev...\n");
//	char jobvl='N', jobvr='N';
//	if (VR!=NULL)
//		jobvr='V';
//	int lda=N, lwork=-1, info;
//
//	double workquery;
//
//	double *VL;
//
//	int ldvl=1, ldvr=1;
//
//	/* get the optimized length of workspace */
//	dgeev_(&jobvl, &jobvr, &N, A, &lda, ER, EI, VL, &ldvl,
//			VR, &ldvr, &workquery, &lwork, &info);
//
//	lwork=(int)workquery;
//
//	double *work=(double *)malloc(lwork*sizeof(double));
//	if (work==NULL)
//	{
//		fprintf(stderr, "# Error: double *work, no enough memory!\n");
//		exit (1);
//	}
//
//	dgeev_(&jobvl, &jobvr, &N, A, &lda, ER, EI, VL, &ldvl,
//			VR, &ldvr, work, &lwork, &info);
//	if (info !=0 )
//	{
//		fprintf(stderr, "# lapack: failure with error %d\n", info);
//		exit (1);
//	}
//
//	free(work);
//}

//inline unsigned char is_zero(double a);
inline BOOL is_zero(double a);
//{
//	if (a<1.0E-6 && a>-1.0E-6)
//		return 1u;
//	else
//		return 0u;
//}

//void dgetrf(int M, int N, double *A, int LDA, int *IPIV, int INFO);
BOOL singularity(int N, double *A);

//int main(int argc, char *argv[])
//{
//	double k, beta1, beta2;
//	if (argc!=4)
//	{
//		printf("dgeev k beta1 beta2");
//		exit(1);
//	}
//	else
//	{
//		k=atof(argv[1]);
//		beta1=atof(argv[2]);
//		beta2=atof(argv[3]);
//	}
//	
//	int n=6;
//	double gamma=1.50, g=-gamma;
//	//double *a=(double *)calloc(n*n, sizeof(double));
//
//	k=-k;
//
//	double a[]={
//		0.00,  0.00, 0.00,  1.00, 0.00, 0.00,
//		0.00,  0.00, 0.00,  0.00, 1.00, 0.00,
//		0.00,  0.00, 0.00,  0.00, 0.00, 1.00,
//		k,     -k,   0.00,  g,       g*beta1,  g*beta2,
//		-k,    2*k,  -k,    g*beta1, g,        g*beta1,
//		0.00, -k,    k,     g*beta2, g*beta1,  g};
//
//
//	/*=================================================================
//	 * calculate the eigenvalues
//	 *=================================================================*/
//
//	/* REAL part */
//	double *ER=(double *)malloc(n*sizeof(double));
//	/* IMAGINARY part */
//	double *EI=(double *)malloc(n*sizeof(double));
//
//	dgeev(a, n, ER, EI);
//
//	printf("k     = %f\t"
//		   "gamma = %f\n"
//		   "beta1 = %f \t"
//		   "beta2 = %f \n"
//		   "intrinsic freq. = %f\n",
//			-k, gamma, beta1, beta2, sqrt(-k));
//
//	printf("Under damping gamma, critical freq. = %f\n",
//			gamma/2);
//
//	int i;
//	for (i=0; i<n; i++)
//		printf("%6.6f\t+ i *\t%6.6f\n", ER[i], EI[i]);
//	
//
//	free(ER);
//	free(EI);
//}

// multiplication of two square matrix: C=AB
void dgemm(int N, double *A, double *B, double *C);

// multiplication of square matrix and vector: y=Ax
void dgemv(int N, double *A, double *x, double *y);

#endif /* DGEEV_H */
