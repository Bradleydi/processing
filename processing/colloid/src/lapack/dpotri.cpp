#include "lapack.h"

#include <cstdio>
#include <cstdlib>


// http://www.netlib.org/lapack/double/dpotri.f
/*  DPOTRI computes the inverse of a real symmetric positive definite
 *  matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
 *  computed by DPOTRF.
 */
extern "C" {
	void dpotrf_(char *UPLO, int *N, double *A,
			int *LDA, int *INFO);
	void dpotri_(char *UPLO, int *N, double *A,
			int *LDA, int *INFO);
}


void dpotri(double *A, int N)
{
	char uplo='U';
	int lda=N, info;

	dpotrf_(&uplo, &N, A, &lda, &info);
	if (info != 0)
	{
		fprintf(stderr, "# lapack/dpotrf: failure with error %d\n", info);
		exit (1);
	}

	dpotri_(&uplo, &N, A, &lda, &info);
	if (info != 0)
	{
		fprintf(stderr, "# lapack/dpotri: failure with error %d\n", info);
		exit (1);
	}
	for (lda=0; lda<N; lda++)
		for (info=lda+1; info<N; info++)
			A[lda*N+info]=A[lda+N*info];
}
