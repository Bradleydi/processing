#include "lapack.h"

#include <cstdio>
#include <cstdlib>


extern "C" {
	void dsyev_(char *JOBZ, char *UPLO, int *N, double *A,
			int *LDA, double *W, double *WORK, int *LWORK,
			int *INFO);
}

/*! E should be allocated
 */

void dsyev(double *A, int N, double *E)
{
	char jobz='V', uplo='L';
	int lda=N, lwork=-1, info;
	
	double workquery;

	/* get the optimized length of WORK */
	dsyev_(&jobz, &uplo, &N, A, &lda, E, &workquery, &lwork, &info);

	lwork=(int)workquery;
	
	double *work=(double *)malloc(lwork*sizeof(double));
	if (work==NULL)
	{
		fprintf(stderr, "# Error: double *work initialized wrong!\n");
		exit (1);
	}

	dsyev_(&jobz, &uplo, &N, A, &lda, E, work, &lwork, &info);
	if (info != 0)
	{
		fprintf(stderr, "# lapack/dsyev: failure with error %d\n", info);
		exit (1);
	}
	free(work);
}
