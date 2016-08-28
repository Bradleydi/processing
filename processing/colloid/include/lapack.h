#ifndef LAPACK_H
#define LAPACK_H

void dsyev(double *, int , double *);

void dpotri(double *, int);
	

extern "C" {
	// A = alpha * xx^T + A
	void dsyr_(char *UPLO, int *N, double *ALPHA, double *X, int *INCX,
				double *A, int *LDA );
}


#endif /* LAPACK_H */
