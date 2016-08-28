#include "fftw.h"

#include <cstdlib>

double* fft1dr(int& n, double *data)
{
	int i, m=n/2+1;
	fftw_complex* ft=(fftw_complex *)fftw_malloc(m*sizeof(fftw_complex));

	/* Create plans */
	fftw_plan plan=fftw_plan_dft_r2c_1d(n, data, ft, FFTW_ESTIMATE);

	/* Compute forward DFT */
	fftw_execute(plan);

	/* Free memory */
	fftw_destroy_plan(plan);

	double *dft=(double *)malloc(2*m*sizeof(double));
	double *pointer=dft;
	for (i=0; i<m; i++)
	{
		*(pointer++)=ft[i][0];
		*(pointer++)=ft[i][1];
	}
	fftw_free(ft);
	return dft;
}
