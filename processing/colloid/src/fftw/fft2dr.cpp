#include "fftw.h"

#include <cstdlib>

void fft2dr(int row, int col, double *data, fftw_complex **pdft)
{
	int i, m=row*(col/2+1);
	fftw_complex* ft=(fftw_complex *)fftw_malloc(m*sizeof(fftw_complex));
	*pdft = ft;
	/* Create plans */
	fftw_plan plan=fftw_plan_dft_r2c_2d(row, col, data, ft, FFTW_ESTIMATE);

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
}
