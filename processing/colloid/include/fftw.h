#ifndef FFTW_H
#define FFTW_H

#include <fftw3.h>

double * fft1dr(int &, double *);

//double * fft2dr(int & /* row */, int & /* col */, double *);
void fft2dr(int row, int col, double *data, fftw_complex **pdft);

#endif /* FFTW_H */
