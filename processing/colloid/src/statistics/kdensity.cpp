#include "statistics.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/**
 * \function kdensity
 * \brief kernel density estimation
 *
 * \param n
 * 		data number
 * \param x
 * 		x coordinate, size=n
 * \param nbin
 * 		bin number, if set to 0, nbin=sqrt(n)
 * \param density
 * 		estimated density, size=2*nbin. Should be allocate memory before calling
 * 		the function. data is saved as [x0, d0, x1, d1, x2, d2, ...]
 * \param h
 *		bandwidth
 *		if GAUSSIAN kernel, h is optimized according to 
 *		http://en.wikipedia.org/wiki/Kernel_density_estimation
 * \param KERNEL
 * 		KERNEL to choose
 *
 * the algorithm is based on 
 * 		http://en.wikipedia.org/wiki/Kernel_density_estimation.
 */

void kdensity(int n, double *x, int nbin, double *density,
		double h /*bandwidth*/, kernel_t KERNEL)
{
	double (*kernel_function[])(double)={ Uniform, Triangular, 
		Epanechnikov, Quartic, Triweight, Tricube, Gaussian, Cosine };

	double *tmp=Malloc(double, n);
	POINTER_NULL(tmp);

	int i,j;
	if ( KERNEL==GAUSSIAN && h <= 0 )
	{
		double mean=pairwise(n, x)/n;
		for (i=0; i<n; i++)
			tmp[i]=(x[i]-mean)*(x[i]-mean);
		double sd=sqrt(pairwise(n, tmp)/n);
		h=sd/pow(0.75*n, 0.2);
	}

	if (h<=0)
		pERROR("bandwidth should be positive.\n");

	double min=x[0], max=x[0];
	for (i=1; i<n; i++)
	{
		if ( min > x[i] )
			min = x[i];
		else if ( max < x[i] )
			max = x[i];
	}

	if (nbin<=0)
		nbin=(int)(sqrt(n));
	const double binsize=(max-min)/nbin;

	double xi, xk;

	double C=1.0/h/n;
	for (i=0; i<nbin; i++)
	{
		xi=min+(i+0.5)*binsize;
		for (j=0; j<n; j++)
		{
			xk=(xi-x[j])/h;
			tmp[j]=kernel_function[KERNEL](xk);
		}
		density[2*i]=xi;
		density[2*i+1]=C*pairwise(n, tmp);
	}

	free(tmp);
}
