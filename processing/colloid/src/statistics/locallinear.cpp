#include "statistics.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/**
 * \function locallinear
 * \brief local linear regression
 *
 * \param n
 * 		data number
 * \param x
 * 		x coordinate, size=n
 * \param y
 * 		y coordinate, size=n
 * \param reg
 * 		regression result, size=n. Should be allocate memory before calling
 * 		the function. the corresponding x coordinate is x.
 * \param h
 *		bandwidth
 *		if GAUSSIAN kernel, h is optimized according to 
 *		http://en.wikipedia.org/wiki/Kernel_density_estimation
 * \param KERNEL
 * 		KERNEL to choose
 *
 * the algorithm is based on http://en.wikipedia.org/wiki/Kernel_smoother.
 */

void locallinear(int n, double *x, double *y, double *reg,
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
		pERROR("bandwidth should be positive.");

	double *BW=Malloc(double, 2*n);
	POINTER_NULL(BW);
	double *BWB=Malloc(double, 4);
	POINTER_NULL(BWB);
	double *BWy=Malloc(double, 2);
	POINTER_NULL(BWy);

	double xk, det;
	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)
		{
			xk=(x[i]-x[j])/h;
			xk=kernel_function[KERNEL](xk);
			BW[j]=xk;
			BW[n+j]=xk*x[j];
		}
		BWB[0]=pairwise(n, BW);
		BWB[1]=pairwise(n, BW+n);
		BWB[2]=BWB[1];
		for (j=0; j<n; j++)
			tmp[j]=BW[n+j]*x[j];
		BWB[3]=pairwise(n, tmp);

		det=BWB[0]*BWB[3]-BWB[1]*BWB[2];
		xk=BWB[0];
		BWB[0]=BWB[3]/det;
		BWB[1]=-BWB[2]/det;
		BWB[2]=BWB[1];
		BWB[3]=xk/det;    // (B^T)WB^{-1}

		for (j=0; j<n; j++)
			tmp[j]=BW[j]*y[j];
		BWy[0]=pairwise(n, tmp);
		for (j=0; j<n; j++)
			tmp[j]=BW[n+j]*y[j];
		BWy[1]=pairwise(n, tmp);

		reg[i]=(BWB[0]+x[i]*BWB[2])*BWy[0]
			   +(BWB[1]+x[i]*BWB[3])*BWy[1];
	}

	free(tmp);
	free(BW);
	free(BWB);
	free(BWy);
}
