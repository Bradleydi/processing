#include "statistics.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/**
 * \function kdensity
 * \brief kernel density estimation unfolding with Gaussian kernel.
 *
 * \param n
 * 		data number
 * \param x
 * 		point process, {x_i} , size=n, should be uniqued and sorted first.
 * 		The sort should be nondecreasing.
 * 		A program sort(int n, double *) is used for such a sort,
 * 		and 
 * 			int N;
 * 			for (i=1; i<n; i++)
 * 				if ( x[i] > x[i-1] + 1.0e-7 )
 * 					x[N++] = x[i];
 * 			n=N;
 * 		is a unique procedure.
 * \param h
 *		bandwidth
 *		if GAUSSIAN kernel, h is optimized according to 
 *		http://en.wikipedia.org/wiki/Kernel_density_estimation
 * \return X
 * 		unfolded levels, size=n.
 *
 * the algorithm is based on 
 * 		http://en.wikipedia.org/wiki/Kernel_density_estimation.
 */

double *unfolding(int n, double *x, double h /*bandwidth*/, 
		unfolding_t UNFOLDING)
{
	int i,j;
	// check whether x is unique and sorted
	for (i=1; i<n; i++)
		if ( x[i] < x[i-1] + 1.0e-7 )
			pERROR("input array \"x\" is not sorted or uniqued,\n"
				   "A program sort(int n, double *) is used for such a sort,\n"
				   "and\n"
				   "    int N;\n"
				   "    for (i=1; i<n; i++)\n"
				   "        if ( x[i] > x[i-1] + 1.0e-7 )\n"
				   "            x[N++] = x[i];\n"
				   "    n=N;\n"
				   "is a unique procedure.");
	
	double *tmp=Malloc(double, n); POINTER_NULL(tmp);
	if ( h <= 0 )
	{
		double mean=pairwise(n, x)/n;
		for (i=0; i<n; i++)
			tmp[i]=(x[i]-mean)*(x[i]-mean);
		double sd=sqrt(pairwise(n, tmp)/n);
		h=sd/pow(0.75*n, 0.2);
		printf("# bandwidth is not positive, reset it to %f\n", h);
	}

	// unfolding
	double *X=Malloc(double, n); POINTER_NULL(X);
	double *a=Malloc(double, n); POINTER_NULL(a);

	if ( UNFOLDING == KDE )
	{
		h*=sqrt(2.0);
		for (i=0; i<n; i++)
			a[i]=x[i]/h;

		for (i=0; i<n; i++)
		{
			for (j=0; j<n; j++)
				tmp[j]=(1+erf(a[i]-a[j]))/2;
			X[i]=pairwise(n, tmp);
		}
	}
	else
	{
		for (i=0; i<n; i++)
			a[i]=(double)(i+1);
		locallinear(n, x, a, X, h, TRICUBE);
	}
	free(a);
	free(tmp);
	return X;
}
