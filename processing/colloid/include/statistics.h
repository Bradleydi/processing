// version 1.0


#ifndef STATISTICS_H
#define STATISTICS_H

#include "colloid_base.h"
#include "miscellaneous.h"
#include <cmath>

struct statistics_base
{
	float mean;
	float sd;

	float min;
	float max;

	int nbin;
	float binsize;
	float * bin;  // center of bin
	int * histogram;
};

void statistics(statistics_base&, const colloid_base&, const int&, const int& =0);


//void pt_info(const colloid_base&, int =0, const bool& =false);
void pt_info(const char *file, int dim, bool show_histogram);

/*!
 * get the data information for a given tracked data.
 *  @param ptid
 *  	tracked data.
 *  @param dim
 *  	dimension of the data, default as ncol-2 (set to 0)
 *  @see data_info
 *  @see pt_info
 */
//void ptid_info(const colloid_base&, int =0, const bool& =false);
void ptid_info(const char *file, int dim, bool show_histogram);


/* ===================================================================
 * kernel functions
 * Those kernels are defined as
 * http://en.wikipedia.org/wiki/Uniform_kernel#Kernel_functions_in_common_use
 */
inline double Uniform(double x)
{
	if ( x >= -1 && x <= 1 )
		return 0.5;
	else
		return 0;
}

inline double Triangular(double x)
{
	if (x<0) x=-x;
	if ( x <= 1 )
		return 1-x;
	else
		return 0;
}

inline double Epanechnikov(double x)
{
	if ( x >= -1 && x <= 1 )
		return 0.75*(1-x*x);
	else
		return 0;
}

inline double Quartic(double x)
{
	if ( x >= -1 && x <= 1 )
		return 15.0*(1-x*x)*(1-x*x)/16.0;
	else
		return 0;
}

inline double Triweight(double x)
{
	if ( x >= -1 && x <= 1 )
		return 35.0*POW3(1-x*x)/32.0;
	else
		return 0;
}

inline double Tricube(double x)
{
	if (x<0) x=-x;
	if ( x <= 1 )
	{
		x=1-x*x*x;
		return 70.0*x*x*x/81.0;
	}
	else
		return 0;
}

#ifndef PI_M   // PI macro
#define PI_M 3.141592654
#endif
#ifndef inv_sqrt2pi
#define inv_sqrt2pi 0.3989422803
#endif

inline double Gaussian(double x)
{
	return inv_sqrt2pi*exp(-0.5*x*x);
}

inline double Cosine(double x)
{
	if ( x >= -1 && x <= 1 )
		return PI_M*cos(x*PI_M/2.0)/4.0;
	else
		return 0;
}

typedef enum { UNIFORM, TRIANGULAR, EPANECHNIKOV, QUARTIC,
				TRIWEIGHT, TRICUBE, GAUSSIAN, COSINE} kernel_t;

typedef enum {NOT_ORDER, ORDER} order_t;

// Compiler will typically not inline the function call as it could do 
// anywhere else. So using function pointers may add up to be significantly 
// slower than using regular function calls, and be avoided as a way 
// to gain performance.
// But in my code, use function pointer will be nice.
// However, the following is an assignment, which should not be in header file,
// or will cause multiple definition error when compiling.
//double (*kernel_function[])(double)={ Uniform, Triangular, 
//	Epanechnikov, Quartic, Triweight, Tricube, Gaussian, Cosine };

/* ===================================================================
 * end of kernel functions */

/* kernel density estimation */
void kdensity(int n, double *x, int nbin, double *density,
		        double h /*bandwidth*/, kernel_t KERNEL);

/* local linear regression */
void locallinear(int n, double *x, double *y, double *reg,
		double h /*bandwidth*/, kernel_t KERNEL);

typedef enum { KDE, LLR } unfolding_t;

double *unfolding(int n, double *x, double h /*bandwidth*/, 
		unfolding_t UNFOLDING);

/* Windsorize data */
void windsorize(int N, double *y, double margin,
		double *pBt, double *pUp);
void windsorize_kernel(int N, double *y, double Btmargin, double Upmargin,
		double *pBt, double *pUp);
void windsorize_wcp(int N, double *y, double *copy, double margin,
		double *pBt, double *pUp);
void windsorize_wcp_kernel(int N, double *y, double *copy, 
		double Btmargin, double Upmargin,
		double *pBt, double *pUp);
#endif /* STATISTICS_H */
