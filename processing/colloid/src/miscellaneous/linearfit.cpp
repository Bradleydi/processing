#include "miscellaneous.h"

using namespace std;

// x=[0, 1, 2, 3, 4, 5,...]
// y=ax+b
// y has been changed after calling the function
void linearfit(int n, double *y, double &a, double &b)
{
	double avx=(double)(n-1)/2.0;
	double ssx=(double)((n-1)*(2*n-1))/6.0-avx*avx;

	double avy=pairwise(n, y)/n;
	int i;
	for (i=0; i<n; i++)
		y[i]=(y[i]-avy)*(i-avx);
	double cxy=pairwise(n, y)/n;
	
	a=cxy/ssx;
	b=avy-a*avx;
}
