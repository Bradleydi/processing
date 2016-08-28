#include "miscellaneous.h"

using namespace std;

double pairwise(int n, double *a)
{
	int i, m;
	double s;
	if (n<=2)
	{
		s=*a;
		for (i=1; i<n; i++)
			s+=*(a+i);
	}
	else
	{
		m=n/2;
		s=pairwise(m, a)+pairwise(n-m, a+m);
	}
	return s;
}
