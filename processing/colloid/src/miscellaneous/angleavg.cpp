#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

double * angleavg(int &nx, int &ny, double *data, int &nr)
{
	double xc=0.5*(double)nx;
	double yc=0.5*(double)ny;
	const int rmax=( nx < ny ) ? ( nx/2 ) : ( ny/2 );
	nr=rmax+1;

	int i, j;
	int xmin, xmax, ymin, ymax;
	if ( nx < ny)
	{
		xmin = 0; xmax = nx;
		ymin = (int)yc - rmax; ymax = (int)yc + rmax;
	}
	else
	{
		xmin = (int)xc - rmax; xmax = (int)xc + rmax;
		ymin = 0; ymax = ny;
	}
	
	double *avg=(double *)calloc(nr, sizeof(double));
	if ( avg==NULL )
		pERROR("double *avg initialized wrong!");
	
	double *count=(double *)calloc(nr, sizeof(double));
	if ( count==NULL )
		pERROR("double *count initialized wrong!");

	double r, dR;
	int R;
	for (i=xmin; i<xmax; i++)
	{
		for (j=ymin; j<ymax; j++)
		{
			r=sqrt((i-xc)*(i-xc)+(j-yc)*(j-yc));
			R=(int)r;
			dR=r-(double)R;
			if (R<nr)
			{
				avg[R]+=data[i+j*nx]*(1-dR);
				avg[R+1]+=data[i+j*nx]*dR;
				count[R]+=1-dR;
				count[R]+=dR;
			}
		}
	}

	for (i=0; i<nr; i++)
		if ( count[i] > 1.0e-6 )
			avg[i]/=count[i];

	return avg;
}
