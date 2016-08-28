#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void Kr(const char *file)
{
	// read K
	int Kn;
	double *K;
	readK(Kn, &K, file);

	// read mean position
	int total_p;
	float *mp_ptr;
	readmp(total_p, &mp_ptr, file);
	
	if ( Kn != 2*total_p)
		pERROR("read files error!");



	int i, j, k;
	const int dim=2;
	

	float xmin=mp_ptr[0], xmax=mp_ptr[0];
	float ymin=mp_ptr[1], ymax=mp_ptr[1];
	float *pointer=mp_ptr+dim;
	for (i=1; i<total_p; i++)
	{
		if ( xmin > *pointer)
			xmin = *pointer;
		else if ( xmax < *pointer)
			xmax = *pointer;
		++pointer;

		if ( ymin > *pointer)
			ymin = *pointer;
		else if ( ymax < *pointer)
			ymax = *pointer;
		++pointer;
	}

	double dx, dy, d;

	const int nbin=100;
	const double binsize=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin))/2/nbin;

	double *Krxx=(double *)calloc(nbin, sizeof(double));
	double *Krxy=(double *)calloc(nbin, sizeof(double));
	double *Kryx=(double *)calloc(nbin, sizeof(double));
	double *Kryy=(double *)calloc(nbin, sizeof(double));
	POINTER_NULL(Krxx);
	POINTER_NULL(Krxy);
	POINTER_NULL(Kryx);
	POINTER_NULL(Kryy);

	unsigned int *counter=(unsigned int *)calloc(nbin, sizeof(unsigned int));
	POINTER_NULL(counter);


	for (i=0; i<total_p; i++)
	{
		for (j=i; j<total_p; j++)
		{
			dx=(double)(mp_ptr[j*dim]-mp_ptr[i*dim]);
			dy=(double)(mp_ptr[j*dim+1]-mp_ptr[i*dim+1]);
			d=sqrt(dx*dx+dy*dy);

			k=(int)(d/binsize);

			if (k<nbin)
			{
				Krxx[k]+=K[2*i*Kn+2*j];
				Krxy[k]+=K[2*i*Kn+2*j+1];
				Kryx[k]+=K[2*i*Kn+Kn+2*j];
				Kryy[k]+=K[2*i*Kn+Kn+2*j+1];
				++counter[k];
			}
		}
	}
	free(mp_ptr);
	
	char *filename=getfilename(file, ".Kr");

	FILE *Krf=fopen(filename, "w");
	FILE_NULL(Krf, filename);

	fprintf(Krf, "# <K(r'+r)> vs. distance\n");
	fprintf(Krf, "# r Kxx Kxy Kyx Kyy total_t\n");

	for (i=0; i<nbin; i++)
	{
		if ( counter[i] > 100 )
		{
			Krxx[i]/=counter[i];
			Krxy[i]/=counter[i];
			Kryx[i]/=counter[i];
			Kryy[i]/=counter[i];
			fprintf(Krf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%d\n",
					(i+0.5)*binsize, Krxx[i], Krxy[i],
					Kryx[i], Kryy[i], counter[i]);
		}
	}
	fclose(Krf);

	free(Krxx);
	free(Krxy);
	free(Kryx);
	free(Kryy);
	free(counter);
	free(filename);
}
