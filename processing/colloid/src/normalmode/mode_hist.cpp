#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void mode_hist_w(int Gn, double *G, int mode, const char *file)
{
	const int total_p=Gn/2;
	if (mode>=Gn)
		pERROR("mode number exceed the total number of modes.");
	else if (mode<0)
		pERROR("mode is negative.");

	double *pointer=G+mode*Gn;
	// rescale the mode so that the result will be standard normal
	// e' = e*sqrt(2N)
	double scale=sqrt(Gn);
	int i;
	for(i=0; i<Gn; i++)
		*(pointer++) *= scale;
	pointer -= Gn;
	double xmin= *pointer, xmax= *(pointer++); 
	double ymin= *pointer, ymax= *(pointer++);
	for(i=1; i<total_p; i++)
	{
		if ( xmin > *pointer )
			xmin = *pointer;
		else if ( xmax < *pointer )
			xmax = *pointer;
		++pointer;
	
		if ( ymin > *pointer )
			ymin = *pointer;
		else if ( ymax < *pointer )
			ymax = *pointer;
		++pointer;
	}

	const int nbin=(int)sqrt(Gn);
	double xbinsize=(1+1.0e-6)*(xmax-xmin)/nbin;
	double ybinsize=(1+1.0e-6)*(ymax-ymin)/nbin;

	int *xhist=(int *)calloc(nbin, sizeof(int));
	if (xhist==NULL)
		pERROR("int *xhist initialized wrong!");
	int *yhist=(int *)calloc(nbin, sizeof(int));
	if (yhist==NULL)
		pERROR("int *yhist initialized wrong!");

	pointer=G+mode*Gn;
	for (i=0; i<total_p; i++)
	{
		++xhist[(int)((*(pointer++) - xmin)/xbinsize)];
		++yhist[(int)((*(pointer++) - ymin)/xbinsize)];
	}


	char suffix [20];
	sprintf(suffix, "_%04d.mh", mode);
	char *filename=getfilename(file, suffix);

	FILE *mhf=fopen(filename, "w");
	if (mhf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}
	
	fprintf(mhf, "# histogram of displacement in mode %d\n", mode);
	fprintf(mhf, "# e' = e*sqrt(2N)\n"
				 "# e_x e'_x probability hist e_y e'_y probability hist\n");
	for (i=0; i<nbin; i++)
	{
		fprintf(mhf, "%6.6f\t%6.6f\t%6.6f\t%d\t%6.6f\t%6.6f\t%6.6f\t%d\t\n",
				(xmin+(i+0.5)*xbinsize)/scale,
				xmin+(i+0.5)*xbinsize,
				(double)(xhist[i])/total_p/xbinsize,
				xhist[i],
				(ymin+(i+0.5)*ybinsize)/scale,
				ymin+(i+0.5)*ybinsize,
				(double)(yhist[i])/total_p/ybinsize,
				yhist[i]);
	}
	fclose(mhf);
	free(filename);
	free(xhist);
	free(yhist);
}


void mode_hist(colloid_base &ptid, int mode, 
		const char *file, bool Remove_Drift)
{
	int Gn;
	double *E, *G;
	if ( readev(Gn, &E, &G, file) != 0 )
	{
		normal_mode(ptid, file, Remove_Drift);
		if ( readev(Gn, &E, &G, file) != 0 )
			pERROR("Can not generate .ev file!");
	}
	free(E);

	mode_hist_w(Gn, G, mode, file);
	free(G);
}
