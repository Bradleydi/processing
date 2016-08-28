#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void mode_spcorr(int n, const char *file)
{
	// get mean position from .mp file
	int total_p;
	float *mp;
	char *mpfilename=getfilename(file, ".mp");
	if ( readmp(total_p, &mp, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file \"%s\"\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);
	
	// get n'th normal mode from .ev file
	char *evfilename=getfilename(file, ".ev");
	int Gn;
	double E;
	double *G;
	if ( readmode(n, Gn, E, &G, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	E=1/sqrt(E);

	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");


	//==================================================
	int i, j, k;
	float *pmp=mp;
	double *pG=G;
	float x=*pmp, y=*(pmp+1), xmin=x, xmax=xmin, ymin=y, ymax=ymin;
	for (i=0; i<total_p; i++)
	{
		// min & max
		x=*pmp; y=*(pmp+1);
		if ( x < xmin ) xmin = x;
		else if ( x > xmax ) xmax = x;

		if ( y < ymin ) ymin = y;
		else if ( y > ymax ) ymax = y;

		pmp+=2;
	}

	float dx, dy, d;
	float w=xmax-xmin, h=ymax-ymin;
	const int nbin=(int)(sqrt(total_p));
	const float binsize=sqrt(w*w+h*h)/2.0/nbin;
	double *corr=Calloc(double, nbin); POINTER_NULL(corr);
	unsigned int *counter=Calloc(unsigned int, nbin); POINTER_NULL(counter);

	pmp=mp; pG=G;
	float *pmp2; double *pG2;
	for (i=0; i<total_p; i++)
	{
		x=*pmp; y=*(pmp+1);
		//pmp2=mp+i*2; pG2=G+i*2; a bug here
		pmp2=mp+i*2+2; pG2=pG+2;
		for (j=i+1; j<total_p; j++) // doesn't accout itself
		{
			dx = *pmp2 - x;
			dy = *(pmp2+1) - y;
			d = sqrt(dx*dx + dy*dy);

			k=(int)(d/binsize);
			if (k<nbin)
			{
				corr[k] +=  (*pG)*(*pG2) + (*(pG+1))*(*(pG2+1)) ;
				++counter[k];
			}
			pmp2 += 2; pG2 += 2;
		}
		pmp+=2; pG+=2;
	}
	//===============================================
	/*
	float dx, dy, d, area;
	float w=xmax-xmin, h=ymax-ymin;
	float density=total_p/w/h;
	const int nbin=sqrt(total_p)/2.0;
	const float binsize=sqrt(w*w+h*h)/2.0/nbin;
	float rmax2=nbin*binsize-0.01; rmax2*=rmax2;
	double corr=Malloc(double, nbin); POINTER_NULL(corr);
	pmp=mp; pG=G;
	float *pmp2, *pG2;
	for (i=0; i<total_p; i++)
	{
		x=*pmp; y=*(pmp+1);
		pmp2=mp+i*2; pG2=G+i*2;
		for (j=i; j<total_p; j++)
		{
			dx = *pmp2 - x;
			d = dx*dx;
			if ( d < rmax2 )
			{
				dy = *(pmp2+1) - y;
				d += dy*dy;
				if ( d < rmax2 )
				{
					if (dx<0) dx=-dx;
					if (dy<0) dy=-dy;
					area = (w-dx)*(h-dy);

					corr[(int)(sqrt(d)/binsize)] += 
						( (*pG)*(*pG2) + (*(pG+1))*(*(pG2+1)) )/density/area;
				}
			}
			pmp2 += 2;
			pG2 += 2;
		}
		pmp+=2;
		pG+=2;
	}
	*/

	free(mp);
	free(G);

	//===================================================
	// print out the data
	char subfix [15];
	sprintf(subfix, "_%d.ecr", n);
	char *filename = getfilename(file, subfix);
	FILE *outfile=fopen(filename, "w");
	FILE_NULL(outfile, filename);

	fprintf(outfile,"# eigenvector correlation vs. r\n"
					"# %d'th normal mode\n"
					"# eigenvalue = %f\n", n, E);
	double norm=1.0/total_p; // total(e_i*e_i)=1
	for (i=0; i<nbin; i++)
		if ( counter[i] != 0 )
			fprintf(outfile, "%f %f %d\n", (i+0.5)*binsize,
							 corr[i]/counter[i]/norm, counter[i]);
	fclose(outfile);
	free(corr);
	free(counter);

	//=========================================================
	// plot by GNUplot, using pipe

	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set xrange [%f:%f]\n"
				   "set yrange [:%f]\n"
				   "set title \"%d'th normal mode\t\tomega=%f\"\n"
				   ,0.0, nbin*binsize, 1.0, n, E);

	fprintf(pipe, "plot '%s' using 1:2 with lp lw 1.5, 0\n", filename);

	fclose(pipe);

	free(filename);
}
