#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

// here the algorithm is 1/2 * ptcl * ptcl * mode = O(N^3)

void mode_spcorr_all(const char *file, double omegaRange)
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
	
	// get all normal modes from .ev file
	char *evfilename=getfilename(file, ".ev");
	int Gn;
	double *E;
	double *G;
	if ( readev(Gn, &E, &G, evfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", 
				evfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");

	int i, j, k, m, total_m=Gn-2;
	// Gnuplot only support float
	float *omega=Malloc(float, Gn-2); POINTER_NULL(omega);
	for (i=0; i<Gn-2; i++)
		omega[i]=1.0/(float)sqrt(E[i+2]);
	free(E);
	if (omegaRange < 0.0 || omegaRange > omega[0] )
		omegaRange = omega[0];

	//==================================================
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

	float dx, dy;
	float w=xmax-xmin, h=ymax-ymin;
	const int nbin=(int)(sqrt(total_p));
	const float binsize=sqrt(w*w+h*h)/2.0/nbin;
	double *corr=Calloc(double, nbin*total_m); POINTER_NULL(corr);
	unsigned int *counter=Calloc(unsigned int, nbin); 
	POINTER_NULL(counter);

	pmp=mp; pG=G+2*Gn;
	float *pmp2; double *pG2, *tmp;
	for (i=0; i<total_p; i++)
	{
		x=*pmp; y=*(pmp+1);
		pmp2=pmp+2; pG2=pG+2; // pmp2 should be point to i+1, not i
		for (j=i+1; j<total_p; j++) // doesn't accout itself
		{
			dx = *pmp2 - x;
			dy = *(pmp2+1) - y;
			//d = sqrt(dx*dx + dy*dy);

			k=(int)(sqrt(dx*dx + dy*dy)/binsize);
			if (k<nbin)
			{
				tmp=corr+k*total_m;
				for (m=0; m<total_m; m++)
				{
					//pG=G+(m+2)*Gn+i*2;
					//pG2=G+(m+2)*Gn+j*2;
					//corr[k*total_m+m] +=  (*pG)*(*pG2) + (*(pG+1))*(*(pG2+1)) ;
					*(tmp++) +=  (*pG)*(*pG2) + (*(pG+1))*(*(pG2+1)) ;
					pG += Gn; pG2 += Gn;
				}
				++counter[k];
				pG -= Gn*total_m;
				pG2 -= Gn*total_m;
			}
			pmp2 += 2; pG2 += 2;
		}
		pmp+=2; pG+=2;
	}

	free(mp);
	free(G);

	//===================================================
	// print out the data
	char *filename = getfilename(file, ".ecr");
	FILE *outfile=fopen(filename, "wb");
	FILE_NULL(outfile, filename);

	//fprintf(outfile,"# eigenvector correlation vs. r\n"
	//				"# %d'th normal mode\n"
	//				"# eigenvalue = %f\n", n, E);
	
	// print as binary matrix format
	// Gnuplot only support binary matrix as single float
	// we put y as Omega, and x as r
	x=(float)total_m;
	if ( (int)fwrite(&x, sizeof(float), 1, outfile) != 1 )
	{
		fprintf(stderr, "# Error: writing to file '%s' failed.\n", filename);
		exit (1);
	}
	if ( (int)fwrite(omega, sizeof(float), total_m, outfile) != total_m )
	{
		fprintf(stderr, "# Error: writing to file '%s' failed.\n", filename);
		exit (1);
	}
	double norm=1.0/total_p; // total(e_i*e_i)=1
	for (i=0; i<nbin; i++)
	{
		if ( counter[i] != 0 )
		{
			x=(float)((i+0.5)*binsize);
			if ( (int)fwrite(&x, sizeof(float), 1, outfile) != 1 )
			{
				fprintf(stderr, "# Error: writing to file '%s' failed.\n",
						filename);
				exit (1);
			}
			for (m=0; m<total_m; m++)
				omega[m]=(float)((*(corr++))/counter[i]/norm);
			if ( (int)fwrite(omega, sizeof(float), total_m, outfile) != 
					total_m )
			{
				fprintf(stderr, "# Error: writing to file '%s' failed.\n", 
						filename);
				exit (1);
			}
		}
		else // a bug that I forgot the case of corr is defined for all nbin
			corr += total_m; // without this, free corr will lead error
	}
	fclose(outfile);
	corr -= total_m*nbin;
	free(corr);
	free(counter);
	free(omega);

	//=========================================================
	// plot by GNUplot, using pipe
	//FILE *pipe=popen("gnuplot -persist", "w");
	char *ecrfile = getfilename(file, "_ecr.plt");
	FILE *ecrf=fopen(ecrfile, "w");
	FILE_NULL(ecrf, ecrfile);
	free(ecrfile);

	fprintf(ecrf,  "reset\n"
				   "unset key\n"
				   "set terminal postscript eps color enhanced\n"
				   "set output ''\n\n"
				   "set xrange [0:%f]\n"
				   "set yrange [0:%f]\n"
				   "set cbrange [:1.0]\n"
				   "set pm3d map\n"
				   "set xlabel 'r'\n"
				   "set ylabel 'Omega'\n"
				   "#set logscale y\n"
				   "set palette rgbformulae 33,13,10\n"
				   ,nbin*binsize, omegaRange);

	fprintf(ecrf, "splot '%s' binary using 2:1:3 with image palette\n", 
			filename);
	fclose(ecrf);
	free(filename);
}
