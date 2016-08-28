#include "normalmode.h"

#include "colloid_base.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
void plot_mode(int& n, const char *file)
{
	// get mean position from .mp file
	int total_p;
	float *mp;
	const char mpsubfix[]=".mp";
	char *mpfilename=getfilename(file, mpsubfix);
	if ( readmp(total_p, &mp, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file \"%s\"\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);
	
	// get n'th normal mode from .ev file
	const char evsubfix[]=".ev";
	char *evfilename=getfilename(file, evsubfix);
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


	//===================================================
	// print out the data
	char modesubfix [15];
	sprintf(modesubfix, "_%d.mode", n);
	char *modefilename = getfilename(file, modesubfix);
	FILE *modefile=fopen(modefilename, "w");
	if (modefile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", modefilename);
		exit (1);
	}

	fprintf(modefile, "# %d'th normal mode\n"
			"# eigenvalue = %f\n", n, E);


	int i;
	float *pmp=mp;
	double *pG=G;
	double ratio=sqrt((double)(Gn))/3.0;
	float x=*pmp, y=*(pmp+1), xmin=x, xmax=xmin, ymin=y, ymax=ymin;
	for (i=0; i<total_p; i++)
	{
		// min & max
		x=*pmp; y=*(pmp+1);
		if ( x < xmin )
			xmin = x;
		else if ( x > xmax )
			xmax = x;
		if ( y < ymin )
			ymin = y;
		else if ( y > ymax )
			ymax = y;

		//fprintf(modefile, "%f %f %f %f\n", *pmp, *(pmp+1), 
		fprintf(modefile, "%f %f %f %f\n", x, y, 
				(*pG)*ratio, (*(pG+1))*ratio);
		pmp+=2;
		pG+=2;
	}
	fclose(modefile);

	free(mp);
	free(G);


	//=========================================================
	// plot by GNUplot, using pipe

	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}
	char NULLsubfix [2];
	NULLsubfix[0]='\0';
	char *NULLfilename = getfilename(file, NULLsubfix);

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set size ratio -1\n"
				   //"set autoscale\n"
				   "set xrange [%f:%f]\n"
				   "set yrange [%f:%f]\n"
				   "set title \"%d'th normal mode\t\tomega=%f\"\n"
				   ,xmin-1.0, xmax+1.0, ymin-1.0, ymax+1.0,
				   n, E);


	fprintf(pipe, "plot '%s' using 1:2:3:4 with vectors lc rgb 'black' lw 1.5\n", modefilename);

	fclose(pipe);

	free(modefilename);
	free(NULLfilename);
}


void plot_mode_eps(int& n, const char *file)
{
	// get mean position from .mp file
	int total_p;
	float *mp;
	const char mpsubfix[]=".mp";
	char *mpfilename=getfilename(file, mpsubfix);
	if ( readmp(total_p, &mp, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file \"%s\"\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);
	
	// get n'th normal mode from .ev file
	const char evsubfix[]=".ev";
	char *evfilename=getfilename(file, evsubfix);
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
	E=1.0/sqrt(E);

	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");


	//===================================================
	// print out the data
	char modesubfix [15];
	sprintf(modesubfix, "_%d.mode", n);
	char *modefilename = getfilename(file, modesubfix);
	FILE *modefile=fopen(modefilename, "w");
	if (modefile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", modefilename);
		exit (1);
	}

	fprintf(modefile, "# %d'th normal mode\n"
			"# eigenvalue = %f\n", n, E);


	int i;
	float *pmp=mp;
	double *pG=G;
	double ratio=sqrt((double)(Gn))/3.0;
	float x=*pmp, y=*(pmp+1), xmin=x, xmax=xmin, ymin=y, ymax=ymin;
	for (i=0; i<total_p; i++)
	{
		// min & max
		x=*pmp; y=*(pmp+1);
		if ( x < xmin )
			xmin = x;
		else if ( x > xmax )
			xmax = x;
		if ( y < ymin )
			ymin = y;
		else if ( y > ymax )
			ymax = y;

		//fprintf(modefile, "%f %f %f %f\n", *pmp, *(pmp+1), 
		fprintf(modefile, "%f %f %f %f\n", x, y, 
				(*pG)*ratio, (*(pG+1))*ratio);
		pmp+=2;
		pG+=2;
	}
	fclose(modefile);

	free(mp);
	free(G);


	//=========================================================
	// plot by GNUplot, using pipe
	char epsTLsubfix [17];
	sprintf(epsTLsubfix, "_mode_%d.eps", n);
	char *eps = getfilename(file, epsTLsubfix);


	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}
	//char NULLsubfix [2];
	//NULLsubfix[0]='\0';
	//char *NULLfilename = getfilename(file, NULLsubfix);

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set terminal postscript eps color enhanced\n"
				   "set output \"%s\"\n"
				   "set size ratio -1\n"
				   //"set autoscale\n"
				   "set xrange [%f:%f]\n"
				   "set yrange [%f:%f]\n"
				   "set title \"%d'th normal mode\t\tomega=%f\"\n"
				   "unset key\n"
				   "unset xtics\n"
				   "unset ytics\n"
				   "set border 0\n"
				   "#set xlabel 'x/a_0'\n"
				   "#set ylabel 'y/a_0'\n", eps, 
				   xmin-1.0, xmax+1.0, ymin-1.0, ymax+1.0,
				   n, E);


	fprintf(pipe, "plot '%s' using 1:2:3:4 with vectors lc rgb 'black' lw 1.5\n", modefilename);

	fclose(pipe);

	free(modefilename);
	//free(NULLfilename);
	free(eps);
}
