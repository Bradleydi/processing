#include "normalmode.h"

#include "colloid_base.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
void softspot(const char *file, int Nmode, int Nsmall)
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
	double *E;
	double *G;
	if ( readev(Gn, &E, &G, evfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", evfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);

	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");

	int i,j;
	//for (i=2; i<Gn; i++)
	//	E[i]=1/sqrt(E[i]);
	free(E);

	//===================================================
	// find soft particles
	double *pointer, e0;
	double *x2=Malloc(double, total_p); POINTER_NULL(x2);
	double *tmp=Malloc(double, total_p); POINTER_NULL(tmp);
	double *sortwork=Malloc(double, total_p); POINTER_NULL(sortwork);
	unsigned char *softspot=Calloc(unsigned char, total_p);
	POINTER_NULL(softspot);
	for (j=Gn-Nmode; j<Gn; j++)
	{
		pointer=G+j*Gn;
		for (i=0; i<total_p; i++)
		{
			e0=*(pointer++);
			x2[i]=e0*e0;
			e0=*(pointer++);
			x2[i]+=e0*e0;
			tmp[i]=x2[i];
		}
		Merge_sort(tmp, sortwork, 0, total_p-1);
		e0=tmp[Nsmall-1];
		for (i=0; i<total_p; i++)
		{
			if (x2[i]<=e0)
				softspot[i]=1;
		}
	}
	free(tmp);
	free(sortwork);
	free(x2);
	free(G);

	//===================================================
	// plot out softspot
	char *softfile=getfilename(file, ".ss");
	FILE *f=fopen(softfile, "w");
	FILE_NULL(f, softfile);
	for (i=0; i<total_p; i++)
		fprintf(f, "%f %f %u\n", mp[2*i], mp[2*i+1], softspot[i]);
	fclose(f);
	free(mp);
	free(softspot);

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
				   //"#set terminal postscript eps color enhanced\n"
				   //"#set output '%s'\n"
				   "unset key\n"
				   "unset tics\n"
				   "unset colorbox\n"
				   "set size ratio -1\n"
				   //"set pm3d map\n"
				   "set palette defined (0 \"blue\", 1 \"red\")\n"
				   "set title \"%s\\n"
				   "soft spots: Nmode=%d Nsmall=%d\"\n"
				   //"set xlabel 'x'\n"
				   //"set ylabel 'y'\n"
				   , NULLfilename, Nmode, Nsmall);


	//fprintf(pipe, "splot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n", 
	//		softfile, 30.0/sqrt((float)total_p));
	fprintf(pipe, "plot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n", 
			softfile, 30.0/sqrt((float)total_p));

	fclose(pipe);

	if ( remove(softfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", softfile);


	free(softfile);
	free(NULLfilename);
}
