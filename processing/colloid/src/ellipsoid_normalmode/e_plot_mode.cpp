#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

void e_plot_mode(const char *file, int n)
{
	if (n<0) pERROR("wrong mode id");
	// get mean position from .mp file
	int total_p;
	float *mp;
	char *mpfilename=getfilename(file, ".mp");
	e_readmp(&total_p, &mp, mpfilename);
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);

	// check mode n whether valid
	if (n>=3*total_p) {
		fprintf(stderr, "# Error: max mode id = %d\n", 3*total_p-1);
		exit (1);
	}
	
	// get n'th normal mode from .ev file
	char *evfilename=getfilename(file, ".ev");
	int Gn;
	double E;
	double *G;
	//n=total_p*3-1-n;
	e_readmode(total_p*3-1-n, &Gn, &E, &G, evfilename);
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	E=1.0/sqrt(E);

	if ( Gn != 3*total_p ) pERROR("reading files wrong.");

	//===================================================
	// print out the data
	char modesubfix [20];
	sprintf(modesubfix, "_%05d.mode", n);
	char *modefilename = getfilename(file, modesubfix);
	FILE *modefile=fopen(modefilename, "w");
	FILE_NULL(modefile, modefilename);

	fprintf(modefile, "# %d'th normal mode\n"
			"# eigenvalue = %f\n", n, E);
	int i, total_p2=total_p*2;
	float *pmp=mp;
	double *pG=G+total_p2, eij;
	double *e2ang = Malloc(double, total_p); POINTER_NULL(e2ang);
	double max_eag=*pG, min_eag=*pG;
	for (i=0; i<total_p; i++) {
		eij=*(pG++); // e_ang
		e2ang[i] = eij*eij;
		if (eij > max_eag ) max_eag=eij;
		else if (eij < min_eag ) min_eag=eij;
	}
	if (-min_eag > max_eag) max_eag=-min_eag;
	eij=pairwise(total_p, e2ang); // e_ang^2 total
	free(e2ang);
	double ratio=sqrt((double)(Gn))/0.5;///(1.0-eij);
	pG=G;
	for (i=0; i<total_p; i++)
	{
		fprintf(modefile, "%f %f %f %f %f\n", *pmp, *(pmp+total_p), 
				(*pG)*ratio, (*(pG+total_p))*ratio,
				*(pG+total_p2));
		++pmp; ++pG;
	}
	Fclose(modefile, modefilename);
	free(mp); free(G);

	//=========================================================
	// plot by GNUplot, using pipe
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}
	char *NULLfilename = getfilename(file, "");
	sprintf(modesubfix, "_mode_%05d.eps", n);
	char *epsfilename = getfilename(file, modesubfix);

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set terminal postscript eps color enhanced\n"
				   "set output \"%s\"\n"
				   "set size ratio -1\n"
				   "set autoscale\n"
				   "set title \"%s\\n"
				   " %d'th normal mode\\nomega=%f\"\n"
				   "unset key\n"
				   "unset xtics\n"
				   "unset ytics\n"
				   "set border 0\n"
				   "set cbrange [%f:%f]\n"
				   "set palette rgbformulae 33,13,10\n"
				   "#set xlabel 'x/a_0'\n"
				   "#set ylabel 'y/a_0'\n", 
				   epsfilename, NULLfilename, n, E, -max_eag, max_eag);
	const float pointsize=30.0/sqrt((float)total_p);
	//fprintf(pipe, "plot '%s' using 1:2:3:4 with vectors, "
	//		      "\"\" using 1:2:5 with points palette pt 6 ps %.2f\n\n", 
	fprintf(pipe, "plot '%s' using 1:2:5 with points palette pt 7 ps %.2f, "
			      "\"\" using 1:2:3:4 with vectors lt -1 lc rgb 'black' lw 1.5\n\n", 
				  modefilename, pointsize);

	fclose(pipe);

	//if ( remove(modefilename) != 0 )
	//	fprintf(stderr, "# Error: file '%s' cannot be removed.\n", modefilename); 
	free(modefilename);
	free(NULLfilename);
	free(epsfilename);
}
