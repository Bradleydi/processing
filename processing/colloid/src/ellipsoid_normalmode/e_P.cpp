#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

void e_P(const char* file)
{
	char *filename=getfilename(file, ".ev");
	FILE *evf=fopen(filename, "rb");
	FILE_NULL(evf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	FreadN(header, 2, evf, filename);
	if ( header[0] != ('e'*256+'v'+'e') ) {
		fprintf(stderr, "# Error: '%s' is not a .ev file!\n", filename);
		exit (1);
	}
	int Gn=header[1];
	int total_p=Gn/3;
	free(header);
	
	double *E=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(E);
	double *G=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(G);
	double *xG=G, *yG=xG+total_p, *OG=yG+total_p;
	// read eigenvalue
	FreadN(E, Gn, evf, filename);
	int i, j;
	char *Pfilename=getfilename(file, ".P");
	FILE *Pfp = fopen(Pfilename, "w");
	FILE_NULL(Pfp, Pfilename);
	fprintf(Pfp, "# Omega Px Py PO\n");
	for (i=0; i<Gn; i++)
	{
		// read eigen vector
		FreadN(G, Gn, evf, filename);
		for (j=0; j<Gn; j++) G[j] *= G[j];
		fprintf(Pfp, "%g %g %g %g\n", 1.0/sqrt(E[i]), pairwise(total_p, xG),
				pairwise(total_p, yG), pairwise(total_p, OG));
	}
	Fclose(evf, filename);
	Fclose(Pfp, Pfilename);
	free(filename);
	free(Pfilename);
	free(E); free(G);
}
