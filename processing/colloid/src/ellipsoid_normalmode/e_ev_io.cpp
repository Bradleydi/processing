#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>


void e_writeev(int Gn, double *E, double *G, const char *file)
{
	char *filename=getfilename(file, ".ev");
	FILE *evf=fopen(filename, "wb");
	FILE_NULL(evf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	header[0]='e'*256+'v'+'e';
	header[1]=Gn;
	
	FwriteN(header, 2, evf, filename);
	free(header);
	FwriteN(E, Gn, evf, filename);

	Gn *= Gn;
	FwriteN(G, Gn, evf, filename);

	Fclose(evf, filename);
	free(filename);
}


void e_readev(int *Gnp, double **Ep, double **Gp, const char *file)
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
	*Gnp=header[1];

	double *E=(double *)malloc((*Gnp)*sizeof(double)); POINTER_NULL(E);
	*Ep=E;
	FreadN(E, *Gnp, evf, filename);

	header[1] *= *Gnp; // Gn*Gn
	double *G=(double *)malloc(header[1]*sizeof(double)); POINTER_NULL(G);
	*Gp=G;
	FreadN(G, header[1], evf, filename);

	Fclose(evf, filename);
	printf("# '%s' loaded\n", filename);
	free(filename);
	free(header);
}


// read N'th mode only
// to avoid read whole .ev file to save memory and add speed
void e_readmode(int N, int *Gnp, double *Ep, double **Gp, const char *file)
{
	if ( N <= 0 ) pERROR("mode index should be greater than 0\n");

	char *filename=getfilename(file, ".ev");
	FILE *evf=fopen(filename, "rb");
	FILE_NULL(evf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	FreadN(header, 2, evf, filename);
	if ( header[0] != ('e'*256+'v'+'e') ) {
		fprintf(stderr, "# Error: '%s' is not a .ev file!\n", filename);
		exit (1);
	}
	*Gnp=header[1];
	free(header);

	int Gn = *Gnp;
	if ( N >= Gn ) {
		fprintf(stderr, "# Error: mode index should be smaller than %d\n", Gn);
		exit (1);
	}

	//========================================================
	// read eigenvalue
	long int offset=(long int)(N*sizeof(double));
	if ( fseek( evf, offset, SEEK_CUR ) != 0 )
		pERROR("reading eigenvalue failed.");

	FreadN(Ep, 1, evf, filename);
	//========================================================
	// read eigen vector
	double *G=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(G);
	*Gp=G;
	offset=(long int)(2*sizeof(int))+(long int)(Gn*(N+1)*sizeof(double));
	if ( fseek( evf, offset, SEEK_SET ) != 0 )
		pERROR("reading eigenvector failed.");
	FreadN(G, Gn, evf, filename);

	Fclose(evf, filename);
	free(filename);
}


// read eigenvalues only
void e_readE(int *Gnp, double **Ep,const char *file)
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
	*Gnp=header[1];
	free(header);

	double *E=(double *)malloc((*Gnp)*sizeof(double)); POINTER_NULL(E);
	*Ep=E;
	FreadN(E, *Gnp, evf, filename);

	Fclose(evf, filename);
	printf("# '%s' loaded\n", filename);
	free(filename);
}
