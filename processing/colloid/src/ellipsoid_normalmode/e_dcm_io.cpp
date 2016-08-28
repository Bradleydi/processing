#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>


void e_writedcm(int Gn, double *G, const char *file)
{
	char *filename=getfilename(file, ".dcm");
	FILE *dcmf=fopen(filename, "wb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	header[0]='d'*256+'c'+'m'+'e';
	header[1]=Gn;
	
	FwriteN(header, 2, dcmf, filename);

	header[1] = Gn*Gn;
	FwriteN(G, header[1], dcmf, filename);

	Fclose(dcmf, filename);
	free(filename);
	free(header);
}


void e_readdcm(int *Gnp, double **Gp, const char *file)
{
	char *filename=getfilename(file, ".dcm");
	FILE *dcmf=fopen(filename, "rb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	
	FreadN(header, 2, dcmf, filename);

	if ( header[0] != ('d'*256+'c'+'m'+'e') )
	{
		fprintf(stderr, "# Error: '%s' is not a .dcm file!\n", filename);
		exit (1);
	}
	*Gnp=header[1];

	header[1] *= *Gnp; // Gn*Gn
	double *G=(double *)malloc(header[1]*sizeof(double)); POINTER_NULL(G);

	*Gp=G;
	FreadN(G, header[1], dcmf, filename);

	Fclose(dcmf, filename);
	free(filename);
	free(header);
}
