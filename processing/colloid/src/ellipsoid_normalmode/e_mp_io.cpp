#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>


void e_writemp(const int &N, float *mp,const char *file)
{
	char *filename=getfilename(file, ".mp");
	FILE *mpf=fopen(filename, "wb");
	FILE_NULL(mpf, filename);

	int * header=(int *)malloc(2*sizeof(int));
	header[0]='m'*256+'p'+'e'; // to distinguish from spherical data
	header[1]=N;
	
	FwriteN(header, 2, mpf, filename);
	FwriteN(mp, 3*N, mpf, filename);
	Fclose(mpf, filename);

	free(header);
	free(filename);
}


void e_readmp(int *Np, float **mp_p, const char *file)
{
	char *filename=getfilename(file, ".mp");
	FILE *mpf=fopen(filename, "rb");
	FILE_NULL(mpf, filename);

	int *header=Malloc(int, 2); POINTER_NULL(header);
	FreadN(header, 2, mpf, filename);
	if ( header[0] != ('m'*256+'p'+'e') )
	{
		fprintf(stderr, "# Error: %s is not a .mp file!\n", filename);
		exit (1);
	}
	*Np=header[1];
	header[1] *= 3;
	float *mp=Malloc(float, header[1]); POINTER_NULL(mp);
	*mp_p=mp;
	FreadN(mp, header[1], mpf, filename);
	Fclose(mpf, filename);

	free(header);
	free(filename);
}
