#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;


void writedcm(int &Gn, double *G, const char *file)
{
	char *filename=getfilename(file, ".dcm");
	
	FILE *dcmf=fopen(filename, "wb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int));
	POINTER_NULL(header);
	header[0]='d'*256+'c'+'m';
	header[1]=Gn;
	
	if ( fwrite(header, sizeof(header[0]), 2, dcmf) != 2 )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	free(header);


	int count=Gn*Gn;
	if ( (int)fwrite(G, sizeof(G[0]), count, dcmf) != count )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	if (fclose(dcmf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		exit (1);
	}

	free(filename);
}


void readdcm(int &Gn, double **Gp, const char *file)
{
	char *filename=getfilename(file, ".dcm");
	FILE *dcmf=fopen(filename, "rb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int));
	POINTER_NULL(header);
	
	if ( fread(header, sizeof(header[0]), 2, dcmf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		exit (1);
	}

	if ( header[0] != ('d'*256+'c'+'m') )
	{
		fprintf(stderr, "# Error: %s is not a .dcm file!\n", filename);
		exit (1);
	}

	Gn=header[1];
	free(header);


	int count=Gn*Gn;
	double *G=(double *)malloc(count*sizeof(double));
	POINTER_NULL(G);

	*Gp=G;
	if ( (int)fread(G, sizeof(double), count, dcmf) != count )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		exit (1);
	}

	if (fclose(dcmf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		exit (1);
	}

	free(filename);
}
