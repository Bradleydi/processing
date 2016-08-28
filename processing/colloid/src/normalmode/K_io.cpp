#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;


void writeK(int &Kn, double *K, const char *file)
{
	char *filename=getfilename(file, ".K");
	
	FILE *dcmf=fopen(filename, "wb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int));
	POINTER_NULL(header);
	header[0]='K'*256;
	header[1]=Kn;
	
	if ( fwrite(header, sizeof(header[0]), 2, dcmf) != 2 )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	free(header);


	int count=Kn*Kn;
	if ( (int)fwrite(K, sizeof(K[0]), count, dcmf) != count )
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


void readK(int &Kn, double **Kp, const char *file)
{
	char *filename=getfilename(file, ".K");
	FILE *dcmf=fopen(filename, "rb");
	FILE_NULL(dcmf, filename);

	int * header=(int *)malloc(2*sizeof(int));
	POINTER_NULL(header);
	
	if ( fread(header, sizeof(header[0]), 2, dcmf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		exit (1);
	}

	if ( header[0] != ('K'*256) )
	{
		fprintf(stderr, "# Error: %s is not a .K file!\n", filename);
		exit (1);
	}

	Kn=header[1];
	free(header);


	int count=Kn*Kn;
	double *K=(double *)malloc(count*sizeof(double));
	POINTER_NULL(K);

	*Kp=K;
	if ( (int)fread(K, sizeof(double), count, dcmf) != count )
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
