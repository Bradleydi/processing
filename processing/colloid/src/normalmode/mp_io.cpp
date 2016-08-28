#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;


void writemp(const int &N, float *mp,const char *file)
{
	const char subfix []=".mp";
	char *filename=getfilename(file, subfix);
	
	FILE *mpf=fopen(filename, "wb");
	if (mpf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	int * header=(int *)malloc(2*sizeof(int));
	header[0]='m'*256+'p';
	header[1]=N;
	
	if ( fwrite(header, sizeof(header[0]), 2, mpf) != 2 )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	free(header);

	if ( (int)fwrite(mp, sizeof(mp[0]), 2*N, mpf) != 2*N )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	if (fclose(mpf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		exit (1);
	}

	free(filename);
}


int readmp(int &N, float **mp_p, const char *file)
{
	const char subfix []=".mp";
	char *filename=getfilename(file, subfix);

	FILE *mpf=fopen(filename, "rb");
	if (mpf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		return 1;
	}

	int * header=(int *)malloc(2*sizeof(int));
	
	if ( fread(header, sizeof(header[0]), 2, mpf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if ( header[0] != ('m'*256+'p') )
	{
		fprintf(stderr, "# Error: %s is not a .mp file!\n", filename);
		return 1;
	}

	N=header[1];
	free(header);

	float *mp=(float *)malloc(2*N*sizeof(float));
	if (mp==NULL)
	{
		fprintf(stderr, "# Error: double *mp initialized wrong!");
		return 1;
	}
	*mp_p=mp;

	if ( (int)fread(mp, sizeof(float), 2*N, mpf) != 2*N )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if (fclose(mpf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		return 1;
	}

	free(filename);
	return 0;
}
