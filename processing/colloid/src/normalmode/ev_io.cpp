#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <ctime>
#include <cstdlib>

using namespace std;

/* #define WHERESTR  "# Error: [%s, line %d]: "
 * #define WHEREARG  __FILE__, __LINE__
 * #define DEBUGPRINT2(...)       fprintf(stderr, __VA_ARGS__)
 * #define pERROR(_message)  \
 * { DEBUGPRINT2(WHERESTR _message"\n", WHEREARG); exit(1); }
*/

void writeev(int &Gn, double *E, double *G, const char *file)
{
	const char subfix []=".ev";
	char *filename=getfilename(file, subfix);
	
	FILE *evf=fopen(filename, "wb");
	if (evf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	int * header=(int *)malloc(2*sizeof(int));
	header[0]='e'*256+'v';
	header[1]=Gn;
	
	if ( fwrite(header, sizeof(header[0]), 2, evf) != 2 )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	free(header);

	if ( (int)fwrite(E, sizeof(E[0]), Gn, evf) != Gn )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	int count=Gn*Gn;
	if ( (int)fwrite(G, sizeof(G[0]), count, evf) != count )
	{
		fprintf(stderr, "# Error: write to file %s failed!\n", filename);
		exit (1);
	}

	if (fclose(evf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		exit (1);
	}

	free(filename);
}


int readev(int &Gn, double **Ep, double **Gp, const char *file)
{
	const char subfix []=".ev";
	char *filename=getfilename(file, subfix);

	FILE *evf=fopen(filename, "rb");
	if (evf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		return 1;
	}

	int * header=(int *)malloc(2*sizeof(int));
	
	if ( fread(header, sizeof(header[0]), 2, evf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if ( header[0] != ('e'*256+'v') )
	{
		fprintf(stderr, "# Error: %s is not a .ev file!\n", filename);
		return 1;
	}

	Gn=header[1];
	free(header);

	double *E=(double *)malloc(Gn*sizeof(double));
	if (E==NULL)
	{
		fprintf(stderr, "# Error: double *E initialized wrong!");
		return 1;
	}
	*Ep=E;

	if ( (int)fread(E, sizeof(double), Gn, evf) != Gn )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	int count=Gn*Gn;
	double *G=(double *)malloc(count*sizeof(double));
	if (G==NULL)
	{
		fprintf(stderr, "# Error: double *G initialized wrong!");
		return 1;
	}
	*Gp=G;
	if ( (int)fread(G, sizeof(double), count, evf) != count )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if (fclose(evf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		return 1;
	}

	free(filename);
	return 0;
}


// read N'th mode only
// to avoid read whole .ev file to save memory and add speed
int readmode(int &N, int &Gn, double &E, double **Gp, const char *file)
{
	if ( N <= 0 )
	{
		fprintf(stderr, "# Error: mode index should be greater than 0\n");
		return 1;
	}
	const char subfix []=".ev";
	char *filename=getfilename(file, subfix);

	FILE *evf=fopen(filename, "rb");
	if (evf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		return 1;
	}

	int * header=(int *)malloc(2*sizeof(int));
	
	if ( fread(header, sizeof(header[0]), 2, evf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if ( header[0] != ('e'*256+'v') )
	{
		fprintf(stderr, "# Error: %s is not a .ev file!\n", filename);
		return 1;
	}

	Gn=header[1];
	free(header);

	if ( N >= Gn )
	{
		fprintf(stderr, "# Error: mode index should be smaller than %d\n", Gn);
		return 1;
	}

	//========================================================
	// read eigenvalue
	long int offset=(long int)(N*sizeof(double));
	if ( fseek( evf, offset, SEEK_CUR ) != 0 )
		pERROR("reading eigenvalue failed.");
	double * Ep=&E;

	if ( (int)fread(Ep, sizeof(double), 1, evf) != 1 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}
	//========================================================
	// read eigen vector

	double *G=(double *)malloc(Gn*sizeof(double));
	if (G==NULL)
	{
		fprintf(stderr, "# Error: double *G initialized wrong!");
		return 1;
	}
	*Gp=G;
	// bug here 
	/*
	offset=(long int)((2+Gn*(N+1))*sizeof(double));
	*/
	offset=(long int)(2*sizeof(int))+(long int)(Gn*(N+1)*sizeof(double));
	if ( fseek( evf, offset, SEEK_SET ) != 0 )
		pERROR("reading eigenvector failed.");

	if ( (int)fread(G, sizeof(double), Gn, evf) != Gn )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		return 1;
	}

	if (fclose(evf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		return 1;
	}

	free(filename);
	return 0;
}


// read eigenvalues only
void readE(int &Gn, double **Ep,const char *file)
{
	const char subfix []=".ev";
	char *filename=getfilename(file, subfix);

	FILE *evf=fopen(filename, "rb");
	if (evf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	int * header=(int *)malloc(2*sizeof(int));
	
	if ( fread(header, sizeof(header[0]), 2, evf) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		exit (1);
	}

	if ( header[0] != ('e'*256+'v') )
	{
		fprintf(stderr, "# Error: %s is not a .ev file!\n", filename);
		exit (1);
	}

	Gn=header[1];
	free(header);

	double *E=(double *)malloc(Gn*sizeof(double));
	if (E==NULL)
	{
		fprintf(stderr, "# Error: double *E initialized wrong!");
		exit (1);
	}
	*Ep=E;

	if ( (int)fread(E, sizeof(double), Gn, evf) != Gn )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", filename);
		exit (1);
	}

	if (fclose(evf)!=0)
	{
		fprintf(stderr, "# Error: closing file %s failed!\n", filename);
		exit (1);
	}

	printf("# File \"%s\" loaded!\n", filename);
	free(filename);
}
