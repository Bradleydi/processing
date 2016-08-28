#include "normalmode.h"

#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"
#include "statistics.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void participation_ratio_w(int Gn, double *E,
		double *G, const char *file)
{	
	const char subfix[]=".pr";
	char *filename=getfilename(file, subfix);

	FILE *prf=fopen(filename, "w");
	if (prf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}
	
	fprintf(prf, "# participation ratio\n");
	const int total_p=Gn/2;
	double *e4=(double *)malloc(total_p*sizeof(double));
	if (e4==NULL)
		pERROR("double *e4 initialized error!");
	double eij;
	fprintf(prf, "# lambda\tomega\tparticipation_ratio\n");
	int i, j;
	/* ======== Aug 9, 2011 ==============
	for (i=2; i<Gn; i++)
	{
		for (j=0; j<total_p; j++)
		{
			eij=G[i*Gn+2*j];
			eij*=eij;
			e4[j]=eij;
			eij=G[i*Gn+(2*j+1)];
			eij*=eij;
			e4[j]+=eij;
			e4[j]=e4[j]*e4[j];
		}
		fprintf(prf, "%6.6f\t", E[i]);
		E[i]=1/sqrt(E[i]);
		fprintf(prf, "%6.6f\t%6.6f\n", 
				E[i], 1/pairwise(total_p, e4)/total_p);
	}
	*/
	double *Gptr=G+2*Gn;
	for (i=2; i<Gn; i++)
	{
		for (j=0; j<total_p; j++)
		{
			eij=*(Gptr++); // e_xi
			e4[j]=eij*eij; // (e_xi)^2
			eij=*(Gptr++); // e_yi
			e4[j]+=eij*eij; // (e_xi)^2+(e_yi)^2=(e_i)^2
			e4[j]=e4[j]*e4[j];
		}
		fprintf(prf, "%6.6f\t%6.6f\t%6.6f\n", E[i],
				1/sqrt(E[i]), 1/pairwise(total_p, e4)/total_p);
	}

	fclose(prf);
	free(e4);
	free(filename);
}


void participation_ratio(colloid_base &ptid, const char *file,
		bool Remove_Drift)
{
	// Participation ratio
	int Gn;
	double *E, *G;
	if ( readev(Gn, &E, &G, file) != 0 )
	{	
		printf("# Calling normal_mode!\n");
		normal_mode(ptid, file, Remove_Drift);
	}
	else
	{
		participation_ratio_w(Gn, E, G, file);
		free(E);
		free(G);
	}
}


void participation_ratio_smooth(const char *file, double h)
{
	// read ev
	int Gn;
	double *E, *G;
	if ( readev(Gn, &E, &G, file) != 0 )
		pERROR("reading error.");

	char *filename=getfilename(file, ".prs");
	FILE *prf=fopen(filename, "w");
	FILE_NULL(prf, filename);
	free(filename);
	const int total_p=Gn/2;
	double *e4=(double *)malloc(total_p*sizeof(double));
	POINTER_NULL(e4);
	double eij;
	fprintf(prf, "# omega\tsmoothed_pr\tparticipation_ratio\tlambda\n");
	int i, j;
	double *Gptr=G+2*Gn;
	double *pr=Malloc(double, Gn-2); POINTER_NULL(pr);
	for (i=2; i<Gn; i++)
	{
		for (j=0; j<total_p; j++)
		{
			eij=*(Gptr++); // e_xi
			e4[j]=eij*eij; // (e_xi)^2
			eij=*(Gptr++); // e_yi
			e4[j]+=eij*eij; // (e_xi)^2+(e_yi)^2=(e_i)^2
			e4[j]=e4[j]*e4[j];
		}
		pr[i-2]=1/pairwise(total_p, e4)/total_p;
	}
	free(e4);
	free(G);

	double *prs=Malloc(double, Gn-2); POINTER_NULL(prs);
	locallinear(Gn-2, E+2, pr, prs, h, GAUSSIAN);

	for (i=Gn-1; i>=2; i--)
	{
		fprintf(prf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\n", 
				1/sqrt(E[i]), prs[i-2], pr[i-2], E[i]);
	}

	fclose(prf);
	free(E);
	free(pr);
	free(prs);
}
