#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"
#include "statistics.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
	/*======================================================
	 * Eigenvalue manipulation
	 *====================================================*/

// cumulative
void DOS_w(int Gn, double *E, const char *file)
{
	const char subfix[]=".dos";
	char *filename=getfilename(file, subfix);

	FILE *dosf=fopen(filename, "w");
	if (dosf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	free(filename);

	// get omega
	// since the first two are flopy modes, which are translation
	// so total mode number is  Gn-2
	double *omega = (double *) malloc ((Gn-2)*sizeof(double));
	POINTER_NULL(omega);
	
	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);

	// since E is increasing as i increasing, so
	// omega decreases

	fprintf(dosf, "# Cummulative density of states:\n"
				  "# omega\tcum_DOS\tcum_DOS/omega\tcum_DOS/omega^2"
				  "\tcum_DOS/omega^3\n");
	double cum_DOS, omg;
	double step=1.0/(double)(Gn-2);
	for (i=0; i<Gn-2; i++)
	{
		//DOS=(double)(hist[i])/binsize[i]/(Gn-2);
		//cum+=(double)(hist[i])/(Gn-2);
		omg=omega[Gn-3-i];
		cum_DOS=(i+1)*step;
		//fprintf(dosf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n", omg,
		fprintf(dosf, "%g\t%g\t%g\t%g\t%g\n", omg,
				cum_DOS, cum_DOS/omg, cum_DOS/omg/omg, 
				cum_DOS/omg/omg/omg);
	}
	
	fclose(dosf);
	free(omega);
}

void DOS(colloid_base &ptid, const char *file, bool Remove_Drift)
{
	int Gn;
	double *E, *G;
	if ( readev(Gn, &E, &G, file) != 0 )
	{	
		printf("# Calling normal_mode!\n");
		normal_mode(ptid, file, Remove_Drift);
	}
	else
	{
		free(G);

		DOS_w(Gn, E, file);
		free(E);
	}
}


// logarithm binning of DOS
void DOSlog(const char *file)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);

	double *omega = (double *) malloc ((Gn-2)*sizeof(double));
	POINTER_NULL(omega);
	
	int i, j;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);
	free(E);

	const int Nmeanspacing=(int)sqrt(Gn-2);
	const double binbase=Nmeanspacing*(omega[0]-omega[Gn-3])/(Gn-2);

	// get dts
	int bin [150];
	bin[0]=0;   // bin[0] is set to 0, not 1
	int count=1;
	for (i=0; i<150; i++)
	{
		j=(int)(pow(1.15,i)+0.5);
		if ( j*Nmeanspacing>Gn-2 )
			break;
		if (j==bin[count-1])
			continue;
		bin[count++]=j;
	}
	
	char *filename=getfilename(file, ".dosl");
	FILE *dosf=fopen(filename, "w");
	FILE_NULL(dosf, filename);
	free(filename);

	fprintf(dosf, "# omega DOS histogram\n");

	double *sqrtbin=(double *)malloc(count*sizeof(double));
	POINTER_NULL(sqrtbin);
	for (i=0; i<count; i++)
		sqrtbin[i]=sqrt((double)bin[i]);

	double omg;
	unsigned int hist=0;
	j=1;
	for (i=Gn-3; i>=0; i--)
	{
		omg=(omega[i]-omega[Gn-3])/binbase;
		if ( omg > bin[count-1] )
		{
			fprintf(dosf, "%.6g\t%.6g\t%d\n", 
					omega[Gn-3]+sqrtbin[j-1]*sqrtbin[j]*binbase,
					(double)(hist)/(Gn-2)/(bin[j]-bin[j-1])/binbase,
					hist);
			break;
		}

		if ( omg < bin[j] )
			++hist;
		else
		{
			fprintf(dosf, "%g\t%g\t%d\n", 
					omega[Gn-3]+sqrtbin[j-1]*sqrtbin[j]*binbase,
					(double)(hist)/(Gn-2)/(bin[j]-bin[j-1])/binbase,
					hist);
			hist=1;
			while ( omg >= bin[j] )
				++j;
		}
	}
	fclose(dosf);
	free(omega);
	free(sqrtbin);
}


// kernel density estimation
void DOS_kdensity(const char *file, double h, int nbin)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);

	double *omega = (double *) malloc ((Gn-2)*sizeof(double));
	POINTER_NULL(omega);
	
	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);
	free(E);
	
	if (nbin<=0)
		nbin=(int)sqrt(Gn-2);
	double *density = Malloc (double, 2*nbin);
	kdensity(Gn-2, omega, nbin, density, h, GAUSSIAN);
	free(omega);

	char *filename=getfilename(file, ".dosd");
	FILE *dosf=fopen(filename, "w");
	FILE_NULL(dosf, filename);
	free(filename);

	fprintf(dosf, "# omega DOS_get_from_kernel_density_estimation\n");
	for (i=0; i<nbin; i++)
		fprintf(dosf, "%6.6f %6.6f\n", density[2*i], density[2*i+1]);

	fclose(dosf);
	free(density);
}
