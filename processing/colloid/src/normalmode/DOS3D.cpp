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
void DOS3D(const char *file)
{
	int Dn;
	double *D;
	readdcm3D(Dn, &D, file);

	double *E=(double *)malloc(Dn*sizeof(double));
	POINTER_NULL(E);
	dsyev(D, Dn, E);
	if (E[0]<0)
	{
		printf("# First eigenvalue: %6.10f\n", E[0]);
		//pERROR("negative eigenvalue(s)");
	}

	char *filename=getfilename(file, ".dos3");

	FILE *dosf=fopen(filename, "w");
	if (dosf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	free(filename);

	// get omega
	// since the first 3 are flopy modes, which are translation
	// so total mode number is  Dn-3
	double *omega = (double *) malloc ((Dn-3)*sizeof(double));
	POINTER_NULL(omega);
	
	int i;
	for (i=0; i<Dn-3; i++)
		omega[i]=1/sqrt(E[i+3]);

	// since E is increasing as i increasing, so
	// omega decreases

	fprintf(dosf, "# Cummulative density of states:\n"
				  "# omega\tcum_DOS\tcum_DOS/omega\tcum_DOS/omega^2"
				  "\tcum_DOS/omega^3\n");
	double cum_DOS, omg;
	double step=1.0/(double)(Dn-3);
	for (i=0; i<Dn-3; i++)
	{
		//DOS=(double)(hist[i])/binsize[i]/(Gn-2);
		//cum+=(double)(hist[i])/(Gn-2);
		omg=omega[Dn-4-i];
		cum_DOS=(i+1)*step;
		//fprintf(dosf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n", omg,
		fprintf(dosf, "%g\t%g\t%g\t%g\t%g\n", omg,
				cum_DOS, cum_DOS/omg, cum_DOS/omg/omg, 
				cum_DOS/omg/omg/omg);
	}
	
	fclose(dosf);
	free(omega);
}
