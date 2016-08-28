#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

void DOS_log(const char *file)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);

	// get omega
	// eigenvalues are in increase order, hence
	// the last two are flopy modes, which are translation
	// freq = 0 // so total mode number is  Gn-2
	int i, flopy=2, total_m=Gn-flopy; // total_mode
	double *omega = Malloc(double, total_m); POINTER_NULL(omega);
	for (i=0; i<total_m; i++) omega[i]=1.0/sqrt(E[i]);
	free(E);

	// log scale omega
	double *logOmega = Malloc(double, total_m); POINTER_NULL(logOmega);
	for (i=0; i<total_m; i++) logOmega[i] = log(omega[i]);
	free(omega);

	const int nbin = (int)sqrt(total_m);
	// since E is increasing as i increasing, so omega decreases
	double binsize = (1.0+1.0e-6)*(logOmega[0]-logOmega[total_m-1])/nbin;
	unsigned int *hist = Calloc(unsigned int, nbin); POINTER_NULL(hist);
	for (i=0; i<total_m; i++) ++hist[(int)((logOmega[i]-logOmega[total_m-1])/binsize)];

	char *filename=getfilename(file, ".dosl");
	FILE *dosf=fopen(filename, "w");
	FILE_NULL(dosf, filename);
	fprintf(dosf, "# density of states:\n"
				  "# omega\tDOS\tDOS/omega\tDOS/omega^2\n");
	double DOS, omg, omg1, omg2;
	for (i=0; i<nbin; i++) {
		omg1=exp(logOmega[total_m-1]+i*binsize);
		omg2=exp(logOmega[total_m-1]+(i+1)*binsize);
		omg=sqrt(omg1*omg2);
		DOS=((double)hist[i])/(total_m)/(omg2-omg1);
		fprintf(dosf, "%g\t%g\t%g\t%g\n", omg, DOS, DOS/omg, DOS/omg/omg);
	}
	free(logOmega);
	Fclose(dosf, filename);
	free(filename);
}
