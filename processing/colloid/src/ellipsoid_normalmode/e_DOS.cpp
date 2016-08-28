#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

/*======================================================
 * cummulative density of state
 *====================================================*/
void e_CDOS(const char *file)
{
	int Gn;
	double *E;
	e_readE(&Gn, &E, file);

	// get omega
	// eigenvalues are in increase order, hence
	// the last two are flopy modes, which are translation
	// freq = 0 // so total mode number is  Gn-2
	int i, flopy=2, total_m=Gn-flopy; // total_mode
	double *omega = Malloc(double, total_m); POINTER_NULL(omega);
	for (i=0; i<total_m; i++) omega[i]=1.0/sqrt(E[i]);
	free(E);

	char *filename=getfilename(file, ".cdos");
	FILE *dosf=fopen(filename, "w");
	FILE_NULL(dosf, filename);
	// since E is increasing as i increasing, so omega decreases
	fprintf(dosf, "# Cummulative density of states:\n"
				  "# omega\tcum_DOS\tcum_DOS/omega\tcum_DOS/omega^2"
				  "\tcum_DOS/omega^3\n");
	double cum_DOS, omg;
	double step=1.0/(double)total_m;
	for (i=0; i<total_m; i++) {
		omg=omega[total_m-i-1];
		cum_DOS=(i+1)*step;
		//fprintf(dosf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n", omg,
		fprintf(dosf, "%g\t%g\t%g\t%g\t%g\n", omg,
				cum_DOS, cum_DOS/omg, cum_DOS/omg/omg, 
				cum_DOS/omg/omg/omg);
	}
	Fclose(dosf, filename);
	free(omega);
	free(filename);
}


void e_DOS(const char *file)
{
	int Gn;
	double *E;
	e_readE(&Gn, &E, file);

	// get omega
	// eigenvalues are in increase order, hence
	// the last two are flopy modes, which are translation
	// freq = 0 // so total mode number is  Gn-2
	int i, flopy=2, total_m=Gn-flopy; // total_mode
	double *omega = Malloc(double, total_m); POINTER_NULL(omega);
	for (i=0; i<total_m; i++) omega[i]=1.0/sqrt(E[i]);
	free(E);

	// average separation
	double avsp = (omega[0]-omega[total_m-1])/(total_m-1);
	
	// find outliers
	// if two freq separate each other a lot, indentify to be a potential
	// outlier. If they forms a cluster, then view as outliers.
	// If one large (>5) cluster exists in high freq part, then all
	// freq. larger than this is viewd as outliers if only a few modes exists
	char *outlier=Calloc(char, total_m-1); POINTER_NULL(outlier);
	for (i=0; i<total_m-1; i++)
		if (omega[i]-omega[i+1] > 2.0*avsp) outlier[i]=1;
	int j;
	for (i=0; i<total_m-6; i++)
	{
		j=0;
		while (outlier[i+j] == 0 && j<5) ++j;
		if (j==5) break;
	}
	const int outlier_bound = i;
	free(outlier);

	const int nbin = (int)sqrt(total_m-outlier_bound);
	// since E is increasing as i increasing, so omega decreases
	double binsize = (1.0+1.0e-6)*omega[outlier_bound]/nbin;
	unsigned int *hist = Calloc(unsigned int, nbin); POINTER_NULL(hist);
	//printf("%d %d %d\n", outlier_bound, nbin, (int)(omega[i]/binsize));
	for (i=outlier_bound; i<total_m; i++) ++hist[(int)(omega[i]/binsize)];
	//printf("OK\n");
	free(omega);

	char *filename=getfilename(file, ".dos");
	FILE *dosf=fopen(filename, "w");
	FILE_NULL(dosf, filename);
	fprintf(dosf, "# density of states:\n"
				  "# outlier_bound = %d\n"
				  "# omega\tDOS\tDOS/omega\tDOS/omega^2\n", outlier_bound);
	double DOS, omg;
	for (i=0; i<nbin; i++) {
		omg=(i+0.5)*binsize;
		DOS=((double)hist[i])/(total_m-outlier_bound)/binsize;
		//fprintf(dosf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n", omg,
		fprintf(dosf, "%g\t%g\t%g\t%g\n", omg, DOS, DOS/omg, DOS/omg/omg);
	}
	Fclose(dosf, filename);
	free(filename);
}
