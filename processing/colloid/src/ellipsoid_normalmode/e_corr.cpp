#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

/* correlations for:
 * (1) translational and orientational parts
 * (2) translational and other local order parameters, (nematic ...)
 * (3) orientational and ...
 * Here only (1) is implemented
 *
 * Definition of (1)
 * C_TO=< ( T - <T> ) * ( O - <O> ) > / sqrt(corr(T) * corr(O))
 * here T is the length of the trans. part, and O is the orien. part
 */

void e_corr(const char *file)
{	
	int Gn;
	double *E, *G;
	e_readev(&Gn, &E, &G, file);

	char *filename=getfilename(file, ".corr");
	FILE *fp=fopen(filename, "w");
	FILE_NULL(fp, filename);
	
	fprintf(fp, "# mode correlation between trans. and orien. parts\n"
			"# freq. C_TO, stdev(T), stdev(O)\n");
	const int total_p=Gn/3;
	int i, j, total_m=Gn-2;
	double *xG=G, *yG=xG+total_p, *OG=yG+total_p;
	double avL, avO;
	for (i=0; i<total_m; i++)
	{
		// calculated the length of trans. part
		for (j=0; j<total_p; j++)
		{
			xG[j]=sqrt(xG[j]*xG[j]+yG[j]*yG[j]); // save result to G[i]
			if ( OG[j] < 0.0 ) OG[j] = -OG[j];
		}
		// calculate the average
		avL = pairwise(total_p, xG)/total_p;
		avO = pairwise(total_p, OG)/total_p;

		// calculate the correlation
		for (j=0; j<total_p; j++) {
			yG[j] = (xG[i]-avL)*(OG[i]-avO);
			xG[j] = (xG[i]-avL)*(xG[i]-avL);
			OG[j] = (OG[i]-avO)*(OG[i]-avO);
		}
		avL=sqrt(pairwise(total_p, xG)/total_p);
		avO=sqrt(pairwise(total_p, OG)/total_p);

		fprintf(fp, "%g %g %g %g\n", 1.0/sqrt(E[i]),
				pairwise(total_p, yG)/total_p/avL*avO,
				avL, avO);
		xG += Gn;
		yG += Gn;
		OG += Gn;
	}

	free(E); free(G);
	Fclose(fp, filename);
	free(filename);
}
