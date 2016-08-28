#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

/* participation ratio */

void e_PR(const char *file)
{	
	int Gn;
	double *E, *G;
	e_readev(&Gn, &E, &G, file);

	char *filename=getfilename(file, ".pr");
	FILE *prf=fopen(filename, "w");
	FILE_NULL(prf, filename);
	
	fprintf(prf, "# participation ratio (PR)\n");
	const int total_p=Gn/3;
	// translational and orientational participation ratio
	// definition e4(trl)/e2(trl), e4(ang)/e2(ang)
	// e2(trl)+e2(ang)=1
	double *e4trl=(double *)malloc(total_p*sizeof(double)); POINTER_NULL(e4trl);
	double *e2ang=(double *)malloc(total_p*sizeof(double)); POINTER_NULL(e2ang);
	double *e4ang=(double *)malloc(total_p*sizeof(double)); POINTER_NULL(e4ang);
	double eij;
	fprintf(prf, "# lambda\tomega\te2ang\tPR(translational)\tPR(orien)\n");
	int i, j, total_m=Gn-2, total_p2=total_p*2;
	double *Gptr=G;
	for (i=0; i<total_m; i++) {
		for (j=0; j<total_p; j++) {
			eij=*Gptr; // e_xi
			e4trl[j]=eij*eij; // (e_xi)^2
			eij=*(Gptr+total_p); // e_yi
			e4trl[j] += eij*eij; // (e_xi)^2+(e_yi)^2=(e_i)^2
			e4trl[j] *= e4trl[j];

			eij=*((Gptr++)+total_p2); // e_ang
			e2ang[j] = eij*eij;
			e4ang[j] = e2ang[j]*e2ang[j];
		}
		Gptr+=total_p2; // be careful here, since Gptr-->x
		eij=pairwise(total_p, e2ang);
		fprintf(prf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\n", E[i],
				1.0/sqrt(E[i]), eij,
				//(1.0-eij)/pairwise(total_p, e4trl)/total_p,
				//eij/pairwise(total_p, e4ang)/total_p);
				(1.0-eij)*(1.0-eij)/pairwise(total_p, e4trl)/total_p,
				eij*eij/pairwise(total_p, e4ang)/total_p);
		// Here is a bug
	}

	free(e4trl); free(e2ang); free(e4ang);
	free(E); free(G);
	Fclose(prf, filename);
	free(filename);
}
