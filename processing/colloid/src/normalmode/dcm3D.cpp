#include "normalmode.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void dcm3D(const char *file)
{
	int Gn;
	double *G;
	readdcm(Gn, &G, file);

	// construct 3D displacement covariance matrix
	int Dn=3*Gn/2;
	double *D = Malloc(double, Dn*Dn);
	POINTER_NULL(D);

	int i, j, total_p=Gn/2, Gij=0/*x_ix_j position*/, Dij=0;
	for (i=0; i<total_p; i++)
	{
		for (j=0; j<total_p; j++)
		{
			D[Dij]=G[Gij]; 			// D_xx = G_xx
			D[Dij+1]=G[Gij+1]; 		// D_xy = G_xy
			D[Dij+2]=G[Gij+1]; 		// D_xz = G_xy

			Dij+=Dn;
			D[Dij]=G[Gij+1]; 		// D_yx = G_xy
			D[Dij+1]=G[Gij+Gn+1]; 	// D_yy = G_yy
			D[Dij+2]=G[Gij+1]; 		// D_yz = G_xy

			Dij+=Dn;
			D[Dij]=G[Gij+1]; 		// D_zx = G_xy
			D[Dij+1]=G[Gij+1]; 		// D_zy = G_xy
			D[Dij+2]=(G[Gij]+G[Gij+Gn+1])/2.0; // D_zz = ( G_xx + G_yy ) / 2
			
			Dij+=Dn;
			Gij+=2*Gn;
		}
		Dij -= Dn*Dn-3;
		Gij -= Gn*Gn-2;
	}
	free(G);
	writedcm3D(Dn, D, file);
	free(D);
}
