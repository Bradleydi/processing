#include "elasticity.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void el_modulus(const char*file, int b, double a0, int dim)
{
	if (dim<2 || dim>3)
		pERROR("Wrong dimension.");
	// get system volumn
	colloid_base ptid;
	char *gdffile = getfilename (file, ".gdf");
	readgdf(ptid, gdffile);
	free(gdffile);
	
	int i;

	int *size = ptid.get_size();
	float *ptid_ptr = ptid.get_array_pointer();

	float *pointer = ptid_ptr;
	float xmin=*pointer, xmax=xmin;
	float ymin=*(pointer+1), ymax=ymin;

	float tmp;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		tmp = *pointer;
		if ( xmin > tmp )
			xmin = tmp;
		else if ( xmax < tmp )
			xmax = tmp;

		tmp = *(pointer+1);
		if ( ymin > tmp )
			ymin = tmp;
		else if ( ymax < tmp )
			ymax = tmp;
	}
	const float V=(xmax-xmin)*(ymax-ymin);
	ptid.free_memory();
	delete [] size;

	// get coarse-grained strains
	char *tmpfile=getfilename(file, "_cg.gdf");
	char *str = Malloc(char, 30); POINTER_NULL(str);
	sprintf(str, "elasticity_results/%d_el_", b);
	char *cgfile = getfilename(str, tmpfile);
	free(tmpfile);
	free(str);
	
	colloid_base cg;
	readgdf(cg, cgfile);
	free(cgfile);

	int *cgsize = cg.get_size();
	float *cg_ptr=cg.get_array_pointer();

	printf("# check NaN...\n");
	int cc=1;
	for (i=0; i<cgsize[2]; i++)
		if (cg_ptr[i]!=cg_ptr[i])
			printf("# NaN encountered for %d time(s)\n", cc++);
	if (cc!=1)
		exit (1);
	else
		printf("# No NaN found.\n");

	double * dV = Malloc(double, cgsize[1]); POINTER_NULL(dV);
	// twice of orientational fluctuation
	double * _2theta = Malloc(double, cgsize[1]); POINTER_NULL(_2theta);
	pointer=cg_ptr+2;
	for (i=0; i<cgsize[1]; i++)
	{
		// e_xx+e_yy=u_xx+u_yy
		dV[i]=(double)(*pointer+*(pointer+3));
		// 2*theta=(u_yx-u_xy)
		_2theta[i]=(double)((*(pointer+2)-*(pointer+1)));
		pointer+=cgsize[0];
	}
	cg.free_memory();
	
	double meandV = pairwise (cgsize[1], dV)/cgsize[1];
	double mean2dt = pairwise (cgsize[1], _2theta)/cgsize[1]; // mean 2*theta


	double * dV2 = Malloc(double, cgsize[1]); POINTER_NULL(dV2);
	// twice of orientational fluctuation
	double * _2theta2 = Malloc(double, cgsize[1]); POINTER_NULL(_2theta2);
	
	double * dV4 = Malloc(double, cgsize[1]); POINTER_NULL(dV4);
	// twice of orientational fluctuation
	double * _2theta4 = Malloc(double, cgsize[1]); POINTER_NULL(_2theta4);

	for (i=0; i<cgsize[1]; i++)
	{
		dV2[i]=(dV[i]-meandV)*(dV[i]-meandV);
		dV4[i]=dV2[i]*dV2[i];

		_2theta2[i]=(_2theta[i]-mean2dt)*(_2theta[i]-mean2dt);
		_2theta4[i]=_2theta2[i]*_2theta2[i];
	}

	double sddV = pairwise (cgsize[1], dV2)/cgsize[1];
	double sd2dt = pairwise (cgsize[1], _2theta2)/cgsize[1]; // mean (2*theta)^2
	double sdV = pairwise (cgsize[1], dV4)/cgsize[1];
	double s2dt = pairwise (cgsize[1], _2theta4)/cgsize[1]; // mean (2*theta)^4

	free(dV);
	free(dV2);
	free(dV4);
	free(_2theta);
	free(_2theta2);
	free(_2theta4);
	delete [] cgsize;

	// output data
	char *ofile = getfilename (file, ".el");
	FILE * out = fopen(ofile, "a");
	FILE_NULL(out, ofile);
	free(ofile);

	double mu, K, E, nu;
	if (dim==2)
	{
		mu=b*b*a0*a0/V/sd2dt;         // shear modulus
		K =2*b*b*a0*a0/V/sddV - mu;     // bulk modulus
		E=4*K*mu/(K+mu);            // Young's modulus
		nu=(K-mu)/(K+mu);     // Poisson ratio
	}
	else if (dim==3)
	{
		mu=b*b*a0*a0/V/sd2dt;         // shear modulus
		K =2*b*b*a0*a0/V/sddV - 4*mu/3;     // bulk modulus
		E=9*K*mu/(3*K+mu);            // Young's modulus
		nu=(3*K-2*mu)/2/(3*K+mu);     // Poisson ratio
	}
	else
		pERROR("Wrong dimension.");

	fprintf(out, "# b K mu E nu (dV)^2 (2dtheta)^2 dVnonGaussian dthetanonG\n");
	fprintf(out, "%d %.6f %.6f %.6f %.6f %.6f %g %g %.6f %.6f\n",
				 b, 1.0/(double)b, K, mu, E, nu, sddV, sd2dt,
				 sdV/3/sddV/sddV-1, s2dt/3/sd2dt/sd2dt-1);
	fclose(out);
}
