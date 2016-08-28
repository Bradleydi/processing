#include "normalmode.h"

#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void U(const char *file, double a0, bool Remove_Drift)
{
	int i, j, k, t, p;
	// read gdf 
	char *gdffilename=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffilename);
	free(gdffilename);

	int *size = ptid.get_size();
	float *ptid_ptr = ptid.get_array_pointer();
	
	const int total_p = (int)ptid_ptr[size[2]-1]+1;

	const int dim = 2;
	if (dim+2>size[0])
		pERROR("Not a 2D data or tracked data!");
	
	const int ti=size[0]-2; //time index in the array
	const int ii=ti+1; // particle index in the array

	float *pointer=ptid_ptr+ti;
	float tt = *pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( tt < *pointer )
			tt = *pointer;
	}
	const int total_t = (int)tt + 1 ;
	
	if ( total_t < 2*total_p )
	{
		fprintf(stderr, 
			"# Error: there should be enough independent configurations!\n"
			"#        frame_number >= 2*particle_number\n");
		exit (1);
	}
	
	if ((total_p*total_t)!=size[1])
		pERROR("bad tracking, there are broken trajectories.");

	// rescale 
	if (a0>0)
	{
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			for (k=0; k<dim; k++)
				*(pointer+k)/=a0;
			pointer+=size[0];
		}
	}
	
	if (Remove_Drift)
		remove_drift(ptid, dim);

	//=====================================================================
	// get mean position
	float * mp_ptr = (float *)calloc(dim*total_p, sizeof(float));
	POINTER_NULL(mp_ptr);

	float *pointer1;
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		for (k=0; k<dim; k++)
			*(pointer1+k)+=*(pointer+k);
		pointer+=size[0];
	}

	// mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
		for (k=0; k<dim; k++)
			*(pointer++)/=total_t;
	
	//=====================================================================
	// get displacement
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(pointer1+k);
		pointer+=size[0];
	}
	free(mp_ptr);

	//=====================================================================
	// read K
	int Kn;
	double *K;
	readK(Kn, &K, file);

	if ( Kn != 2*total_p)
		pERROR("read files error!");
	
	
	//======================================================================
	double *U = (double *) malloc (total_t*sizeof(double));
	POINTER_NULL(U);

	const int nKuu = (total_p*(total_p+1))/2;
	double *Kuu = (double *) malloc (nKuu*sizeof(double));
	POINTER_NULL(Kuu);

	double kuu;
	double *Kuu_ptr=Kuu;
	const int TS=total_t*size[0];

	printf("# calculating the potential energy for each frame ...\n");
	for (t=0; t<total_t; t++)
	{
		pointer=ptid_ptr+t*size[0];
		for (i=0; i<total_p; i++)
		{
			//pointer=ptid_ptr+(i*total_t+t)*size[0];
			pointer1=ptid_ptr+(i*total_t+t)*size[0];
			for (j=i; j<total_p; j++)
			{
				//pointer1=ptid_ptr+(j*total_t+t)*size[0];

				p=2*i*Kn+2*j; // xx
				kuu=K[p]*(*pointer)*(*pointer1);
				++p;         // xy
				//p=2*i*Kn+2*j+1
				kuu+=K[p]*(*pointer)*(*(pointer1+1));
				//p+=Kn-1 - 2*i-1;     // yx
				p+=Kn;     // yy
				//p=(2*i+1)*Kn+Kn+2*j+1
				kuu+=K[p]*(*(pointer+1))*(*(pointer1+1));
				--p;         // yx
				//p=(2*i+1)*Kn+Kn+2*j
				kuu+=K[p]*(*(pointer+1))*(*pointer1);
				if (j==i)
					*(Kuu_ptr++)=kuu/2;
				else
					*(Kuu_ptr++)=kuu;
				pointer1 += TS;
			}
			pointer += TS;
		}
		U[t]=pairwise(nKuu, Kuu)/Kn;  // potential energy per freedom
		Kuu_ptr=Kuu;
	}
	free(K);
	free(Kuu);
	delete [] size;

	char *Ufile=getfilename(file, ".U");
	FILE *Uf=fopen(Ufile, "w");
	FILE_NULL(Uf, Ufile);
	free(Ufile);

	fprintf(Uf, "# potential energy of each frame\n");
	for (t=0; t<total_t; t++)
		fprintf(Uf, "%d %6.6f\n", t, U[t]);
	fclose(Uf);

	//======================================================================
	// distribution of U
	double minU=U[0], maxU=U[0];
	double meanU=pairwise(total_t, U)/total_t;
	double *U2 = (double *) malloc (total_t*sizeof(double));
	POINTER_NULL(U2);
	U2[0]=(U[0]-meanU)*(U[0]-meanU);
	for (t=1; t<total_t; t++)
	{
		if ( minU > U[t] )
			minU = U[t];
		else if ( maxU < U[t] )
			maxU = U[t];
		U2[t]=(U[t]-meanU)*(U[t]-meanU);
	}
	double meanU2=pairwise(total_t, U2)/total_t;
	free(U2);

	const int nbin=100;
	const double binsize = (maxU-minU)*(1+1.0e-6)/nbin;

	unsigned int *hist = (unsigned int *)calloc(nbin, sizeof(unsigned int));
	POINTER_NULL(hist);

	for (t=0; t<total_t; t++)
		++hist[(int)((U[t]-minU)/binsize)];
	free(U);

	char *dUfile=getfilename(file, ".dU");
	FILE *dUf=fopen(Ufile, "w");
	FILE_NULL(dUf, dUfile);
	free(dUfile);

	fprintf(dUf, "# mean potential energy:           %6.6f k_BT\n", meanU);
	fprintf(dUf, "# fluctuation of potential energy: %6.6f k_BT\t"
				 " (per freedom)\n", meanU2);
	fprintf(dUf, "# specific heat (per freedom)    : %6.6f k_B\n", 2*meanU2*Kn);
	fprintf(dUf, "# distribution of potential energy of each frame\n");
	for (i=0; i<nbin; i++)
		if ( hist[i]!=0 )
			fprintf(dUf, "%6.6f %6.6f %u\n", 
					(i+0.5)*binsize+minU, ((double)hist[i])/total_t/binsize,
					hist[i]);
	fclose(dUf);
	free(hist);
}
