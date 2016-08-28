#include "normalmode.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;


void normal_mode(colloid_base& ptid, const char *file, double a0,
		bool Remove_Drift)
{
	
	int i, j, k, t, p;
	int *size=ptid.get_size();
	const int dim=2; // dimension
	if (dim+2>size[0])
		pERROR("Not a 2D data or tracked data!");
	
	const int ti=size[0]-2; //time index in the array
	const int ii=ti+1; // particle index in the array
	
	float *ptid_ptr=ptid.get_array_pointer();
	float *pointer, *pointer1;

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

	// total particle number
	// = last element of ptid+1
	const int total_p=(int)ptid_ptr[size[2]-1]+1;

	pointer=ptid_ptr+ti;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t=(int)tt+1;
	
	printf("# total particle\t\t:\t%d\n", total_p);
	printf("# total frame   \t\t:\t%d\n", total_t);

	if ( total_t < 2*total_p )
	{
		fprintf(stderr, "# Error: there should be enough independent configurations!\n");
		fprintf(stderr, "#        frame_number >= 2*particle_number\n");
		exit (1);
	}
	
	if ((total_p*total_t)!=size[1])
		pERROR("bad tracking, there are broken trajectories.");

	
	if (Remove_Drift)
		remove_drift(ptid, dim);

	//=====================================================================
	// get mean position
	
	float * mp_ptr = (float *)calloc(dim*total_p, sizeof(float));
	POINTER_NULL(mp_ptr);

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
			*(pointer++)/=total_t;  //frame_ptcl[p];

	
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

	writemp(total_p, mp_ptr, file);
	free(mp_ptr);

	//======================================================================
	// get index
	float **index=(float **)malloc(total_t*total_p*sizeof(float *));
	POINTER_NULL(index);

	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		p=(int)(*(pointer+ii));
		index[t+p*total_t]=pointer;
		pointer+=size[0];
	}

	//======================================================================
	// get covariance matrix
	printf("# getting covariance matrix...\n");
	int Gn=total_p*2;
	double *G=(double *)malloc(Gn*Gn*sizeof(double));
	double *CXX=(double *)malloc(total_t*sizeof(double));
	double *CXY=(double *)malloc(total_t*sizeof(double));
	double *CYX=(double *)malloc(total_t*sizeof(double));
	double *CYY=(double *)malloc(total_t*sizeof(double));
	if (G==NULL)
		pERROR("double *G initialized error!");
	if (CXX==NULL)
		pERROR("double *CXX initialized error!");
	if (CXY==NULL)
		pERROR("double *CXY initialized error!");
	if (CYX==NULL)
		pERROR("double *CYX initialized error!");
	if (CYY==NULL)
		pERROR("double *CYY initialized error!");

	int itt, jtt;
	for (i=0; i<total_p; i++)
	{
		itt=i*total_t;
		for (j=i; j<total_p; j++)
		{
			jtt=j*total_t;
			for (t=0; t<total_t; t++)
			{
				pointer=index[itt+t];
				pointer1=index[jtt+t];

				CXX[t]=(double)((*pointer)*(*pointer1));
				CXY[t]=(double)((*pointer)*(*(pointer1+1)));
				CYX[t]=(double)((*(pointer+1))*(*pointer1));
				CYY[t]=(double)((*(pointer+1))*(*(pointer1+1)));
			}
			G[2*i*Gn+2*j]=pairwise(total_t, CXX)/total_t;
			G[(2*i+1)*Gn+2*j+1]=pairwise(total_t, CYY)/total_t;
			G[2*i*Gn+2*j+1]=pairwise(total_t, CXY)/total_t;
			if(j!=i)
				G[(2*i+1)*Gn+2*j]=pairwise(total_t, CYX)/total_t;
		}
	}

	for(i=0; i<Gn; i++)
		for (j=0; j<i; j++)
			G[i*Gn+j]=G[j*Gn+i];

	free(index);
	free(CXX);
	free(CXY);
	free(CYX);
	free(CYY);

	ptid.free_memory();
	
	//printf("# DONE!\n");
	
	// write G to file.dcm
	writedcm(Gn, G, file);

	double *E=(double *)malloc(Gn*sizeof(double));
	if (E==NULL)
		pERROR("double *E initialized error!");
	printf("# Calculating eigenvalues...\n");
	dsyev(G, Gn, E);

	if (E[0]<0)
	{
		printf("# First eigenvalue: %6.10f\n", E[0]);
		pERROR("negative eigenvalue(s)");
	}

	// write E and G to file.ev
	writeev(Gn, E, G, file);

	participation_ratio_w(Gn, E, G, file);
	free(G);

	DOS_w(Gn, E, file);
	level_separation_w(Gn, E, file);

	free(E);

	delete [] size;
}
