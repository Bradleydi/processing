#include "normalmode.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;

void uur(colloid_base& ptid, const char *file, double a0, bool Remove_Drift)
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

	// frame number per particle
	int *frame_ptcl =(int *)calloc(total_p, sizeof(int));
	POINTER_NULL(frame_ptcl);

	pointer=ptid_ptr+ti;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		
		++frame_ptcl[(int)(*(pointer+1))];
		pointer+=size[0];
	}
	const int total_t=(int)tt+1;
	
	if ((total_p*total_t)!=size[1])
		pERROR("bad tracking, there are broken trajectories.");

	printf("# total particle\t\t:\t%d\n", total_p);
	printf("# total frame   \t\t:\t%d\n", total_t);

	if ( total_t < 2*total_p )
	{
		fprintf(stderr, "# Error: there should be enough independent"
				        " configurations!\n"
						"#        frame_number >= 2*particle_number\n");
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
	// get covariance
	printf("# getting covariance ...\n");
	double *CXX=Malloc(double, total_t); POINTER_NULL(CXX);
	double *CXY=Malloc(double, total_t); POINTER_NULL(CXY);
	double *CYX=Malloc(double, total_t); POINTER_NULL(CYX);
	double *CYY=Malloc(double, total_t); POINTER_NULL(CYY);

	float xmin=mp_ptr[0], xmax=mp_ptr[0];
	float ymin=mp_ptr[1], ymax=mp_ptr[1];
	pointer=mp_ptr+dim;
	for (i=1; i<total_p; i++)
	{
		if ( xmin > *pointer)
			xmin = *pointer;
		else if ( xmax < *pointer)
			xmax = *pointer;
		++pointer;

		if ( ymin > *pointer)
			ymin = *pointer;
		else if ( ymax < *pointer)
			ymax = *pointer;
		++pointer;
	}

	float dx, dy, d;
	const int nbin=100;
	const double binsize=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin))/2/nbin;

	double *uurxx=(double *)calloc(nbin, sizeof(double));
	POINTER_NULL(uurxx);
	double *uurxy=(double *)calloc(nbin, sizeof(double));
	POINTER_NULL(uurxy);
	double *uuryx=(double *)calloc(nbin, sizeof(double));
	POINTER_NULL(uuryx);
	double *uuryy=(double *)calloc(nbin, sizeof(double));
	POINTER_NULL(uuryy);

	unsigned int *counter=(unsigned int *)calloc(nbin, sizeof(unsigned int));
	POINTER_NULL(counter);

	int itt, jtt;
	for (i=0; i<total_p; i++)
	{
		itt=i*total_t;
		for (j=i+1; j<total_p; j++) // not accout i itself
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
			dx=mp_ptr[j*dim]-mp_ptr[i*dim];
			dy=mp_ptr[j*dim+1]-mp_ptr[i*dim+1];
			d=sqrt(dx*dx+dy*dy);

			k=(int)(d/binsize);
			if (k<nbin)
			{
				uurxx[k]+=pairwise(total_t, CXX)/total_t;
				uurxy[k]+=pairwise(total_t, CXY)/total_t;
				uuryx[k]+=pairwise(total_t, CYX)/total_t;
				uuryy[k]+=pairwise(total_t, CYY)/total_t;
				++counter[k];
			}
		}
	}

	// i itself (r=0)
	double C0=0;
	for (i=0; i<total_p; i++)
	{
		itt=i*total_t;
		for (t=0; t<total_t; t++)
		{
			pointer=index[itt+t];
			CXX[t]=(double)((*pointer)*(*pointer));
			CYY[t]=(double)((*(pointer+1))*(*(pointer+1)));
		}
		C0+=pairwise(total_t, CXX)+pairwise(total_t, CYY);
	}
	C0 /= total_p*total_t;

	free(mp_ptr);
	free(index);
	free(CXX);
	free(CXY);
	free(CYX);
	free(CYY);
	
	char *filename=getfilename(file, ".uur");
	FILE *uurf=fopen(filename, "w");
	FILE_NULL(uurf, filename);

	fprintf(uurf, "# <u(r')*u(r'+r)> vs. distance\n");
	fprintf(uurf, "# r uxux uxuy uyux uyuy (uxux+uyuy)/(uxux+uyuy)(r=0) statistics\n");
	fprintf(uurf, "# <u(r)^2> = %f\n", C0);

	int min_counter = total_t/10;
	for (i=0; i<nbin; i++)
		if (counter[i] >= min_counter )
		{
			uurxx[i]/=counter[i];
			uurxy[i]/=counter[i];
			uuryx[i]/=counter[i];
			uuryy[i]/=counter[i];
			fprintf(uurf, "%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%d\n",
					(i+0.5)*binsize, uurxx[i], uurxy[i],
					uuryx[i], uuryy[i], (uurxx[i]+uuryy[i])/C0, counter[i]);
		}
	fclose(uurf);

	free(uurxx);
	free(uurxy);
	free(uuryx);
	free(uuryy);
	free(counter);
	free(filename);
	delete [] size;
}
