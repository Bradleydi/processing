#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;


void covariance(colloid_base& ptid, const char *file, bool Remove_Drift)
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

	// total particle number
	// = last element of ptid+1
	const int total_p=(int)ptid_ptr[size[2]-1]+1;

	// frame number per particle
	int *frame_ptcl =(int *)calloc(total_p, sizeof(int));
	if (frame_ptcl==NULL)
		pERROR("int *frame_ptcl initialized error!");


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
	
	int *ptcl_frame;
	float *drift;
	if (Remove_Drift)
	{
		//particle number per frame
		ptcl_frame=(int *)calloc(total_t, sizeof(int));
		if (ptcl_frame==NULL)
			pERROR("int *ptcl_frame initialized error!");

		pointer=ptid_ptr+ti;
		for (i=0; i<size[1]; i++)
		{
			++ptcl_frame[(int)(*pointer)];
			pointer+=size[0];
		}

		// get drift sum
		drift=(float *)calloc(dim*total_t, sizeof(float));
		if (drift==NULL)
			pERROR("float *drift initialized error!");
	
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			t=(int)(*(pointer+ti));
			pointer1=drift+t*dim;
			for (k=0; k<dim; k++)
				*(pointer1+k)+=*(pointer+k);
			pointer+=size[0];
		}

		// average
		pointer=drift;
		for (t=0; t<total_t; t++)
		{
			for (k=0; k<dim; k++)
				*(pointer++)/=ptcl_frame[t];
		}

		// remove drift
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			t=(int)(*(pointer+ti));
			pointer1=drift+dim*t;
			for (k=0; k<dim; k++)
				*(pointer+k)-=*(pointer1+k);
			pointer+=size[0];
		}

		free(drift);
		free(ptcl_frame);
	}

	//=====================================================================
	// get mean position
	
	float * mp_ptr = (float *)calloc(dim*total_p, sizeof(float));
	if (mp_ptr==NULL)
		pERROR("float *mp_ptr initialized error!");

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
			*(pointer++)/=frame_ptcl[p];

	
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


	free(frame_ptcl);
	

		//printf("test\n");
	//======================================================================
	// get index
	float **index=(float **)malloc(total_t*total_p*sizeof(float *));
	if (index==NULL)
		pERROR("float **index initialized error!");

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
	printf("# getting covariance ...\n");
	double *CXX=(double *)malloc(total_t*sizeof(double));
	double *CXY=(double *)malloc(total_t*sizeof(double));
	double *CYX=(double *)malloc(total_t*sizeof(double));
	double *CYY=(double *)malloc(total_t*sizeof(double));
	
	double *CXX2=(double *)malloc(total_t*sizeof(double));
	double *CXY2=(double *)malloc(total_t*sizeof(double));
	double *CYX2=(double *)malloc(total_t*sizeof(double));
	double *CYY2=(double *)malloc(total_t*sizeof(double));
	double meanXX, meanXY, meanYX, meanYY, d, dd;
	if (CXX==NULL)
		pERROR("double *CXX initialized error!");
	if (CXY==NULL)
		pERROR("double *CXY initialized error!");
	if (CYX==NULL)
		pERROR("double *CYX initialized error!");
	if (CYY==NULL)
		pERROR("double *CYY initialized error!");

	if (CXX2==NULL)
		pERROR("double *CXX2 initialized error!");
	if (CXY2==NULL)
		pERROR("double *CXY2 initialized error!");
	if (CYX2==NULL)
		pERROR("double *CYX2 initialized error!");
	if (CYY2==NULL)
		pERROR("double *CYY2 initialized error!");

	const char subfix[]=".cv";
	char *filename=getfilename(file, subfix);

	FILE *cvf=fopen(filename, "w");
	if (cvf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}

	fprintf(cvf, "# covariance vs. distance\n");
	fprintf(cvf, "# d CXX CXY CYX CYY SDXX SDXY SDYX SDYY total_t\n");

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
			dd=mp_ptr[j*dim]-mp_ptr[i*dim];
			d=dd*dd;
			dd=mp_ptr[j*dim+1]-mp_ptr[i*dim+1];
			d+=dd*dd;
			d=sqrt(d);

			meanXX=pairwise(total_t, CXX)/total_t;
			meanXY=pairwise(total_t, CXY)/total_t;
			meanYX=pairwise(total_t, CYX)/total_t;
			meanYY=pairwise(total_t, CYY)/total_t;
			for (t=0; t<total_t; t++)
			{
				CXX2[t]=CXX[t]-meanXX;
				CXX2[t]=CXX2[t]*CXX2[t];

				CXY2[t]=CXY[t]-meanXY;
				CXY2[t]=CXY2[t]*CXY2[t];
				
				CYX2[t]=CYX[t]-meanYX;
				CYX2[t]=CYX2[t]*CYX2[t];

				CYY2[t]=CYY[t]-meanYY;
				CYY2[t]=CYY2[t]*CYY2[t];
			}
			fprintf(cvf, 
					"%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%6.6f\t%d\n", 
					d, meanXX, meanXY, meanYX, meanYY, 
					pairwise(total_t, CXX2)/total_t,
					pairwise(total_t, CXY2)/total_t,
					pairwise(total_t, CYX2)/total_t,
					pairwise(total_t, CYY2)/total_t,
					total_t);
			
		}
	}
	fclose(cvf);

	free(index);
	free(CXX);
	free(CXY);
	free(CYX);
	free(CYY);
	free(CXX2);
	free(CXY2);
	free(CYX2);
	free(CYY2);
	free(mp_ptr);
	free(filename);

	delete [] size;
}
