#include "image_analysis.h"

#include <iostream>
#include <cmath>

using namespace std;

void bpass(colloid_base& bp, 
		   colloid_base& image, 
		   float& lnoise,
		   float& lobject)
{
	int i,j,t;
	int * size=image.get_size();
	float * pixels=image.get_array_pointer();
	
	int w=(lobject>(2*lnoise)) ? int(lobject+0.5) : int(2*lnoise+0.5);
	int N=2*w+1;

	//float * r=new float [N];
	//for (i=0; i<N; i++)
	//	r[i]=float(i-w)/2/lnoise;

	float * xpt=new float [N];
	float r, total=0;
	for (i=0; i<N; i++)
	{
		r=float(i-w)/2/lnoise;
		xpt[i]=exp(-r*r);
		total+=xpt[i];
	}

	for (i=0; i<N; i++)
		xpt[i]/=total;

	total=0;
	for (int i=0; i<N; i++)
		total+=xpt[i]*xpt[i];
	//float factor=total-1.0/N;

	// convolution
	int col=size[0], row=size[1];

	float * g=new float [size[2]];
	float * b=new float [size[2]];
	for (i=0; i<size[2]; i++)
	{
		g[i]=0; b[i]=0;
	}
	for (j=0; j< row; j++)
		for (t=N/2; t<col-N/2; t++)
		{
			for (i=0; i<N; i++)
				g[j*col+t]+=float(pixels[j*col+t+i-N/2])*xpt[i];
		}
	float *tmp=new float [size[2]];
	for (i=0; i<size[2]; i++)
	{
		tmp[i]=g[i];
		g[i]=0;
	}
	for (j=0; j< col; j++)
		for (t=N/2; t<row-N/2; t++)
		{
			for (i=0; i<N; i++)
				g[j+t*col]+=float(tmp[j+col*(t+i-N/2)])*xpt[i];
		}
	//=======================================================

	for (j=0; j< row; j++)
		for (t=N/2; t<col-N/2; t++)
		{
			for (i=0; i<N; i++)
				b[j*col+t]-=float(pixels[j*col+t+i-N/2])/N;
		}
	for (i=0; i<size[2]; i++)
	{
		tmp[i]=b[i];
		b[i]=0;
	}
	for (j=0; j< col; j++)
		for (t=N/2; t<row-N/2; t++)
		{
			for (i=0; i<N; i++)
				b[j+t*col]-=float(tmp[j+col*(t+i-N/2)])/N;
		}


	delete [] size;
	//delete [] r;
	delete [] xpt;
	delete [] g;
	delete [] b;
	delete [] tmp;
}
