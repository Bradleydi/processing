#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

double * points2image(int &nxy, double *xydata, 
					  int &w, int &h, 
					  double *values, bool exact)
{
	// max and min of xydata
	double *pointer=xydata;
	double xmin=*(pointer++);
	double ymin=*(pointer++);
	int i;
	for (i=0; i<nxy; i++)
	{
		if ( xmin > *pointer )
			xmin = *pointer;

		++pointer;
		
		if ( ymin > *pointer )
			ymin = *pointer;

		++pointer;
	}

	pointer=xydata;
	*pointer-=xmin;
	double xmax=*(pointer++);
	*pointer-=ymin;
	double ymax=*(pointer++);

	for (i=0; i<nxy; i++)
	{
		*pointer-=xmin;
		if ( xmax < *pointer )
			xmax = *pointer;

		++pointer;
		
		*pointer-=ymin;
		if ( ymax < *pointer )
			ymax = *pointer;

		++pointer;
	}

	int Ixmax=(int)xmax;
	int Iymax=(int)ymax;

	if (w<=0 || w>Ixmax)
		w=Ixmax+1;
	if (h<=0 || h>Iymax)
		h=Iymax+1;


	char *good=(char *)calloc(nxy, sizeof(char));
	int count=0;
	pointer=xydata;
	for (i=0; i<nxy; i++)
	{
		if ( (int)(*(pointer++)) >= w )
		{
			good[i]=1;
			++count;
		}
		if ( (int)(*(pointer++)) >= h )
		{
			good[i]=1;
			++count;
		}
	}
	
	if (count >= nxy)
		pERROR("setted w and h is not valid.");

	int xi, yi, index;
	double xf, yf, xfyf;
	double *image=(double *)calloc(w*h, sizeof(double));
	pointer=xydata;
	if (values==NULL)
	{
		if (exact)
		{
			for (i=0; i<nxy; i++)
			{
				if (good[i]!=0)
				{
					xi=(int)*pointer;
					xf=*(pointer++)-xi;
					yi=(int)*pointer;
					yf=*(pointer++)-yi;
					xfyf=xf*yf;

					index=xi+yi*w;
					image[index]+=1.0-xf-yf-xfyf;
					image[index+1]+=xf-xfyf;
					image[index+w]+=yf-xfyf;
					image[index+w+1]+=xfyf;
				}
			}
		}
		else
		{
			for (i=0; i<nxy; i++)
			{
				if (good[i]!=0)
				{
					xi=(int)(*(pointer++));
					yi=(int)(*(pointer++));

					image[xi+yi*w]+=1.0;
				}
			}
		}
	}
	else
	{
		if (exact)
		{
			for (i=0; i<nxy; i++)
			{
				if (good[i]!=0)
				{
					xi=(int)*pointer;
					xf=*(pointer++)-xi;
					yi=(int)*pointer;
					yf=*(pointer++)-yi;
					xfyf=xf*yf;

					index=xi+yi*w;
					image[index]+=values[i]*(1.0-xf-yf-xfyf);
					image[index+1]+=values[i]*(xf-xfyf);
					image[index+w]+=values[i]*(yf-xfyf);
					image[index+w+1]+=values[i]*xfyf;
				}
			}
		}
		else
		{
			for (i=0; i<nxy; i++)
			{
				if (good[i]!=0)
				{
					xi=(int)(*(pointer++));
					yi=(int)(*(pointer++));

					image[xi+yi*w]+=values[i];
				}
			}
		}
	}
	free(good);
	return image;
}
