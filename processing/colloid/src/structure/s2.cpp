#include "structure.h"
#include "colloid_base.h"
#include "miscellaneous.h"
#include "io.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;


void s2(const char *file, const int dim, bool tracked)
{
	char *gdffile=getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);
	free(gdffile);

	const double PI = acos(-1);

	int i, j, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();
	int tindex=size[0]-1;
	if (tracked)
		--tindex;

	float *pointer=p_ptr+tindex;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if (tt < *pointer)
			tt = *pointer;
		pointer += size[0];
	}
	const int total_t=(int)tt+1;

	int *ptcl_frame = (int *) calloc (total_t, sizeof(int));
	POINTER_NULL(ptcl_frame);
	pointer=p_ptr+tindex;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer += size[0];
	}
	
	int *cum_ptcl_frame = (int *) malloc((total_t+1)*sizeof(int));
	POINTER_NULL(cum_ptcl_frame);
	cum_ptcl_frame[0]=0;
	for (i=0; i<total_t; i++)
		cum_ptcl_frame[i+1]=cum_ptcl_frame[i]+ptcl_frame[i];

	float **newt = (float **)malloc(size[1]*sizeof(float *));
	POINTER_NULL(newt);
	int *index_now = (int *)calloc(total_t, sizeof(int));
	POINTER_NULL(index_now);
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+tindex));
		newt[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer += size[0];
	}
	free(index_now);
	free(cum_ptcl_frame);

	//=======================================================
	pointer=p_ptr;
	float xmin=*pointer, xmax=*pointer;
	float ymin=*(pointer+1), ymax=*(pointer+1);
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if (xmin > *pointer)
			xmin=*pointer;
		else if (xmax < *pointer)
			xmax=*pointer;
		
		if (ymin > *(pointer+1))
			ymin=*(pointer+1);
		else if (ymax < *(pointer+1))
			ymax=*(pointer+1);
	}
	double N=double(size[1])/total_t;
	double w=xmax-xmin, h=ymax-ymin;
	double V=w*h;
	double density=N/V;
	//=======================================================

	/*
	float rmin=0, rmax=30;
	int nbin=300;
	float binsize=(rmax-rmin)/nbin;
	*/
	//float rmin=0; 
	float rmax=( w < h ) ? (float)w/2 : (float)h/2;
	//if (rmax>150)
	//	rmax=150;
	float binsize=0.1;  // resolution of particle track technic
	//const int nbin=(int)((rmax-rmin)/binsize);
	const int nbin=(int)(rmax/binsize);
	rmax=nbin*binsize-0.01; // the last bin may be inaccurate, but others
	                        // are OK.
	float rmax2=rmax*rmax;

	//unsigned long * g = new unsigned long [nbin];
	double * g = (double *)calloc(nbin, sizeof(double));
	POINTER_NULL(g);

	float *pointer1, *pointer2;
	float dx, dy, d;
	double area;
	float **newt_pointer1=newt, **newt_pointer2;
	for (t=0; t<total_t; t++)
	//for (t=0; t<tmax; t++)
	{
		//newt_pointer=newt+cum_ptcl_frame[t];
		//newt_pointer1=newt+cum_ptcl_frame[t];
		for (i=0; i<ptcl_frame[t]; i++)
		{
			//pointer1=newt[cum_ptcl_frame[t]+i];
			//pointer1=*(newt_pointer+i);
			pointer1=*(newt_pointer1++);     // i
			newt_pointer2=newt_pointer1;     // i+1
			for (j=i+1; j<ptcl_frame[t]; j++)
			{
				//pointer2=newt[cum_ptcl_frame[t]+j];
				//pointer2=*(newt_pointer+j);
			
				pointer2=*(newt_pointer2++);     // j
				
				dx=*pointer1-*pointer2;
				d=dx*dx;
				if ( d < rmax2 )
				{
					dy=*(pointer1+1)-*(pointer2+1);
					d+=dy*dy;
					if ( d < rmax2 )
					{
						if (dx<0)
							dx=-dx;
						if (dy<0)
							dy=-dy;

						area=(w-dx)*(h-dy);
						
						g[(int)(sqrt(d)/binsize)]+=2/density/area;
					}
				}
				/*
				dx=*pointer1-*pointer2;
				dy=*(pointer1+1)-*(pointer2+1);
				d=dx*dx+dy*dy;
				d=sqrt(d);
				*/

				/*
				//k=(int)((d-rmin)/binsize);
				k=(int)(d/binsize);
				if (k<nbin) 
				{*/
				/*
				if ( d < rmax )
				{
					if (dx<0)
						dx=-dx;
					if (dy<0)
						dy=-dy;

					area=(w-dx)*(h-dy);
					
					g[(int)(d/binsize)]+=2/density/area;
				}
				*/
			}
		}
	}

	free(ptcl_frame);
	//free(cum_ptcl_frame);
	free(newt);


	double *s=Malloc(double, nbin); POINTER_NULL(s);
	for (i=0; i<nbin; i++)
	{
		d=(i+0.5)*binsize;
		//cout << d << ' ' << g[i]/PI/2/d/binsize/density/total_t << endl;
		//printf("%6.6f %6.6f\n", d, g[i]/PI/2/d/binsize/density/total_t);
		g[i]/=PI*2.0*d*binsize*density*total_t;
		if ( dim == 2 )
		{
			if ( g[i] <= 0.0 )
				s[i]=d*PI*density*binsize;
			else
				s[i]=(g[i]*log(g[i])-g[i]+1.0)*d*PI*density*binsize;
		}
		else if ( dim == 3 )
		{
			if ( g[i] <= 0.0 )
			s[i]=d*d*4.0*PI*density*binsize;
			else
			s[i]=(g[i]*log(g[i])-g[i]+1.0)*d*d*4.0*PI*density*binsize;
		}
	}
	// find a0
	for ( i=10; i<nbin; i++) // search from 1 pixel
	{
		if ( g[i] > 1.0 )
		{
			for (j=-10; j<=10; j++)
				if ( g[i+j] > g[i] )
					break;
			if ( j==11 )
				break;
		}
	}
	double a0=g[i];

	char *filename = getfilename (file, ".s2");
	FILE *infile = fopen(filename, "w");
	FILE_NULL(infile, filename);
	free(filename);
	fprintf(infile, "# a0 = %6.6f\n", a0);
	if ( dim == 2 )
		fprintf(infile, "# exess entropy s2 = -%6.6f\n", pairwise(nbin, s));
	else if ( dim == 3 )
		fprintf(infile, "# exess entropy s2 = -%6.6f\n", pairwise(nbin, s)/a0);
	fprintf(infile, "# r s2 g\n");
	double s2=0.0;
	for (i=0; i<nbin; i++)
	{
		d=(i+0.5)*binsize;
		s2 += s[i];
		if ( dim == 2)
			fprintf(infile, "%6.6f %6.6f %6.6f\n", d, -s2, g[i]);
		else if ( dim == 3 )
			fprintf(infile, "%6.6f %6.6f %6.6f\n", d, -s2/a0, g[i]);
	}

	free(g);
	free(s);
	delete [] size;
}
