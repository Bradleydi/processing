#include "structure.h"
#include "dct.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void triangle(const char *file, double threshold, bool tracked)
{
	if (threshold<=0)
		 pERROR("threshold should be in (0,1).");

	char *gdffile=getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);
	free(gdffile);

	int i, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
		--ti;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");

	float *pointer=p_ptr+ti;
	float t1=*pointer, t2=t1;
	for (i=0; i<size[1]; i++)
	{
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin=(int)t1;
	// reindex time
	pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		*pointer=(float)((int)(*pointer)-tmin);
		pointer+=size[0];
	}

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(ptcl_frame);
	pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	int *cum_ptcl_frame=(int *)malloc((total_t+1)*sizeof(int));
	POINTER_NULL(cum_ptcl_frame);
	cum_ptcl_frame[0]=0;
	cum_ptcl_frame[1]=ptcl_frame[0];
	int maxptcl=ptcl_frame[0];
	for (i=1; i<total_t; i++)
	{
		cum_ptcl_frame[i+1]=cum_ptcl_frame[i]+ptcl_frame[i];
		if ( maxptcl < ptcl_frame[i] )
			maxptcl = ptcl_frame[i];
	}

	if (maxptcl*total_t > 2*size[1])
		pERROR("too many broken trajectories!");

	int *index_now=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(index_now);

	float **index=(float **)malloc(size[1]*sizeof(float *));
	POINTER_NULL(index);
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	float *points=(float *)malloc(2*maxptcl*sizeof(float));
	POINTER_NULL(pointer);
	float *pointer1;
	int ntri;
	/* max_ triangle: 2*maxptcl
	 * so length of ``tri'': 3*max_nedge = 6*maxptcl 
	 * (3 vertices for each triangle)
	 */
	int *tri=(int *)malloc(6*maxptcl*sizeof(int));
	POINTER_NULL(tri);
	//unsigned char *eq_tri=Malloc(unsigned char, 2*maxptcl);
	//POINTER_NULL(eq_tri);
	unsigned char *count=Malloc(unsigned char, maxptcl);
	POINTER_NULL(count);
	float *r=Malloc(float, maxptcl);
	POINTER_NULL(r);
	float **pindex;
	float *A, *B, *C; // triangle vertics A, B, C
	float area, perimeter, Ax, Ay, Bx, By, Cx, Cy, ratio;
	// area/permeter^2 isosceles right triangle
	const float ratio_ri=1.0/2.0/(2.0+sqrt(2.0))/(2.0+sqrt(2.0)) ; 
	// area/permeter^2 equilateral triangle-isosceles right triangle
	const float diff_ratio_eq_ri=sqrt(3.0)/36.0-ratio_ri; 	
	//const float unfold_thres=diff_ratio_eq_ri*((float)threshold)+ratio_ri;
	int j;
	float min=ratio_ri;
	float max=sqrt(3.0)/36.0;
	const int nbin=100;
	const float binsize=(max-min)*(1.0+1.0e-6)/nbin;
	unsigned int *hist = Malloc(unsigned int, nbin);
	POINTER_NULL(hist);
	//char *tmpfile = tmpnam (NULL);
	char *tmpfile=getfilename(file, ".tri");
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);

	//for (t=0; t<1;t++)
	for (t=0; t<total_t; t++)
	{
		pointer1=points;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			*(pointer1++)=*pointer;
			*(pointer1++)=*(pointer+1);
		}
		// get triangles of Delaunay triangluation
		dct_tri(ptcl_frame[t], points, ntri, tri);

		/*==========================================================*/
		// check those triangles
		pindex=index+cum_ptcl_frame[t];
		for (i=0; i<ptcl_frame[t]; i++)
		{
			count[i]=0;
			r[i]=0.0;
		}
		for (i=0; i<ntri; i++)
		{
			A=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex A
			B=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex B
			C=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex C

			Ax=*A; Ay=*(A+1);
			Bx=*B; By=*(B+1);
			Cx=*C; Cy=*(C+1);
			perimeter=sqrt((Ax-Bx)*(Ax-Bx)+(Ay-By)*(Ay-By));
			perimeter+=sqrt((Ax-Cx)*(Ax-Cx)+(Ay-Cy)*(Ay-Cy));
			perimeter+=sqrt((Cx-Bx)*(Cx-Bx)+(Cy-By)*(Cy-By));

			area=((Ax-Cx)*(By-Ay)-(Ax-Bx)*(Cy-Ay))/2;
			if (area<0)
				area=-area;

			ratio=area/perimeter/perimeter;
			j=*(tri-3);
			++count[j];
			r[j]+=ratio;
			j=*(tri-2);
			++count[j];
			r[j]+=ratio;
			j=*(tri-1);
			++count[j];
			r[j]+=ratio;
		}
		tri-=3*ntri;
		for (i=0; i<nbin; i++)
		{
			hist[i]=0;
		}
		for (i=0; i<ptcl_frame[t]; i++)
		{
			r[i]/=count[i];
			if ( r[i] > min )
				++hist[(int)((r[i]-min)/binsize)];
		}
		for (i=0; i<nbin; i++)
			fprintf(infile, "%d %.6f %.6g %u\n", t+tmin,
					(min+(i+0.5)*binsize-ratio_ri)/diff_ratio_eq_ri, 
					((float)hist[i])*diff_ratio_eq_ri/ptcl_frame[t]/binsize,
					hist[i]);
	}

	free(r);
	free(hist);

	fclose(infile);
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}
	fprintf(pipe, "reset\n"
				  "unset key\n"
				  "set pm3d map\n"
				  "set palette rgbformulae 33,13,10\n"
				  "set xrange [%d:%d]\n"
				  "set yrange [0:1]\n"
				  "splot \"%s\" using 1:2:3 with image palette\n"
				  "exit\n", tmin, tmin+total_t-1, tmpfile);
	fclose(pipe);

	free(count);

	free(points);
	free(tri);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	delete [] size;
}
