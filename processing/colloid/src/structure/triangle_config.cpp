#include "structure.h"
#include "dct.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void triangle_config(const char *file, double threshold, bool tracked)
{
	if ( threshold<=0 || threshold>=1.0 )
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

	float *pointer=p_ptr;
	float xmin=*pointer, xmax=xmin, ymin=*(pointer+1), ymax=ymin;
	float t1=*(pointer+ti), t2=t1;
	for (i=1; i<size[1]; i++)
	{
		if ( xmin > *pointer )
			xmin = *pointer;
		else if ( xmax < *pointer )
			xmax = *pointer;
		
		if ( ymin > *(pointer+1) )
			ymin = *(pointer+1);
		else if ( ymax < *(pointer+1) )
			ymax = *(pointer+1);

		if ( t1 > *(pointer+ti) )
			t1 = *(pointer+ti);
		else if ( t2 < *(pointer+ti) )
			t2 = *(pointer+ti);
		pointer+=size[0];
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin = (int)t1;
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
	unsigned char *nuclear=Malloc(unsigned char, maxptcl);
	POINTER_NULL(nuclear);
	float **pindex;
	float *A, *B, *C; // triangle vertics A, B, C
	float area, perimeter, Ax, Ay, Bx, By, Cx, Cy, ratio;
	// area/permeter^2 isosceles right triangle
	const float ratio_ri=1.0/2.0/(2.0+sqrt(2.0))/(2.0+sqrt(2.0)) ; 
	// area/permeter^2 equilateral triangle-isosceles right triangle
	const float diff_ratio_eq_ri=sqrt(3.0)/36.0-ratio_ri; 	
	const float unfold_thres=diff_ratio_eq_ri*((float)threshold)+ratio_ri;

	char suffix [15];
	char *epsfile, *tmpfile= tmpnam (NULL);
	FILE *infile, *pipe;
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
			nuclear[i]=0;
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

			if ( area/perimeter/perimeter > unfold_thres )
			{
				// is an equilateral triangle
				nuclear[*(tri-3)]=1;
				nuclear[*(tri-2)]=1;
				nuclear[*(tri-1)]=1;
				/*
				nuclear[(A-p_ptr)/size[0]]=1;
				nuclear[(B-p_ptr)/size[0]]=1;
				nuclear[(C-p_ptr)/size[0]]=1;
				*/
			}
		}
		tri-=3*ntri;

		/*==========================================================*/
		// draw config by Xfig
		
		/*==========================================================*/
		// draw config by Gnuplot
		infile=fopen(tmpfile, "w");
		FILE_NULL(infile, tmpfile);
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			fprintf(infile, "%f %f %u\n", *pointer, *(pointer+1),
							nuclear[i-cum_ptcl_frame[t]]);
							//nuclear[(pointer-p_ptr)/size[0]]);
		}
		fclose(infile);
		pipe=popen("gnuplot -persist", "w");
		if (pipe==NULL)
		{
			fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
			exit (1);
		}
		sprintf(suffix, "_%05d.eps", t+tmin);
		epsfile=getfilename(file, suffix);
		fprintf(pipe, "reset\n"
					  "set terminal postscript eps color enhanced\n"
					  "set output '%s'\n"
					  "unset key\n"
					  "unset tics\n"
					  "set size ratio -1\n"
					  "set xrange [%f:%f]\n"
					  "set yrange [%f:%f]\n"
					  "set pm3d map\n"
					  "set palette rgbformulae 33,13,10\n"
					  "unset colorbox\n"
					  "splot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n"
					  "exit\n", epsfile, xmin, xmax, ymin, ymax, tmpfile,
					  50.0/sqrt((float)maxptcl));
		fclose(pipe);
		free(epsfile);

	}
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);

	free(nuclear);
	free(points);
	free(tri);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	delete [] size;
}
