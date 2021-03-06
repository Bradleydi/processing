#include "structure.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void triangle_cluster_plot(const char *file, int frame, bool tracked)
{
	char *gdffile=getfilename(file, "_cluster.gdf");
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
	const int ci=ti-1;

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

	char suffix [15];
	char *epsfile, *tmpfile= tmpnam (NULL);
	FILE *infile, *pipe;
	//for (t=0; t<5;t++)
	for (t=0; t<total_t; t++)
	{
		/*==========================================================*/
		// draw config by Gnuplot
		infile=fopen(tmpfile, "w");
		FILE_NULL(infile, tmpfile);
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			if ( *(pointer+ci) > 0 )
			fprintf(infile, "%f %f %d\n", *pointer, *(pointer+1),
							((int)(*(pointer+ci)))%8+1);
							//nuclear[i-cum_ptcl_frame[t]]);
							//nuclear[(pointer-p_ptr)/size[0]]);
			else if ( *(pointer+ci) < 0 )
			fprintf(infile, "%f %f %d\n", *pointer, *(pointer+1),
							(-(int)(*(pointer+ci)))%8+1);
			else
			fprintf(infile, "%f %f %d\n", *pointer, *(pointer+1),
							0);
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
					  //"splot '%s' using 1:2:(int($3)%%8) with points palette pt 7 ps %.2f\n"
					  "splot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n"
					  "exit\n", epsfile, 
					  xmin-5, xmax+5, ymin-5, ymax+5, tmpfile,
					  50.0/sqrt((float)maxptcl));
		fclose(pipe);
		free(epsfile);

	}
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);

	free(index);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	delete [] size;
}
