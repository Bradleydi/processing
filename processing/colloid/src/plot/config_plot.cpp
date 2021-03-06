#include "plot.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "unistd.h"

void config_plot(const char *file, const char *str_frameindex, 
				 int framerate, bool Tracked)
{
	if (framerate <=0 )
	{
		framerate=2;
		printf("set frame rate to default value: 2fpt\n");
	}
	char *gdffile=getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);
	free(gdffile);

	int i, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	const int dim=2;
	int ti=size[0]-1;
	if (Tracked)
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
	
	int* frameindex=getframeindex(str_frameindex, tmin, (int)t2);
	frameindex[0]-=tmin;
	frameindex[1]-=tmin;
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

	//==================================================================
	// plot
	//==================================================================

	char *tmpfile= tmpnam (NULL);
	FILE *infile, *pipe;
	const float pointsize=50.0/sqrt((float)maxptcl);
	const int sleeptime=1000000/framerate;
	
	//==================================================================
	// open pipe to Gnuplot
	//==================================================================
	pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}
	float Xmargin=(xmax-xmin)*0.05, Ymargin=(ymax-ymin)*0.05;
	fprintf(pipe, //"set terminal postscript eps color enhanced\n"
				  //"set output '%s'\n"
				  "unset key\n"
				  //"unset tics\n"
				  "set size ratio -1\n"
				  "set xrange [%f:%f]\n"
				  "set yrange [%f:%f]\n"
				  //"set cbrange [%f:%f]\n"
				  //"set cbtics\n"
				  //"set pm3d map\n"
				  //"set palette rgbformulae 33,13,10\n"
				  //"unset colorbox\n"
				  ,//epsfile, 
				  xmin-Xmargin, xmax+Xmargin, 
				  ymin-Ymargin, ymax+Ymargin);

	for (t=frameindex[0]; t<=frameindex[1];t+=frameindex[2])
	{
		/*==========================================================*/
		// draw config by Gnuplot
		infile=fopen(tmpfile, "w");
		FILE_NULL(infile, tmpfile);
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			fprintf(infile, "%f %f\n", *pointer, *(pointer+1));
		}
		fclose(infile);

		printf("Frame: %d\n", t);
		fprintf(pipe, "plot '%s' using 1:2 with points pt 7 ps %.2f\n",
					  tmpfile, pointsize);
		fflush(pipe);  // to put them into pipe and leave buffer
		usleep(sleeptime);
		//sleep(2);
	}
	fprintf(pipe, "exit\n");
	fclose(pipe);
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);

	free(index);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	free(frameindex);
	delete [] size;
}
