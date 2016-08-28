#include "plot.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "unistd.h"

// this is a virtual trajectory plot tool, since actually the
// program doesn't need the data being tracked.
//
// 
// In this program, we choose binary image format
void traj_plot(const char *file, int twindow, int tstep,
				 const char *str_frameindex, 
				 int framerate, bool Tracked)
{
	if ( twindow<=0 )
	{
		twindow = 10;
		printf("# set twindow to default value: 10\n");
	}
	if ( tstep<=0 || tstep>=twindow )
	{
		tstep = 1;
		printf("# set tstep to default value: 1\n");
	}
	if (framerate <=0 )
	{
		framerate=2;
		printf("# set frame rate to default value: 2fpt\n");
	}
	char *gdffile=getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);
	free(gdffile);

	int i, t, k;
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
	if ( twindow>=total_t )
	{
		twindow = 10;
		printf("# set twindow to default value: 10\n");
	}
	
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
	// binary image set
	//==================================================================
	const int ixmax=(int)(xmax+1.5);
	const int ixmin=(int)(xmin-0.5);
	const int iymax=(int)(ymax+1.5);
	const int iymin=(int)(ymin-0.5);

	const int col = (ixmax-ixmin);
	const int row = (iymax-iymin);
	const int image_size = col*row;
	const int offset = ixmin+iymin*col;

	unsigned char *image=Malloc(unsigned char, image_size);
	POINTER_NULL(image);

	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer+=size[0];
	}

	//==================================================================
	// plot
	//==================================================================

	char *tmpfile= tmpnam (NULL);
	FILE *infile, *pipe;
	//const float wxtpointsize=20.0/sqrt((float)maxptcl);
	//const float x11pointsize=70.0/sqrt((float)maxptcl);
	const int sleeptime=1000000/framerate;
	unsigned char color;
	
	//==================================================================
	// open pipe to Gnuplot
	//==================================================================
	pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}
	printf("# time window :\t%d\n step :\t\t%d\n", twindow, tstep);
	fprintf(pipe, //"set terminal postscript eps color enhanced\n"
				  //"set output '%s'\n"
				  "unset key\n"
				  //"unset tics\n"
				  "set size ratio -1\n"
				  "set xrange [%d:%d]\n"
				  "set yrange [%d:%d]\n"
				  //"if ( GPVAL_TERM eq 'wxt' ) set style line 1 pt 7 ps %f;"
				  //"else if ( GPVAL_TERM eq 'x11' ) "
				  //"set style line 1 pt 7 ps %f;\n"
				  "set cbrange [0:%d]\n"
				  "unset colorbox\n"
				  "set palette rgbformulae 33,13,10\n"
				  "set pm3d map\n"
				  , ixmin, ixmax, iymin, iymax, 
				  //wxtpointsize, x11pointsize, 
				  (int)ceil(((float)twindow)/tstep));

	for (t=frameindex[0]; t<=frameindex[1];t+=frameindex[2])
	{
		//==========================================================
		// save data
		//==========================================================
		color=1u;
		for (i=0; i<image_size; i++) image[i]=0u;
		for (k=0; k<twindow; k+=tstep)
		{
			for (i=cum_ptcl_frame[t+k]; i<cum_ptcl_frame[t+k+1]; i++)
			{
				pointer=index[i];
				*(image+(int)(*pointer+0.5)+((int)(*(pointer+1)+0.5))*col-offset)=color;
				//*(image+((int)(*pointer+0.5)+
				//((int)(*(pointer+1)+0.5))*col-offset)*factor)=1; 
				//// since I didn't use factor so far, so I remove it.
			}
			++color;// here it's a bug
		}
		infile=fopen(tmpfile, "wb");
		FILE_NULL(infile, tmpfile);
		if ( (int)fwrite(image, sizeof(image[0]), image_size, infile) !=
				image_size )
			pERROR("save image data failed");
		fclose(infile);

		printf("# Frame: %6d - %6d\n", t+tmin, t+k+tmin);
		fflush(stdout);
		fprintf(pipe, "splot '%s' binary array=%dx%d format=\"%%uchar\" with image\n", tmpfile, col, row);
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
