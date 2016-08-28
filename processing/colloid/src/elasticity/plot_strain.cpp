#include "elasticity.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void plot_strain(const char *file, const char *str_frameindex, 
				const char *which_strain, int framerate)
{
	if (framerate <=0 )
	{
		framerate=2;
		printf("set frame rate to default value: 2fpt\n");
	}

	int WS=0;
	if (which_strain[0]=='-')
	{
		if (which_strain[1]=='b')
			WS=0;
		else if (which_strain[1]=='s')
			WS=1;
		else if (which_strain[1]=='o')
			WS=2;
		else
		pERROR("Wrong strain specified:\n\t-b\tbulk strain"
				"\n\t-s\tshear strain\n\t-o\tlocal rotation");
	}
	else
		pERROR("Wrong strain specified:\n\t-b\tbulk strain"
				"\n\t-s\tshear strain\n\t-o\tlocal rotation");

	char *gdffile=getfilename(file, "_stn.gdf");
	colloid_base stn, ptid;
	readgdf(stn, gdffile);
	gdffile=getfilename(file, ".gdf");
	readgdf(ptid, gdffile);
	free(gdffile);

	int i, t;
	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();

	const int dim=2;
	int ti=size[0]-2;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");

	float *pointer=ptid_ptr;
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
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		*pointer=(float)((int)(*pointer)-tmin);
		pointer+=size[0];
	}

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(ptcl_frame);
	pointer=ptid_ptr+ti;
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
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	char *tmpfile= tmpnam (NULL);
	FILE *infile, *pipe;
	float *stn_ptr=stn.get_array_pointer();
	int *stn_size=stn.get_size();
	float *spointer=stn_ptr;

	float *sd=Calloc(float, 3); POINTER_NULL(sd);
	int CC=0;
	for (i=0; i<size[1]; i++)
	{
		if (*spointer > 1.0e5)
			continue;

		t2=*spointer+*(spointer+3);  //bulk strain
		sd[0]+=t2*t2;

		t2=*(spointer+1) + ( *(spointer+2) ); // shear strain
		sd[1]+=t2*t2;

		t2=*(spointer+2) - ( *(spointer+1) ); // local rotation
		sd[2]+=t2*t2;

		++CC;
		spointer+=stn_size[0];
	}
	sd[1]/=4.0; sd[2]/=4.0;
	for (i=0; i<3; i++)
	{
		sd[i]=2*sqrt(sd[i]/CC);
		printf("%f ", sd[i]);
	}
	putchar('\n');

	pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}
	fprintf(pipe, "reset\n");
	if (WS==0)
		fprintf(pipe, "set title 'Bulk strain'\n");
	else if (WS==1)
		fprintf(pipe, "set title 'Shear strain'\n");
	else if (WS==2)
		fprintf(pipe, "set title 'Local rotation'\n");
	else
		pERROR("unknown strain.");

	fprintf(pipe, //"set terminal postscript eps color enhanced\n"
				  //"set output '%s'\n"
				  "unset key\n"
				  "unset tics\n"
				  "set size ratio -1\n"
				  "set xrange [%f:%f]\n"
				  "set yrange [%f:%f]\n"
				  "set cbrange [%f:%f]\n"
				  "set cbtics\n"
				  //"set pm3d map\n"
				  "set palette rgbformulae 33,13,10\n"
				  //"unset colorbox\n"
				  //"plot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n"
				  //"exit\n", epsfile, 
				  ,xmin-5, xmax+5, ymin-5, ymax+5, 
				  -sd[WS], sd[WS]);//, tmpfile,
				  //50.0/sqrt((float)maxptcl));
	
	const float pointsize=50.0/sqrt((float)maxptcl);
	const int sleeptime=1000000/framerate;

	//for (t=0; t<3;t++)
	for (t=frameindex[0]; t<=frameindex[1];t+=frameindex[2])
	//for (t=0; t<total_t; t++)
	{
		/*==========================================================*/
		// draw config by Gnuplot
		infile=fopen(tmpfile, "w");
		FILE_NULL(infile, tmpfile);
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			spointer=stn_ptr+((pointer-ptid_ptr)/size[0])*stn_size[0];
			if ( *spointer > 1.0e5 ) // no neighbor particle
				continue;
			if (WS==0)
				t2=*spointer+*(spointer+3);  //bulk strain
			else if (WS==1)
				t2=(*(spointer+1) + ( *(spointer+2) ))/2.0; // shear strain
			else
				t2=(*(spointer+2) - ( *(spointer+1) ))/2.0; // local rotation

			if ( t2 > sd[WS] )
				t2 = sd[WS];
			else if ( t2 < -sd[WS] )
				t2 = -sd[WS];
			fprintf(infile, "%f %f %g\n", *pointer, *(pointer+1),
							t2);//*(spointer) );
		}
		fclose(infile);

		printf("Frame: %d\n", t);
		fprintf(pipe, "plot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n"
					  ,tmpfile, pointsize);
		fflush(pipe);  // to put them into pipe and leave buffer
		usleep(sleeptime);

	}
	fprintf(pipe, "exit\n");
	fclose(pipe);
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);

	free(index);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	free(sd);
	free(frameindex);
	delete [] size;
	delete [] stn_size;
}
