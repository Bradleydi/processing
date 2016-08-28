#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

// used only for dim=2
void traj_length(const char *file)
{
	char *gdffile = getfilename (file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);

	float *ptid_ptr=ptid.get_array_pointer();
	int *size=ptid.get_size();
	const int ti=size[0]-2;
	const int ii=size[0]-1;
	int i, p;

	const int dim=2;
	if ( dim+2 > size[0] )
		pERROR("data should be tracked and dim = 2.\n");

	const int total_p=(int)(ptid_ptr[size[2]-1])+1;

	// appeared frame number for particle
	int *frame_ptcl=Calloc(int, total_p);
	POINTER_NULL(frame_ptcl);
	/*
	float *pointer=ptid_ptr+(size[0]-1); // id
	int i;
	for (i=0; i<size[1]; i++)
	{
		++frame_ptcl[(int)(*pointer)];
		pointer+=size[0];
	}*/

	/* this is wrong. If no particle exists in all frame, i.e.
	 *  for all i, frame_ptcl[i] < total_t
	 *  so that we can not get the right total_t
	int total_t=frame_ptcl[0];
	for (i=1; i<total_p; i++)
	{
		if ( total_t < frame_ptcl[i] )
			total_t = frame_ptcl[i];
	}
	*/
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

		pointer+=ti;
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
		++pointer;  // id
		++frame_ptcl[(int)(*(pointer++))]; // pointer -> next x
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin = (int)t1;
	printf("# total frame               = %d\n"
		   "# total particle            = %d\n", total_t, total_p);
	if (tmin!=0) printf("# \tframe index [ %d : %d ]\n", tmin, (int)t2);

	// get mean position
	float *mp = Calloc(float, 2*total_p); POINTER_NULL(mp);
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		p=(int)(*(pointer+ii)); // id
		mp[2*p]   += *pointer;    // x
		mp[2*p+1] += *(pointer+1); // y
		pointer += size[0];
	}

	for (i=0; i<total_p; i++)
	{
		*(mp++) /= frame_ptcl[i];       // x
		*(mp++) /= frame_ptcl[i];   // y
	}
	mp -= 2*total_p;
	
	char *tmpfile = tmpnam (NULL);
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);

	fprintf(infile, "# spatial distribution of"
			"trajectory lengthes of \"%s\"\n", gdffile);

	for (i=0; i<total_p; i++)
		fprintf(infile, "%.6f %.6f %.6f\n", mp[2*i], mp[2*i+1],
				        (float)(frame_ptcl[i])/total_t);
	Fclose(infile, tmpfile);
	free(gdffile);
	free(frame_ptcl);
	free(mp);
	delete [] size;

	// plot by Gnuplot
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	float Xmargin=(xmax-xmin)*0.05, Ymargin=(ymax-ymin)*0.05;
	const float pointsize=30.0/sqrt((float)size[1]/(float)total_t);
	fprintf(pipe, "reset\n"
			      "unset key\n"
				  "set pm3d map\n"
				  "set palette rgbformulae 33,13,10\n"
				  "set size ratio -1\n"
				  "set tmargin at screen 0.95\n"
				  "set xrange [%f:%f]\n"
				  "set yrange [%f:%f]\n"
				  "splot '%s' using 1:2:3 with points palette pt 7 ps %.2f\n",
				  xmin-Xmargin, xmax+Xmargin, 
				  ymin-Ymargin, ymax+Ymargin,
				  tmpfile, pointsize);
	fclose(pipe);
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);
}
