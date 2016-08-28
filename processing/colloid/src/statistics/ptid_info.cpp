// version 1.0

#include "statistics.h"
#include "io.h"
#include "miscellaneous.h"

#include <cmath>
#include <cstdlib>


void ptid_info(const char *file, int dim, bool show_histogram)
{
	int i,j;
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	float *ptid_ptr=ptid.get_array_pointer();
	int * size=ptid.get_size();
	int ncol=size[0];
	if (dim==0)
		dim=ncol-2;
	else if (dim>ncol-2)
		pERROR("given dim is larger than the given data.");
	else if (dim <0)
		pERROR("given dim given is negative.");

	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	int * frame_ptcl=Calloc(int, total_p); POINTER_NULL(frame_ptcl);

	float *pointer=ptid_ptr+ncol-1; // id
	for (i=0; i<size[1]; i++)
	{
		++frame_ptcl[(int)(*(pointer))];
		pointer+=ncol;
	}

	int ti=size[0]-2;
	pointer=ptid_ptr+ti;
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
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		*pointer=(float)((int)(*pointer)-tmin);
		pointer+=size[0];
	}

	int * ptcl_frame=Calloc(int, total_t); POINTER_NULL(ptcl_frame);
	
	float * min=Malloc(float, dim); POINTER_NULL(min);
	float * max=Malloc(float, dim); POINTER_NULL(max);
	pointer=ptid_ptr;
	for (j=0; j<dim; j++)
	{
		min[j]=*pointer;
		max[j]=*pointer;
		++pointer;
	}
	pointer+=ncol-dim-2;  // t
	++ptcl_frame[(int)(*pointer)];
	pointer+=2;
	for (i=1; i<size[1]; i++)
	{
		for (j=0; j<dim; j++)
		{
			if ( min[j] > *pointer)
				min[j]=*pointer;
			else if (max[j] < *pointer)
				max[j]=*pointer;
			++pointer;
		}
		pointer+=ncol-dim-2;
		++ptcl_frame[(int)(*pointer)]; // t
		pointer+=2;
	}
	
	int minp=ptcl_frame[0], maxp=minp;
	for (i=1; i<total_t; i++)
	{
		if (minp > ptcl_frame[i] )
			minp=ptcl_frame[i];
		else if (maxp < ptcl_frame[i])
			maxp=ptcl_frame[i];
	}
	
	// histogram of ptcl_frame
	int * hp_t=Calloc(int, maxp-minp+1); POINTER_NULL(hp_t);
	for (i=0; i<total_t; i++)
		++hp_t[ptcl_frame[i]-minp];


	printf("# Dimension:        %d\n", dim);
	printf("# Total particle:   %d\n", total_p);
	printf("# Total frame:      %d\n", total_t);
	if (tmin!=0)
		printf("# frame index:       [ %d : %d ]\n", tmin, (int)t2);
	if (total_t*total_p>=2*size[1])
		printf("# Warning: too many particles, should be tracked better.\n");
	printf("\n#              Min        Max\n");
	
	for (int i=0; i<dim; i++)
		printf("#    dim %d:    %f          %f\n", i, *(min+i), *(max+i));
	printf("\n");

	double ave=(double)(size[1])/total_t;
	int iave=(int)ave;
	unsigned long sd=0;
	for (i=minp; i<maxp+1; i++)
		sd+=((unsigned long)hp_t[i-minp])*(unsigned long)((i-iave)*(i-iave));
	printf("# Average partilces in a frame:        %f\n", ave);
	printf("# Particle fluctuation:                %f\n"
		   "# (standard derivation of particles in a frame)\n", 
		   sqrt((double)(sd)/total_t+2.0*ave*iave-(double)(iave*iave)-ave*ave));


	// histogram of frame_ptcl
	int * ht_p=Calloc(int, total_t+1); POINTER_NULL(ht_p);
	for (i=0; i<total_p; i++)
		++ht_p[frame_ptcl[i]];

	
	ave=(double)(size[1])/total_p;
	iave=(int)ave;
	sd=0;
	for (i=1; i<=total_t; i++)
		sd+=((unsigned long)ht_p[i])*(unsigned long)((i-iave)*(i-iave));
	printf("\n# Average trajectory length:           %f\n", ave);
	printf("# Trajectory length fluctuation:       %f\n"
		   "# (standard derivation of tracjetory length)\n",
		   sqrt((double)(sd)/total_p+2.0*ave*iave-(double)(iave*iave)-ave*ave));
	
	if (show_histogram) {
		FILE *pipe=popen("gnuplot -persist", "w");
		if (pipe==NULL) {
			fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
			exit (1);
		}
		//plot ptcl_frame
		char *tmpfile= tmpnam (NULL);
		FILE *fp = fopen(tmpfile, "w"); FILE_NULL(fp, tmpfile);
		for (i=0; i<total_t; i++)
			fprintf(fp, "%d %d\n", i+tmin, ptcl_frame[i]);
		Fclose(fp, tmpfile);
		fprintf(pipe, "set term x11 0\n"
					  "unset key\n"
					  "set xlabel 'frame index'\n"
					  "set ylabel 'particle number in a frame'\n"
					  "plot '%s' w points pt 7 ps 1.5\n"
					  , tmpfile);
		fflush(pipe);
		sleep(1);

		//printf("\n# Particle number in a frame histogram:\n");
		fp = fopen(tmpfile, "w"); FILE_NULL(fp, tmpfile);
		for (i=0; i<maxp-minp+1; i++)
		{
			if (hp_t[i]!=0)
				fprintf(fp, "%d %d\n", i+minp, hp_t[i]);
		}
		Fclose(fp, tmpfile);
		fprintf(pipe, "set term x11 1\n"
					  "unset key\n"
					  "set xlabel 'particle number in a frame'\n"
					  "set ylabel 'histogram'\n"
					  "plot '%s' w points pt 7 ps 1.5\n"
					  , tmpfile);
		fflush(pipe);
		sleep(1);

		// plot frame_ptcl
		fp = fopen(tmpfile, "w"); FILE_NULL(fp, tmpfile);
		for (i=0; i<total_p; i++)
			fprintf(fp, "%d %d\n", i, frame_ptcl[i]);
		Fclose(fp, tmpfile);
		fprintf(pipe, "set term x11 2\n"
					  "unset key\n"
					  "set xlabel 'particle index'\n"
					  "set ylabel 'particle trajectory length'\n"
					  "plot '%s' w points pt 7 ps 1.5\n"
					  , tmpfile);
		fflush(pipe);
		sleep(1);

		//printf("\n# Trajectory length histogram:\n");
		fp = fopen(tmpfile, "w"); FILE_NULL(fp, tmpfile);
		for (i=0; i<=total_t; i++)
		{
			if (ht_p[i]!=0)
				fprintf(fp, "%d %d\n",  i, ht_p[i]);
		}
		Fclose(fp, tmpfile);
		fprintf(pipe, "set term x11 3\n"
					  "unset key\n"
					  "set xlabel 'particle index'\n"
					  "set ylabel 'particle trajectory length'\n"
					  "plot '%s' w points pt 7 ps 1.5\n"
					  , tmpfile);
		fprintf(pipe, "exit\n");
		fclose(pipe);
		if ( remove(tmpfile) != 0 )
			fprintf(stderr, 
					"# Error: tmpfile '%s' cannot be removed.\n", tmpfile);
	}
	
	delete [] size;
	free(frame_ptcl);
	free(ptcl_frame);
	free(min);
	free(max);
	free(hp_t);
	free(ht_p);
}
