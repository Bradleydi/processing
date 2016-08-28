#include "statistics.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void pt_info(const char *file, int dim, bool show_histogram)
{
	int i,j;
	char *gdffile=getfilename(file, ".gdf");
	colloid_base pt;
	readgdf(pt, gdffile);
	float *pt_ptr=pt.get_array_pointer();
	int * size=pt.get_size();
	int ncol=size[0];
	if (dim==0)
		dim=ncol-1;
	else if (dim>ncol-1)
		pERROR("given dim is larger than the given data.");
	else if (dim <0)
		pERROR("given dim given is negative.");
	
	int ti=size[0]-1;

	float *pointer=pt_ptr+ti;
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
	pointer=pt_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		*pointer=(float)((int)(*pointer)-tmin);
		pointer+=size[0];
	}

	int * ptcl_frame=Calloc(int, total_t); POINTER_NULL(ptcl_frame);

	float * min=Malloc(float, dim); POINTER_NULL(min);
	float * max=Malloc(float, dim); POINTER_NULL(max);
	pointer=pt_ptr;
	for (j=0; j<dim; j++)
	{
		min[j]=*pointer;
		max[j]=*pointer;
		++pointer;
	}
	pointer+=size[0]-dim-1;
	++ptcl_frame[(int)(*pointer)];
	for (i=1; i<size[1]; i++)
	{
		++pointer;
		for (j=0; j<dim; j++)
		{
			if ( min[j] > *pointer )
				min[j] = *pointer;
			else if ( max[j] < *pointer )
				max[j] = *pointer;
			++pointer;
		}
		pt_ptr+=ncol-dim-1;
		++ptcl_frame[(int)(*pointer)];
	}
	
	int minp=ptcl_frame[0], maxp=minp;
	for (i=1; i<total_t; i++)
	{
		if (minp > ptcl_frame[i] )
			minp=ptcl_frame[i];
		else if (maxp < ptcl_frame[i])
			maxp=ptcl_frame[i];
	}
	
	int * ht=Calloc(int, maxp-minp+1); POINTER_NULL(ht);

	for (i=0; i<total_t; i++)
		++ht[ptcl_frame[i]-minp];
	
	printf("# Dimension:        %d\n", dim);
	printf("# Total frame:      %d\n", total_t);
	if (tmin!=0)
		printf("# frame index:       [ %d : %d ]\n", tmin, (int)t2);
	printf("\n#              Min        Max\n");
	
	for (int i=0; i<dim; i++)
		printf("#    dim %d:    %f          %f\n", i, *(min+i), *(max+i));
	printf("\n");

	double ave=(double)(size[1])/total_t;
	int iave=(int)ave;
	unsigned long sd=0;
	for (i=minp; i<maxp+1; i++)
		sd+=((unsigned long)ht[i-minp])*(unsigned long)((i-iave)*(i-iave));
	printf("# Average partilces in a frame:        %f\n", ave);
	printf("# Particle fluctuation:                %f\n"
		   "# (standard derivation of particles in a frame)\n", 
		   sqrt((double)(sd)/total_t+2*ave*iave-(double)(iave*iave)-ave*ave));


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
			if (ht[i]!=0)
				fprintf(fp, "%d %d\n", i+minp, ht[i]);
		}
		Fclose(fp, tmpfile);
		fprintf(pipe, "set term x11 1\n"
					  "unset key\n"
					  "set xlabel 'particle number in a frame'\n"
					  "set ylabel 'histogram'\n"
					  "plot '%s' w points pt 7 ps 1.5\n"
					  , tmpfile);
		fprintf(pipe, "exit\n");
		fclose(pipe);
		if ( remove(tmpfile) != 0 )
			fprintf(stderr, 
					"# Error: tmpfile '%s' cannot be removed.\n", tmpfile);
	}

	delete [] size;
	free(ptcl_frame);
	free(min);
	free(max);
	free(ht);
}
