#include "data_preprocess.h"
#include "colloid_base.h"
#include "miscellaneous.h"
#include "io.h"
#include "statistics.h"

#include <cstdio>
#include <cstdlib>

// used for  nonlinear drift ( This seems to be the case in usual)
// data should be tracked
void drift_llr(const char *file, int dim, int smooth)
{
	printf("# Removing drift by local linear regression...\n");
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);
	
	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();
	const int ti = size[0]-2;
	const int ii = size[0]-1;

	/* test tracked data */
	if ( dim <= 0 )
	{
		dim = ti;
		printf("# Set DIM to be %d\n", dim);
	}
	else if ( dim+2 > size[0] )
		pERROR("not a tracked data!");
	
	/* get total_p  & total_t */
	int i, j, t;
	float *pointer=ptid_ptr+ti;
	float tt=*pointer, pp=pointer[1];
	float t1=tt, t2=tt, p1=pp, p2=pp;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( t1 < *pointer )
			t1 = *pointer;	// tmax
		else if ( t2 > *pointer)
			t2 = *pointer;  // tmin
		
		pp=pointer[1];
		if ( p1 < pp )
			p1 = pp;
		else if ( p2 > pp )
			p2 = pp;
	}
	
	if ( (int)p2 != 0 )
		pERROR("wrong particle index");
	
	const int total_p=(int)(p1)+1;
	const int tmin = (int)t2;
	const int total_t=(int)t1 - tmin + 1;
	printf("# total_p\t=\t%d\n# total_t\t=\t%d\n", total_p, total_t);
	if (tmin!=0)
		printf("# trange [%d:%d]\n", tmin, (int)t1);
	
	// don't reindex time to avoid side effect
	int *ptcl_frame = Calloc(int, total_t); POINTER_NULL(ptcl_frame);
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer) - tmin];
		pointer += size[0];
	}

	int *cum_ptcl_frame = Malloc(int, total_t+1); POINTER_NULL(cum_ptcl_frame);
	cum_ptcl_frame[0]=0;
	for (i=0; i<total_t; i++)
		cum_ptcl_frame[i+1] = cum_ptcl_frame[i]+ptcl_frame[i];
	//printf("%d\n", cum_ptcl_frame[total_t-1]);

	// reindex the data
	float **index=Malloc(float *, size[1]); POINTER_NULL(index);
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		index[cum_ptcl_frame[(int)(pointer[ti])-tmin]++]=pointer;
		pointer += size[0];
	}
	// reconstruct cum_ptcl_frame
	cum_ptcl_frame[0]=0;
	for (i=0; i<total_t; i++)
		cum_ptcl_frame[i+1] = cum_ptcl_frame[i]+ptcl_frame[i];

	free(ptcl_frame);
	
	/*===================================================================
	 * find out drift between successive frames
	 *===================================================================*/
	/* If there is a frame having no particles, then we throw away the data
	 * for the successive drift.
	 * However, this will affect our calculation for the intergrated drift.
	 * The best way seems to intertroplate some data for the missing drift.
	 * Nevertheless, we can try to find the drift between two frames with a gap
	 * and in this way, we need not to intertroplate or guess.
	 */
	// i is for frame 1, and j is for frame 2
	
	// exist is used to store the position in cum_ptcl_frame for
	// further usage.
	// if exist is 0, then it's a particle is not contained in that frames.
	// else the value exist[i]-1 is the poision in cum_ptcl_frame
	// Actually, this method is quite dangerous, since if one forgets to
	// reset good[i] to 0, then exist[i] will still in the array range,
	// and 0 will be used, which is index[0]
	int *exist1=Malloc(int, total_p); POINTER_NULL(exist1);
	int *exist2=Malloc(int, total_p); POINTER_NULL(exist2);
	// good particles are those in both frame
	// A good particle plot should be left to the track program,
	// just to make sure the track is good
	unsigned char *good=Malloc(unsigned char, total_p); POINTER_NULL(good);
	int tnext, Ndrift=0, Ngood;
	// drift has dim+1 columns, 1 for time
	// [t0 t1 t2 ... xdrift1 xdrift2 ... ydrift1 ...]
	double *drift=Calloc(double, total_t*(dim+1)); POINTER_NULL(drift);
	float *pointer2;
	double *pdrift=drift;
	*(pdrift++)=(double)tmin;
	//for ( t=0; t<total_t; t++)
	printf("# Getting drift ...\n");
	t=0;
	while (t<total_t-1)
	{
		if ( cum_ptcl_frame[t+1] == cum_ptcl_frame[t] ) 
		{
			// actually, due to the coding here, this part will never happen.
			// the t=tnext will always bring the data to a frame containing
			// particles
			printf("# Warning: frame %d contains no particles\n"
				   "#          This frame is skipped for calculating drift\n",
				   t+tmin);
			continue;
		}
		tnext = t+1;
		if ( cum_ptcl_frame[tnext+1] == cum_ptcl_frame[tnext] ) 
		{
			printf("# Warning: frame %d contains no particles\n"
				   "#          This frame is skipped for calculating drift\n",
				   tnext+tmin);
			continue;
		}

		for (i=0; i<total_p; i++) {	exist1[i]=0; exist2[i]=0; good[i]=0; }
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i]+ii;
			exist1[(int)(*pointer)]=i+1;
		}
		for (j=cum_ptcl_frame[tnext]; j<cum_ptcl_frame[tnext+1]; j++)
		{
			pointer2=index[j]+ii;
			exist2[(int)(*pointer2)]=j+1;
		}
		Ngood=0;
		for (i=0; i<total_p; i++) 
			if ( exist1[i] && exist2[i]) { good[i]=1; ++Ngood; }
		if ( Ngood < total_p/2 && Ngood < 300 )
			printf("# Warning: particles containing in two successive frames\n"
				   "#          is too few ( Ngood = %d frame_id = %d, %d)\n", 
				   Ngood, t+tmin, tnext+tmin);
		// calculate the mean x, y (and z) for the good particles
		// Assuming that the particle index is ordered
	
		for (i=0; i<total_p; i++)
		{
			if (good[i])
			{
				pointer=index[exist1[i]-1];
				pointer2=index[exist2[i]-1];
				for (j=0; j<dim; j++) 
					pdrift[total_t*(j+1)] += (double)(pointer2[j] - pointer[j]);
			}
		}
		for (j=0; j<dim; j++) 
			pdrift[total_t*(j+1)] /= Ngood;
		*(pdrift++)=(double)(tnext+tmin);
		++Ndrift;

		t=tnext; // directly go to the next time, avoiding the frames contains
				//no particles
	}
	free(exist1); free(exist2); free(good);
	free(cum_ptcl_frame); free(index);
	ptid.free_memory();
	delete [] size;
	printf("OK\n");
	/*================================================================
	 * get drift end ( not the integrated drift)
	 *==================================================================*/

	/*============================================================*/
	/* get integrated drift */
	pdrift=drift+total_t;  // xdrift
	for (j=0; j<dim; j++) 
	{
		for (t=1; t<Ndrift; t++)
			pdrift[t] += pdrift[t-1];
		pdrift += total_t;
	}
	/*================================================================
	 * get drift end ( integrated drift)
	 *==================================================================*/

	/*============================================================*/
	/* local linear regress the drift */
	/* time is just drift[0:Ndrift-1] */
	printf("# smoothing the drift by local linear regression...\n");
	double *reg=Malloc(double, total_t*dim); POINTER_NULL(reg);
	kernel_t KERNEL = EPANECHNIKOV; // next version may be to change kernel
	for (j=0; j<dim; j++) 
		locallinear(Ndrift, drift, drift+(j+1)*total_t, reg+j*total_t, smooth, KERNEL);
	
	/*================================================================
	 * plot drift & regression result
	 *==================================================================*/
	char *tmpfile = tmpnam (NULL);
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);
	pdrift=drift;
	for (t=0; t<Ndrift; t++)
	{
		fprintf(infile, "%f ", *pdrift);
		for (j=0; j<dim; j++) 
			fprintf(infile, "%g ", pdrift[(j+1)*total_t]);
		for (j=0; j<dim; j++) 
			fprintf(infile, "%g ", reg[j*total_t+t]);
		fprintf(infile, "\n");
		++pdrift;
	}
	Fclose(infile, tmpfile);
	free(drift);
	free(reg);


	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe for gnuplot!\n");
		exit (1);
	}
	
	fprintf(pipe,  "reset\n"
				   "set title \"drift of %s\\n"
				   "total particle = %d\"\n"
				   "set xlabel 't'\n"
				   "set ylabel 'drift'\n"
				   "set xrange [%f:%f]\n"
				   , file, total_p, t2, t1);

	if (dim==1)
		fprintf(pipe, "plot \"%s\" using 1:2 title \"x\" w p lc rgb \"red\","
			"\"\" using 1:3 title \"x-drift\" w l lw 3 lc rgb \"red\"\n"
					   ,tmpfile);
	else if (dim ==2 )
		fprintf(pipe, "plot \"%s\" using 1:2 title \"x\" w p lc rgb \"red\","
			"\"\" using 1:4 title \"x-drift\" w l lw 3 lc rgb \"red\","
		    "\"\" using 1:3 title \"y\" w p lc rgb \"blue\","
			"\"\" using 1:5 title \"y-drift\" w l lw 3 lc rgb \"blue\"\n"
					   ,tmpfile);
	else if (dim ==3 )
		fprintf(pipe, "plot \"%s\" using 1:2 title \"x\" w p lc rgb \"red\","
			"\"\" using 1:5 title \"x-drift\" w l lw 3 lc rgb \"red\","
		    "\"\" using 1:3 title \"y\" w p lc rgb \"blue\","
			"\"\" using 1:6 title \"y-drift\" w l lw 3 lc rgb \"blue\"\n"
		    "\"\" using 1:4 title \"z\" w p lc rgb \"black\","
			"\"\" using 1:7 title \"z-drift\" w l lw 3 lc rgb \"black\","
					   ,tmpfile);
	else
		fprintf(stderr, "# Error: plotting has not been implemented yet\n");

	fclose(pipe);

	if( remove( tmpfile ) != 0 )
	{
		fprintf(stderr, "# Error: deleting file \"%s\" failed!", tmpfile);
		exit (1);
	}
}
