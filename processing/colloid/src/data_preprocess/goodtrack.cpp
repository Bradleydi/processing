#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

/* Definition:
 * Good particles (well tracked particles)
 * 		In given time window [t0, t0+L], the number of frames in which 
 * 		a particle exists exceeds a given threshold, we call this particle
 * 		as a good particle.
 *
 * Algorithm:
 * Just by counting the existing frame number for a particle given a time
 * window.
 * 1. get initial trajectory length
 * 2. check old start (s1) and new end (e2)
 * 		if (s1 exist) new length = old length -1
 * 		else		  new length = old length
 * 		if (e2 exist) new length = old length +1
 * 		else		  new length = old length
 * 3. if (length >= threshold_length) good particle
 * 4. account number of good particles, store
 * 		(time window width, number of good particles, start frame)
 * 	  after compared the number of good particle with old one
 * 5. shift the window and repeat 2-5.
 * 6. plot (time window width, max number of good particles, start frame)
 */

/* This program account the number of well tacked particles.
 * Since it only account trajectory length, so position data is useless,
 * in this case, we simply extract the time and particle index to save
 * memory(, and time).
 */

void goodtrack(const char *file, double threshold)
{
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);

	int i, t;
	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();

	//const int dim=2;
	if (size[0]<4)
		pERROR("dimension is not 2, or not a tracked data!");

	// first size[1] part is saved as time index, and the 2nd harf is ptcl index
	int *time_ptcl=Malloc(int, size[1]*2); POINTER_NULL(time_ptcl);
	
	int *pointer = time_ptcl;
	ptid_ptr -= 2;
	for (i=0; i<size[1]; i++)
	{
		ptid_ptr += size[0];
		*pointer = (int)(*ptid_ptr);
		*(pointer+size[1]) = (int)(*(ptid_ptr+1));
		++pointer;
	}

	const int N=size[1];
	// free memory
	ptid.free_memory();
	delete [] size;

	pointer = time_ptcl;
	int t1 = *pointer, t2 = t1;
	for (i=1; i<N; i++)
	{
		++pointer;
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
	}
	const int total_t = (t2-t1) + 1;
	const int tmin = t1;
	const int total_p = time_ptcl[N*2-1]+1;
	
	// reindex time
	pointer = time_ptcl;
	for (i=0; i<N; i++)
		*(pointer++) -= tmin;
	//============================================================
	// index search for time & ptcl.

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(ptcl_frame);
	pointer = time_ptcl;
	for (i=0; i<N; i++)
		++ptcl_frame[*(pointer++)];

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

	if (maxptcl*total_t > 10*N)
		pERROR("too many broken trajectories (10x more than actual number)!");

	int *index_now=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(index_now);

	// store the particle index, not array index
	int *index=(int *)malloc(N*sizeof(int));
	POINTER_NULL(index);
	pointer = time_ptcl;
	for (i=0; i<N; i++)
	{
		t=*pointer;
		index[cum_ptcl_frame[t]+(index_now[t]++)]= *(pointer+N);
		++pointer;
	}
	free(index_now);
	
	//==================================================================
	// find good track
	//==================================================================
	//const float threshold = 0.9; // for good tracked particles
	if (threshold < 0.0)
		threshold = 1.0;
	int window=1000;  // window width
	int threshold_length = (int)(threshold*window);
	int wstep=10; // window width increasment step
	int Nwindow=(total_t-window)/wstep;  // number of window width to calculate

	// initial window frame (start point) for finial result
	int *p0=Malloc(int, Nwindow); POINTER_NULL(p0);
	// max good particle number for finial result
	int *mgp=Malloc(int, Nwindow); POINTER_NULL(mgp);

	// for each window with different initial frame
	// good particles number
	//int *gp=Malloc(int, total_p); POINTER_NULL(gp);
	// length
	int *length=Malloc(int, total_p); POINTER_NULL(length);

	/*
	// storage of start points for each particle (used in loop)
	// if marked as -1, meaning no info for it.
	// So if the particle exists all the time, spp[p] should be
	// t-window;
	// However, this is not enough if we define good tracked particle
	// as particles exists through 90% of the period without broken traj.
	int *spp=Malloc(int, total_p); POINTER_NULL(spp);
	// longest trajectory length
	int *longest=Malloc(int, total_p); POINTER_NULL(longest);

	// the information about the trajectory connected to the window end
	int *spn=Malloc(int, total_p); POINTER_NULL(spn);
	int *length=Malloc(int, total_p); POINTER_NULL(length);
	*/

	int goodtraj;
	//int checkptcl;
	int wi=0, p;

	for (wi=0; wi<Nwindow; wi++)
	{
		window=wi*wstep+1000;
		threshold_length = (int)(threshold*window);

		// initialization
		for (p=0; p<total_p; p++)
			length[p]=0;
		//for (t=0; t<window; t++)
		//{
		//	for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		//	{
		//		p=index[i];
		//		++length[p];
		//	}
		//}
		for (i=0; i<cum_ptcl_frame[window+1]; i++)  // time from 0 to window-1
			++length[index[i]];
		goodtraj=0;
		for (p=0; p<total_p; p++)
			if ( length[p] >= threshold_length)
				++goodtraj;
		mgp[wi]=goodtraj;
		p0[wi]=tmin;
		
		for (t=1; t<total_t-window; t++)
		{
			// check old start
			for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++) 
				--length[index[i]];  // existed s1

			// check new end
			for (i=cum_ptcl_frame[t+window]; i<cum_ptcl_frame[t+window+1]; i++) 
				++length[index[i]];  // existed e2

			goodtraj=0;
			for (p=0; p<total_p; p++)
				if ( length[p] >= threshold_length)
					++goodtraj;
			if ( mgp[wi] < goodtraj )
			{
				mgp[wi]=goodtraj;
				p0[wi]=t;
			}
		}
	}

	free(time_ptcl);


	//==================================================================
	// plot
	//==================================================================

	char *tmpfile= tmpnam (NULL);
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);
	for (wi=0; wi<Nwindow; wi++)
		fprintf(infile, "%d %d %d\n", wi*wstep+1000, mgp[wi], p0[wi]);
	fclose(infile);
	free(mgp);
	free(p0);
	
	//==================================================================
	// open pipe to Gnuplot
	//==================================================================
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}
	fprintf(pipe, "unset key\n"
				  "set xrange [%d:%d]\n"
				  "set title \"total_p=%d\tthreshold=%f\"\n"
				  "set xlabel 'time window width'\n"
				  "set ylabel 'good particle number'\n"
				  "set y2label 'good track initial frame' rotate by -90\n"
				  "set y2tics\n"
				  "set ytics nomirror\n"
				  "plot '%s' axis x1y1, '%s' u 1:3 axis x1y2\n", 
				  1000, total_t, total_p, threshold, tmpfile, tmpfile);

	fclose(pipe);
	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);

	free(index);
	free(ptcl_frame);
	free(cum_ptcl_frame);
}
