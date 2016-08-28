#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

/* Use steady state detection (SSD) to detect irreversible rearrangment (IR)
 * Other methods can also be implemented, but here I just want to make it
 * easier
 */

void IR_SSD(const char *file, int dim, int window)
{
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);

	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();

	if (dim<=0) pERROR("Dimension should be larger than 0\n");
	if (size[0]<dim+2) pERROR("Not a tracked data\n");
	const int ti=size[0]-2; //time index in the array
	const int ii=ti+1; // particle index in the array

	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	/* get total_t */
	int i, p, d, NSS;
	float *pointer=ptid_ptr+ti;
	float t1=*pointer, t2=t1;
	for (i=1; i<size[1]; i++)
	{
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
		pointer += size[0];
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin = (int)t1;
	printf("# total frame               = %d\n"
			"# total particle            = %d\n", total_t, total_p);
	if (tmin!=0) printf("# \tframe index [ %d : %d ]\n", tmin, (int)t2);

	/* reindex the data by ptcl, not time as usual */
	int *frame_ptcl=Calloc(int, total_p); POINTER_NULL(frame_ptcl);
	pointer=ptid_ptr+ii;
	for (i=0; i<size[1]; i++)
	{
		++frame_ptcl[(int)(*pointer)];
		pointer += size[0];
	}

	int *cum_frame_ptcl=Malloc(int, (total_p+1)); POINTER_NULL(cum_frame_ptcl);
	cum_frame_ptcl[0]=0;
	cum_frame_ptcl[1]=frame_ptcl[0];
	int maxframe=frame_ptcl[0];
	for (p=1; p<total_p; p++)
	{
		cum_frame_ptcl[p+1]=cum_frame_ptcl[p]+frame_ptcl[p];
		if ( maxframe < frame_ptcl[p] ) maxframe = frame_ptcl[p];
	}
	if (maxframe*total_p > 2*size[1])
		pERROR("too many broken trajectories!");
	int *index_now=frame_ptcl;
	for (p=0; p<total_p; p++) index_now[p]=0;
	float **index=(float **)malloc(size[1]*sizeof(float *));
	POINTER_NULL(index);
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		p=(int)(*(pointer+ii));
		index[cum_frame_ptcl[p]+(index_now[p]++)]=pointer;
		pointer+=size[0];
	}

		//printf("%d %d %d\n", frame_ptcl[0], frame_ptcl[3], frame_ptcl[100]);
	double *position=Malloc(double, maxframe); POINTER_NULL(position);
	int *SS = Malloc(int, maxframe/2); POINTER_NULL(SS);
	char *workspace = Malloc(char, maxframe); POINTER_NULL(workspace);
	if ( window == 0 )
	{
		window =  maxframe/10;
		fprintf(stderr, "# set window width to %d\n", window);
	}
	for (p=0; p<total_p; p++)
	{
		for (d=0; d<dim; d++)
		{
			for (i=cum_frame_ptcl[p]; i< cum_frame_ptcl[p+1]; i++)
				*(position++)=(double)(*(index[i]+d));
			position -= frame_ptcl[p];
			//printf("calling SSD_batch for %d...\n", p);
			SSD_batch(frame_ptcl[p], position, window, &NSS, SS, workspace);
			//printf("end\n");
			if (NSS!=1)
			{
				printf("# particle %d has %d local equilibrium states\n"
					   "#     SSD_id    SSD_interval\n",
						p, NSS);
				for (i=0; i<NSS; i++)
					printf("#     %d       [ %d : %d ]\n", i, 
							(int)*(index[cum_frame_ptcl[p]+SS[2*i]]+ti),
							(int)*(index[cum_frame_ptcl[p]+SS[2*i+1]]+ti));
			}
		}
	}
	 p= 1921 ;
	 d=0;
	for (i=cum_frame_ptcl[p]; i< cum_frame_ptcl[p+1]; i++)
		printf("%d %f\n", d++, *index[i] );
	free(position);
	free(SS);
	free(workspace);
	free(frame_ptcl);
	free(cum_frame_ptcl);
	free(index);
	delete [] size;
}
