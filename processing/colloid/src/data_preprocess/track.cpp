#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

void track(const char *file, const float maxdisp, const int memory)
{
	const float maxdisp2=maxdisp*maxdisp;

	colloid_base pt;
	readgdf(pt, file);

	int *size=pt.get_size();
	float *pt_ptr=pt.get_array_pointer();

	const int dim=2;
	if (size[0] < dim +1 )
		pERROR("not a 2D data.");
	const int ti=size[0]-1;
	const int total_t = (int)pt_ptr[size[2]-1]+1;

	int *ptcl_frame = Calloc(int, total_t); POINTER_NULL(ptcl_frame);
	int i, j, k;
	float *pointer=pt_ptr+ti; // t
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	// first pointer of each time
	int maxptcl_frame=ptcl_frame[0];
	float **begin = Malloc(float *, total_t); POINTER_NULL(begin);
	begin[0]=pt_ptr;
	for (i=1; i<total_t; i++)
	{
		begin[i]=begin[i-1]+ptcl_frame[i-1];
		if ( maxptcl_frame < ptcl_frame[i] )
			maxptcl_frame = ptcl_frame[i];
	}


	// save edges
	int total_edge;
	int total_cluster;
	const int maxedge=5*maxptcl_frame;
	printf("# maxedge= %d\n", maxedge);
	int *edge = Malloc(int, 2*maxedge); POINTER_NULL(edge);
	int *cluster = Malloc(int, maxedge); POINTER_NULL(cluster);
	int *Ncluster = Malloc(int, maxedge); POINTER_NULL(Ncluster);
	char *calculated = Malloc(char, maxedge); POINTER_NULL(calculated);
	int *id = Malloc(int, size[1]); POINTER_NULL(id);
	for (i=0; i<ptcl_frame[0]; i++)
		id[i]=i;

	float *pointerj, *pointerk;
	float dy, d;
	int t;
	char newcluster;
	for (t=0; t<1;t++)//total_t; i++)
	{
		total_edge=0;
		pointerj=begin[t];
		for (j=0; j<ptcl_frame[t]; j++)
		{
			pointerk=begin[t+1];
			for (k=0; k<ptcl_frame[t+1]; k++)
			{
				d = *pointerj - *pointerk;
				d*=d;
				if ( d < maxdisp2 )
				{
					dy= *(pointerj+1) - *(pointerk+1);
					d+=dy*dy;
					if ( d < maxdisp2 )
					{
						*(edge++)=(pointerj-pt_ptr)/size[0];
						*(edge++)=(pointerk-pt_ptr)/size[0];
						++total_edge;
						if ( total_edge > maxedge )
							pERROR("too big maxdisp!");
					}
				}
				pointerk+=size[0];
			}
			pointerj+=size[0];
		}

		// check isolated edge
		edge-=2*total_edge;
		cluster[0]=0;
		total_cluster=0;
		Ncluster[total_cluster]=1;
		for (i=0; i<total_edge; i++)
		{
			newcluster=1;
			for (j=0; j<i; j++)
			{
				if ( edge[2*i] == edge[2*j] ||
					 edge[2*i+1] == edge[2*j+1] )
				{
					cluster[i]=cluster[j];
					++Ncluster[cluster[i]];
					newcluster=0;
					break;
				}
			}
			if (newcluster)
			{
				Ncluster[total_cluster]=1;
				cluster[i]=(total_cluster++);
			}
		}
		printf("# total cluster : %d\n", total_cluster);
		printf("# total edge : %d\n", total_edge);
		printf("# total particle : %d  %d\n", ptcl_frame[t], ptcl_frame[t+1]);
		for (i=0; i<total_edge; i++)
			calculated[i]=0;
		for (i=0; i<total_edge; i++)
		{
			if (calculated[i])
				continue;

			calculated[i]=1;
			if ( Ncluster[cluster[i]] == 1 ) // isolated
				id[edge[2*i+1]]=id[edge[2*i]];
			else // nontrivial cases
			{
				printf("# total edge in cluster[%d] : %d\n", 
						cluster[i], Ncluster[cluster[i]]);
				//if ( Ncluster[cluster[i]] >= 7 )
				//	pERROR("difficult combinatorics encountered.");
				for(j=i+1; j<total_edge; j++)
					if (Ncluster[cluster[j]]==Ncluster[cluster[i]])
						calculated[j]=1;
			}
		}
	}
	

	free(begin);
	free(edge);
	free(cluster);
	free(id);

	free(ptcl_frame);
	delete [] size;
}
