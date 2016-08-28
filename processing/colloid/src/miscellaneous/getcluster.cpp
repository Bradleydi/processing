#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

/* to speed up the program, since it'll be called for each frame
 * which is a large amount of task, so I think it's a good idea to 
 * allocate memory first.
 * there are potential 2*nedge nodes in the graph, so that 
 * cluster should be allocated by a mount of
 * 	2*nnode<=4*nedge
 * memory.
 *
 * Array cluster will sort 2*Nnode elements, and first Nnode are
 * node indices, and the second half is the corresponding cluster indices.
 *
 * There are two methods to access cluster indices, both based on pointers
 * of nodes:
 * 1. save those pointers in an array indexed from 0 to (Nnode-1)
 * 2. use binary search for each node encountered.
 *
 * Of course, method 2 is memory-saved, While method 1 is time-saved.
 * Since in my program, Nedge ~ O(Nnode) ~ O(10^3), so we implement 
 * method 1 here.
 *
 * I also found that just storing the cluster indices, not the pointers of
 * nodes is better.
 */
void getcluster(int Nedge, int *edges, 
		int *pNnode, int *pNcluster, int *cluster)
{
	int i;
	int N=Nedge*2;
	for (i=0; i<N; i++)
		cluster[i]=edges[i];
	sort_int(N, cluster);
	int Nnode=1, MaxNode;
	for (i=1; i<N; i++)
		if ( cluster[i] != cluster[i-1] )
			cluster[Nnode++]=cluster[i];
	MaxNode=cluster[Nnode-1]+1;

	// working place to storage of cluster indices
	int *work;

	if ( Nnode == MaxNode )
		work=cluster+Nnode;
	else if ( Nnode < MaxNode )// there are some node indices not
		work=Malloc(int, MaxNode); // appeared in the edge list
	else
		pERROR("Something wrong on sorting or uniquing node list.");
	POINTER_NULL(work);

	for (i=0; i<MaxNode; i++)
		work[i]=-1;

	int j, Ncluster=0, p, q, large, small;
	for (i=0; i<Nedge; i++)
	{
		p=*(edges++);
		q=*(edges++);
		if ( work[p]==-1 && work[q]==-1 ) // not assigned a cluster index
		{
			work[p]=Ncluster; work[q]=Ncluster++;
		}
		else if ( work[p]==-1 && work[q]!=-1 ) // p is not assigned, q is
			work[p]=work[q];
		else if ( work[p]!=-1 && work[q]==-1 ) // q is not assigned, p is
			work[q]=work[p];
		else if ( work[p]!=work[q] )  // p, q have been assigned
		{
			--Ncluster;
			// if they have been assigned to a same index, do nothing.
			// so only consider that they are different.
			// Now what we do is to merge to clusters indexed by
			// work[p] & work[q] together.
			// We choose new index as the smaller one.
			// For the remaining indices, the ones > the larger one should
			// be decreasing 1, while the ones < the larger one keep,
			// and = the larger one to be the smaller one
			if ( work[p] < work[q] )
			{
				small = work[p]; large = work[q];
			} else {
				small = work[q]; large = work[p];
			}
			if ( Nnode == MaxNode )
			{
				for (j=0; j<MaxNode; j++)
					if ( work[j] > large )
						--work[j];
					else if ( work[j] == large )
						work[j] = small;
			} else {
				for (j=0; j<Nnode; j++)
					if ( work[cluster[j]] > large )
						--work[cluster[j]];
					else if ( work[cluster[j]] == large )
						work[cluster[j]] = small;
			}
		}
	}
	edges-=2*Nedge;

	if ( Nnode != MaxNode )
	{
		for (i=0; i<Nnode; i++)
			cluster[Nnode+i]=work[cluster[i]];
		free(work);
	}
	
	*pNnode=Nnode;
	*pNcluster=Ncluster;
}
