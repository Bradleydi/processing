#include "structure.h"
#include "dct.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void triangle_cluster(const char *file, double threshold, bool tracked)
{
	if ( threshold<=0 || threshold>=1.0 )
		 pERROR("threshold should be in (0,1).");

	char *gdffile=getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);
	free(gdffile);

	int i, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
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

	float *points=(float *)malloc(2*maxptcl*sizeof(float));
	POINTER_NULL(pointer);
	float *pointer1;
	int ntri;
	/* max_ triangle: 2*maxptcl
	 * so length of ``tri'': 3*max_nedge = 6*maxptcl 
	 * (3 vertices for each triangle)
	 */
	int *tri=(int *)malloc(6*maxptcl*sizeof(int));
	POINTER_NULL(tri);
	int *edges=(int *)malloc(6*maxptcl*sizeof(int));
	POINTER_NULL(edges);
	/* Actually, we should allocate 4*Nedge memory, i.e. 12*maxptcl,
	 * but usually, we have at most maxptcl nodes, and the copy
	 * process in getcluster asks for only 2*Nedge<=6*maxptcl memory,
	 * so I only allocate 6*maxptcl memory for the results from 
	 * function getcluster.
	 * Note that array tri will not be used once edges are found,
	 * so we just use tri as the container to save memory.
	*/
	int *TmpCluster=(int *)malloc(6*maxptcl*sizeof(int));
	POINTER_NULL(TmpCluster);
	//int *TmpCluster=tri;
	float *cluster=Calloc(float, size[1]);
	//for (i=0; i<size[1]; i++)
	//	cluster[i]=-1.0;
	POINTER_NULL(cluster);
	float **pindex;
	float *A, *B, *C; // triangle vertics A, B, C
	float area, perimeter, Ax, Ay, Bx, By, Cx, Cy;
	// area/permeter^2 isosceles right triangle
	const float ratio_ri=1.0/2.0/(2.0+sqrt(2.0))/(2.0+sqrt(2.0)) ; 
	// area/permeter^2 equilateral triangle-isosceles right triangle
	const float diff_ratio_eq_ri=sqrt(3.0)/36.0-ratio_ri; 	
	const float unfold_thres=diff_ratio_eq_ri*((float)threshold)+ratio_ri;
	int Nedge, Ncluster, Nnode, j, pp, qq, ss;
	unsigned char new_edge_pq, new_edge_ps;
	int *check;
	//unsigned char *boundary=Malloc(unsigned char, maxptcl);
	//POINTER_NULL(boundary);
	//unsigned char *lattice=Malloc(unsigned char, maxptcl);
	//POINTER_NULL(lattice);
	unsigned char *right=Malloc(unsigned char, 2*maxptcl); // is right tri?
	POINTER_NULL(right);
	int *boundary=TmpCluster+2*maxptcl;
	int *clusterindex=TmpCluster+3*maxptcl;
	/*
	 *
I found that if I only store whether a particle is crystal-like, it's
not enough to get cluster boundary. The reason is very simple, but hard to
relize, that there are boundaries between two crystal clusters.
In practice, there can be no 'liquid'-like particles between two crystal 
clusters, so my old method will fail for this case.
There are 2 methods to solve this problem:
1. find clusters first, and then find boundary that nodes in a Delaunay 
triangle have different cluster indices.
This method needs storage of triangulation after finding cluster.
2. store whether a Delaunay triangle is right triangle, since there are right
triangles with all nodes denoted as crystal-like, since they are between
two crystal clusters. So those triangles are actual boundaries. The method is
if it's a right triangle, all crystal-like nodes are boundaries.
Method 2 is much easy to implimented to recent program, and uses less memory.
So I would like to use method 2.
>>> After some try, I found that method 2 is not right, since there are 
triangles that are right, but with three vertices in same cluster.
By method 2, the 3 vertices will be identified as boundaries, even they 
are inner points.
So the only choice is method 1.
To implemented method 1, notice that the array TmpCluster should be 
allocated to 6*maxptcl memory. Even all of them may be used in function
getcluster(), however, after the return, at most 2*maxptcl memory will
be used. So we can use the other 4*maxptcl memory to store the boundary
information, and other informations.
Also, a tmpCluster information for current frame should also be stored,
to make the procedure that finding cluster index given a node much easier.
*/

	//for (t=0; t<1;t++)
	for (t=0; t<total_t; t++)
	{
		pointer1=points;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			*(pointer1++)=*pointer;
			*(pointer1++)=*(pointer+1);
		}
		// get triangles of Delaunay triangluation
		dct_tri(ptcl_frame[t], points, ntri, tri);

		/*==========================================================*/
		// check those triangles
		pindex=index+cum_ptcl_frame[t];
		Nedge=0;
		for (i=0; i<ntri; i++)
		{
			pp=*(tri++);
			qq=*(tri++);
			ss=*(tri++);
			A=*(pindex+pp);
			B=*(pindex+qq);
			C=*(pindex+ss);
			/*
			A=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex A
			B=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex B
			C=index[cum_ptcl_frame[t] + (*(tri++))]; // triangle vertex C
			*/

			Ax=*A; Ay=*(A+1);
			Bx=*B; By=*(B+1);
			Cx=*C; Cy=*(C+1);
			perimeter=sqrt((Ax-Bx)*(Ax-Bx)+(Ay-By)*(Ay-By));
			perimeter+=sqrt((Ax-Cx)*(Ax-Cx)+(Ay-Cy)*(Ay-Cy));
			perimeter+=sqrt((Cx-Bx)*(Cx-Bx)+(Cy-By)*(Cy-By));

			area=((Ax-Cx)*(By-Ay)-(Ax-Bx)*(Cy-Ay))/2;
			if (area<0)
				area=-area;

			if ( area/perimeter/perimeter <= unfold_thres )
				right[i]=1;
			else
			{
				right[i]=0;
				//lattice[pp]=1; lattice[qq]=1; lattice[ss]=1;
				// is an equilateral triangle
				new_edge_pq=1;
				new_edge_ps=1;
				check=edges;
				for (j=0; j<Nedge; j++)
				{
					if ( pp == *check && qq == *(check+1) )
					{
						new_edge_pq=0;
						if ( new_edge_ps==0 ) break;
					}
					else if ( qq == *check && pp == *(check+1) )
					{
						new_edge_pq=0;
						if ( new_edge_ps==0 ) break;
					}

					if ( pp == *check && ss == *(check+1) )
					{
						new_edge_ps=0;
						if ( new_edge_pq==0 ) break;
					}
					else if ( ss == *check && pp == *(check+1) )
					{
						new_edge_ps=0;
						if ( new_edge_pq==0 ) break;
					}
					check+=2;
				}
				if ( new_edge_pq )
				{
					edges[2*Nedge]=pp;
					edges[2*(Nedge++)+1]=qq;
				}
				if ( new_edge_ps )
				{
					edges[2*Nedge]=pp;
					edges[2*(Nedge++)+1]=ss;
				}
			}
		}
		tri-=3*ntri;

		getcluster(Nedge, edges, &Nnode, &Ncluster, TmpCluster);
		
		for (i=0; i<ptcl_frame[t]; i++)
		{
			boundary[i]=0;
			clusterindex[i]=-1;
			// boundary == 1 interface
		}
		for (i=0; i<Nnode; i++)
			clusterindex[TmpCluster[i]]=TmpCluster[Nnode+i];
		
		// find boundary
		for (i=0; i<ntri; i++)
		{
			pp=*(tri++);
			qq=*(tri++);
			ss=*(tri++);
			//if ( lattice[pp]*lattice[qq]*lattice[ss] == 0 )
			if ( right[i] && ( clusterindex[pp] != clusterindex[qq] ||
					 clusterindex[pp] != clusterindex[ss] ) )
			{
				/*
				if ( lattice[pp] ) boundary[pp]=1;
				if ( lattice[qq] ) boundary[qq]=1;
				if ( lattice[ss] ) boundary[ss]=1;
				*/
				if ( clusterindex[pp] != -1 ) boundary[pp]=1;
				if ( clusterindex[qq] != -1 ) boundary[qq]=1;
				if ( clusterindex[ss] != -1 ) boundary[ss]=1;
			}
		}
		tri-=3*ntri;

		for (i=0; i<Nnode; i++)
		{
			if ( boundary[TmpCluster[i]] )
				cluster[(*(pindex+TmpCluster[i])-p_ptr)/size[0]]=
					(float)(-TmpCluster[Nnode+i]-1);
			else
				cluster[(*(pindex+TmpCluster[i])-p_ptr)/size[0]]=
					(float)(TmpCluster[Nnode+i]+1);
		}
	}
	free(edges);
	free(points);
	free(tri);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	//free(boundary);
	free(TmpCluster);
	free(index);

	char *outfile = getfilename (file, "_cluster.gdf") ;
	FILE *f=fopen(outfile, "wb");
	FILE_NULL(f, outfile);

	int headersize [6];
	headersize[0]=82991;
	headersize[1]=2;
	headersize[2]=size[0]+1;
	headersize[3]=size[1];
	headersize[4]=4; // IDL float type code
	headersize[5]=(size[0]+1)*size[1];

	if ( (int)fwrite(headersize, sizeof(headersize[0]), 6, f) != 6 )
	{
		fprintf(stderr, "# Error: write to file '%s' failed!\n", outfile);
		exit (1);
	}
	
	t=size[0]+1;
	float *tmp=Malloc(float, t);
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		for (j=0; j<ti; j++)
			tmp[j]=*(pointer+j);
		tmp[ti]=cluster[i];
		// change time back
		tmp[ti+1]=(float)((int)(*(pointer+ti))+tmin);
		for (j=ti+1; j<size[0]; j++)
			tmp[j+1]=*(pointer+j);
		pointer+=size[0];
		
		if ( (int)fwrite(tmp, sizeof(tmp[0]), t, f) != t )
		{
			fprintf(stderr, "# Error: write to file '%s' failed!\n", outfile);
			exit (1);
		}
	}
	fclose(f);
	free(outfile);
	free(tmp);
	free(cluster);
	delete [] size;
}
