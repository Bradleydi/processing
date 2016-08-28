#include "dct.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>

void dct(int &np, float *points, int &nedge, int *edges)
{
	cardinal n=np;

	edge *l_cw, *r_ccw;
	point **p_sorted, **p_temp;
	INDEX i;

  	alloc_memory(n);
	
	for (i=0; i<n; i++)
	{
		p_array[i].x=*(points++);
		p_array[i].y=*(points++);
	}


	/* Initialise entry edge pointers. */
	for (i = 0; i < n; i++)
		p_array[i].entry_pt = NULL;

		/* Sort. */
	p_sorted = (point **)malloc((unsigned)n*sizeof(point *));
	if (p_sorted == NULL)
		panic("triangulate: not enough memory\n");
	p_temp = (point **)malloc((unsigned)n*sizeof(point *));
	if (p_temp == NULL)
		panic("triangulate: not enough memory\n");
	for (i = 0; i < n; i++)
		p_sorted[i] = p_array + i;

	merge_sort(p_sorted, p_temp, 0, n-1);
  
	free((char *)p_temp);


	/* Triangulate. */
	divide(p_sorted, 0, n-1, &l_cw, &r_ccw);

	free((char *)p_sorted);


	// store edges
	edge *e_start, *e;
	point *u, *v;

	nedge=0;
	for (i=0; i<n; i++)
	{
		u=&p_array[i];
		e_start=e=u->entry_pt;
		do
		{
			v=Other_point(e, u);
			if (u<v)
			{
				++nedge;
				*(edges++)=u-p_array;
				*(edges++)=v-p_array;
			}
			e=Next(e, u);
		}while(!Identical_refs(e, e_start));
	}

	free_memory();
}


void dct_RNG(int &np, float *points, int &nedge, int *edges) 
{
	int *tmp_edges=(int *)malloc(6*np*sizeof(int));
	if (tmp_edges==NULL)
		panic("dct_RNG: int *tmp_edges!\n");
	
	int oldnedge, i;
	dct(np, points, oldnedge, tmp_edges);
	
	//==========================================================
	// count neighbor
	int *neighbors=(int *)calloc(np, sizeof(int));
	if (neighbors==NULL)
		panic("dct: int *neighbors!\n");
	
	for (i=0; i<2*oldnedge; i++)
		++neighbors[tmp_edges[i]];
	
	int *cum_neighbors=(int *)malloc((np+1)*sizeof(int));
	if (cum_neighbors==NULL)
		panic("dct: int *cum_neighbors!\n");
	cum_neighbors[0]=0;
	for (i=0; i<np; i++)
		cum_neighbors[i+1]=cum_neighbors[i]+neighbors[i];

	free(neighbors);

	int *neighbor_index=(int *)malloc(cum_neighbors[np]*sizeof(int));
	if (neighbor_index==NULL)
		panic("dct: int *neighbor_index!\n");
	int *index_now=(int *)calloc(np, sizeof(int));
	if (index_now==NULL)
		panic("dct: int *index_now!\n");
	int p, q;
	for (i=0; i<oldnedge; i++)
	{
		p=tmp_edges[2*i];
		q=tmp_edges[2*i+1];
		neighbor_index[cum_neighbors[p]+index_now[p]]=q;
		++index_now[p];
		neighbor_index[cum_neighbors[q]+index_now[q]]=p;
		++index_now[q];
	}
	
	free(index_now);

	int o, j;
	float d, px, py, qx, qy, ox, oy, pq, qo, op;
	char * notgood=(char *)calloc(oldnedge, sizeof(char));
	if (notgood==NULL)
		panic("dct: char *notgood!\n");

	nedge=oldnedge;
	for (i=0; i<oldnedge; i++)
	{
		p=tmp_edges[2*i];
		q=tmp_edges[2*i+1];

		px=points[2*p];
		py=points[2*p+1];
		qx=points[2*q];
		qy=points[2*q+1];

		d=px-qx;
		pq=d*d;
		d=py-qy;
		pq+=d*d;

		for (j=cum_neighbors[p]; j<cum_neighbors[p+1]; j++)
		{
			o=neighbor_index[j];
			ox=points[2*o];
			oy=points[2*o+1];

			d=qx-ox;
			qo=d*d;
			d=qy-oy;
			qo+=d*d;

			d=ox-px;
			op=d*d;
			d=oy-py;
			op+=d*d;

			if (pq>qo && pq>op)
				notgood[i]=1;
		}
		
		for (j=cum_neighbors[q]; j<cum_neighbors[q+1]; j++)
		{
			o=neighbor_index[j];
			ox=points[2*o];
			oy=points[2*o+1];

			d=qx-ox;
			qo=d*d;
			d=qy-oy;
			qo+=d*d;

			d=ox-px;
			op=d*d;
			d=oy-py;
			op+=d*d;

			if (pq>qo && pq>op)
				notgood[i]=1;
		}
		if (notgood[i])
			--nedge;
		
	}

	
	free(cum_neighbors);
	free(neighbor_index);

	int *jj=edges;
	for (i=0; i<oldnedge; i++)
		if (notgood[i]==0)
		{
			*(jj++)=tmp_edges[2*i];
			*(jj++)=tmp_edges[2*i+1];
		}

	free(notgood);
	free(tmp_edges);
}


void dct_draw(int &np, float *points, int &nedge, int *edges,
              const char *filename, unsigned char defect)
{
	int defect1=5, defect2=7;
	if (defect==35)
	{
		defect1=3;
		defect2=5;
	}
	else if (defect!=57)
		panic("dct: unknown defect!\n");

	//==========================================================
	// count neighbor
	int *neighbors=(int *)calloc(np, sizeof(int));
	if (neighbors==NULL)
		panic("dct: int *neighbors!\n");
	
	int i;
	for (i=0; i<2*nedge; i++)
		++neighbors[edges[i]];
	
	//==========================================================
	// draw
	FILE * out=fopen(filename, "w");

	// HEADER
	fprintf(out, "#FIG 3.2 Produced by colloid\n");
	fprintf(out, "Landscape\n");
	fprintf(out, "Center\n");
	fprintf(out, "Inches\n");
	fprintf(out, "Letter\n");
	fprintf(out, "100.00\n");
	fprintf(out, "Single\n");
	fprintf(out, "-2\n");
	fprintf(out, "1200 2\n");

	// boundary ?
	

	// draw node
	int factor=100;
	int ax, ay;
	int r=300;
	for (i=0; i<np; i++)
	{
		if (neighbors[i]==defect1)
		{
			ax=(int)(factor*points[2*i]);
			ay=(int)(factor*points[2*i+1]);
			// blue
			fprintf(out, "1 3 0 3 1 1 50 -1 20 0.000 1 0.0000 ");
			fprintf(out, "%d %d %d %d %d %d %d %d\n", ax, ay, r, r, ax, ay, ax+r, ay+r);
		}
		else if (neighbors[i]==defect2)
		{
			ax=(int)(factor*points[2*i]);
			ay=(int)(factor*points[2*i+1]);
			// red
			fprintf(out, "1 3 0 3 4 4 50 -1 20 0.000 1 0.0000 ");
			fprintf(out, "%d %d %d %d %d %d %d %d\n", ax, ay, r, r, ax, ay, ax+r, ay+r);
		}
	}
	free(neighbors);

	
	// draw edge
	int v, w;
	int vx, vy, wx, wy;

	for (i=0; i<nedge; i++)
	{
		fprintf(out, "2 1 0 2 0 7 51 -1 -1 0.000 0 0 -1 0 0 2\n\t");
		v=edges[2*i];
		w=edges[2*i+1];
		vx=(int)(factor*points[2*v]+0.5);
		vy=(int)(factor*points[2*v+1]+0.5);
		wx=(int)(factor*points[2*w]+0.5);
		wy=(int)(factor*points[2*w+1]+0.5);
		fprintf(out, "%d %d %d %d\n", vx, vy, wx, wy);
	}
	fclose(out);
}
