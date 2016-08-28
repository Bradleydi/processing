#include "dct.h"

#include <cstdio>
#include <cstdlib>


void dct_tri(int &np, float *points, int &ntri, int *tri)
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

	// store triangles
	edge *e_start, *e, *next;
	point *u, *v, *w;
	point *t;

	ntri=0;
	for (i = 0; i < n; i++) 
	{
		u = &p_array[i];
	 	e_start = e = u->entry_pt;
		do
		{
			v = Other_point(e, u);
			if (u < v) 
			{
				next = Next(e, u);
				w = Other_point(next, u);
				if (u < w)
					if (Identical_refs(Next(next, w), Prev(e, v))) 
					{  
						/* Triangle. */
						if (v > w) { t = v; v = w; w = t; }
						++ntri;
						*(tri++)=u-p_array;
						*(tri++)=v-p_array;
						*(tri++)=w-p_array;
					}
			}
			/* Next edge around u. */
			e = Next(e, u);
		} while (!Identical_refs(e, e_start));
	}
	tri-=3*ntri;
	free_memory();
}
