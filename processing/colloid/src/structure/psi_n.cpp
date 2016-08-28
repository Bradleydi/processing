// version 0.9.2

#include "structure.h"
#include "dct.h"
#include "miscellaneous.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;


void psi_n(colloid_base &p, int n, 
		   double *Re_psi, double *Im_psi, bool tracked)
{
	int i, j, t;

	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
		--ti;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");
	float tt=0;

	float *pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)tt + 1;

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	int *cum_ptcl_frame=(int *)malloc((total_t+1)*sizeof(int));
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

	float **index=(float **)malloc(size[1]*sizeof(float *));
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	float *points=(float *)malloc(2*maxptcl*sizeof(float));
	float *pointer1;
	int nedge, node1, node2, ind;
	/* max_ nedge: 3*maxptcl
	 * so length of ``edges'': 2*max_nedge = 6*maxptcl 
	 * (two vertices for each edge)
	 * all neighbor ids = 2*nedge
	 */
	int *edges=(int *)malloc(6*maxptcl*sizeof(int));
	int *neighbors=(int *)malloc(maxptcl*sizeof(int));
	int *cum_neighbors=(int *)malloc((maxptcl+1)*sizeof(int));
	index_now=(int *)malloc(maxptcl*sizeof(int));
	int *neighbors_id=(int *)malloc(6*maxptcl*sizeof(int));
	double re_psi, im_psi, dx, dy, theta;
	float **pindex;
	for (t=0; t<total_t; t++)
	{
		pointer1=points;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			*(pointer1++)=*pointer;
			*(pointer1++)=*(pointer+1);
		}
		if (n==6)
			dct(ptcl_frame[t], points, nedge, edges);
		else if (n==4)
			dct_RNG(ptcl_frame[t], points, nedge, edges);

		for (i=0; i<ptcl_frame[t]; i++)
		{
			neighbors[i]=0;
			index_now[i]=0;
		}
		for (i=0; i<2*nedge; i++)
			++neighbors[edges[i]];

		cum_neighbors[0]=0;
		for (i=0; i<ptcl_frame[t]; i++)
			cum_neighbors[i+1]=cum_neighbors[i]+neighbors[i];
		
		for (i=0; i<nedge; i++)
		{
			node1=edges[2*i];
			node2=edges[2*i+1];
			neighbors_id[cum_neighbors[node1]+(index_now[node1]++)]=node2;
			neighbors_id[cum_neighbors[node2]+(index_now[node2]++)]=node1;
		}
		
		/*==========================================================*/
		/* psi_n */
		pindex=index+cum_ptcl_frame[t];
		for (i=0; i<ptcl_frame[t]; i++)
		{
			pointer=*(pindex+i);
			re_psi=0; im_psi=0;
			for (j=cum_neighbors[i]; j<cum_neighbors[i+1]; j++)
			{
				pointer1=*(pindex+neighbors_id[j]);
				dx=(double)(*pointer1-*pointer);
				dy=(double)(*(pointer1+1)-*(pointer+1));
				theta=n*atan2(dy, dx);
				re_psi+=cos(theta);
				im_psi+=sin(theta);
			}
			ind=(pointer-p_ptr)/size[0];
			Re_psi[ind]=re_psi/neighbors[i];
			Im_psi[ind]=im_psi/neighbors[i];
		}
	}

	free(neighbors);
	free(cum_neighbors);
	free(neighbors_id);
	free(index_now);
	free(points);
	free(edges);
	free(ptcl_frame);
	free(cum_ptcl_frame);

	delete [] size;
}


void draw_psi_n(colloid_base &p, int n, const char *filename, bool tracked)
{
	int i, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	double *Re_psi=(double *)malloc(size[1]*sizeof(double));
	double *Im_psi=(double *)malloc(size[1]*sizeof(double));

	psi_n(p, n, Re_psi, Im_psi, tracked);

	
	/* |psi|^2 */
	double *abs_psi_2=(double *)malloc(size[1]*sizeof(double));
	for (i=0; i<size[1]; i++)
		abs_psi_2[i]=Re_psi[i]*Re_psi[i]+Im_psi[i]*Im_psi[i];
	
	free(Re_psi);
	free(Im_psi);

	int nstr=0;
	while (filename[nstr]!='\0')
		++nstr;

	/* filename_00000.fig\0, so nstr+11 */
	char *figfile=(char *)malloc((nstr+11)*sizeof(char));
	for (i=0; i<nstr; i++)
		figfile[i]=filename[i];

	figfile[nstr]='_';
	figfile[nstr+6]='.';
	figfile[nstr+7]='f';
	figfile[nstr+8]='i';
	figfile[nstr+9]='g';
	figfile[nstr+10]='\0';

	
	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
		--ti;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");
	float tt=0;

	float *pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)tt + 1;
	
	char str [5];
	int factor=100;
	int ax, ay;
	int r=300;

	FILE *out;

	double threshold=0.25; /* 0.5*0.5 */

	for (t=0; t<total_t; t++)
	{
		sprintf(str, "%05d", t);
		for (i=0; i<5; i++)
		figfile[nstr+i+1]=str[i];

		out=fopen(figfile, "w");
		/* HEADER */
		fprintf(out, "#FIG 3.2 Produced by colloid\n");
	    fprintf(out, "Landscape\n");
	    fprintf(out, "Center\n");
	    fprintf(out, "Inches\n");
	    fprintf(out, "Letter\n");
	    fprintf(out, "100.00\n");
	    fprintf(out, "Single\n");
    	fprintf(out, "-2\n");
	    fprintf(out, "1200 2\n");

		pointer=p_ptr;
		for (i=0; i<size[1]; i++)
		{
			if ((int)(*(pointer+ti))==t)
			{
				ax=(int)(factor*(*pointer));
				ay=(int)(factor*(*(pointer+1)));

				/* crystal-like particle */
				/* blue */
				if (abs_psi_2[i]>threshold)
				{
					fprintf(out, "1 3 0 3 1 1 50 -1 20 0.000 1 0.0000 ");
		            fprintf(out, "%d %d %d %d %d %d %d %d\n", 
							ax, ay, r, r, ax, ay, ax+r, ay+r);
				}
				/* liquid-like particle */
				/* red */
				else
				{
					fprintf(out, "1 3 0 3 4 4 50 -1 20 0.000 1 0.0000 ");
		            fprintf(out, "%d %d %d %d %d %d %d %d\n", 
							ax, ay, r, r, ax, ay, ax+r, ay+r);
				}
			}
			pointer+=size[0];
		}

		fclose(out);
	}

	free(figfile);
	free(abs_psi_2);

	delete [] size;
}


void draw_psi_n_with_bond(colloid_base &p, int n, 
		const char *filename, int frame_number, bool tracked)
{
	int i, j, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();


	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
		--ti;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");
	float tt=0;

	float *pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)tt + 1;
	
	if (frame_number>total_t)
		pERROR("set frame_number larger than total_t!");

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	int *cum_ptcl_frame=(int *)malloc((total_t+1)*sizeof(int));
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

	float **index=(float **)malloc(size[1]*sizeof(float *));
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	float *points=(float *)malloc(2*maxptcl*sizeof(float));
	float *pointer1;
	int nedge, node1, node2;
	/* max_ nedge: 3*maxptcl
	 * so length of ``edges'': 2*max_nedge = 6*maxptcl 
	 * (two vertices for each edge)
	 * all neighbor ids = 2*nedge
	 */
	int *edges=(int *)malloc(6*maxptcl*sizeof(int));
	int *neighbors=(int *)malloc(maxptcl*sizeof(int));
	int *cum_neighbors=(int *)malloc((maxptcl+1)*sizeof(int));
	index_now=(int *)malloc(maxptcl*sizeof(int));
	int *neighbors_id=(int *)malloc(6*maxptcl*sizeof(int));
	double re_psi, im_psi, dx, dy, theta;
	float **pindex;

	int v, w;
	int vx, vy, wx, wy;

	/* |psi|^2 */
	double abs_psi_2;

	int nstr=0;
	while (filename[nstr]!='\0')
		++nstr;

	/* filename_00000.fig\0, so nstr+11 */
	char *figfile=(char *)malloc((nstr+11)*sizeof(char));
	for (i=0; i<nstr; i++)
		figfile[i]=filename[i];

	figfile[nstr]='_';
	figfile[nstr+6]='.';
	figfile[nstr+7]='f';
	figfile[nstr+8]='i';
	figfile[nstr+9]='g';
	figfile[nstr+10]='\0';
	
	char str [5];
	int factor=100;
	int ax, ay;
	int r=300;

	FILE *out;

	double threshold=0.25; /* 0.5*0.5 */
	
	for (t=0; t<frame_number; t++)
	{
		sprintf(str, "%05d", t);
		for (i=0; i<5; i++)
		figfile[nstr+i+1]=str[i];

		out=fopen(figfile, "w");
		/* HEADER */
		fprintf(out, "#FIG 3.2 Produced by colloid\n");
	    fprintf(out, "Landscape\n");
	    fprintf(out, "Center\n");
	    fprintf(out, "Inches\n");
	    fprintf(out, "Letter\n");
	    fprintf(out, "100.00\n");
	    fprintf(out, "Single\n");
    	fprintf(out, "-2\n");
	    fprintf(out, "1200 2\n");


		pointer1=points;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			*(pointer1++)=*pointer;
			*(pointer1++)=*(pointer+1);
		}
		if (n==6)
			dct(ptcl_frame[t], points, nedge, edges);
		else if (n==4)
			dct_RNG(ptcl_frame[t], points, nedge, edges);

		for (i=0; i<ptcl_frame[t]; i++)
		{
			neighbors[i]=0;
			index_now[i]=0;
		}
		for (i=0; i<2*nedge; i++)
			++neighbors[edges[i]];

		cum_neighbors[0]=0;
		for (i=0; i<ptcl_frame[t]; i++)
			cum_neighbors[i+1]=cum_neighbors[i]+neighbors[i];
		
		for (i=0; i<nedge; i++)
		{
			node1=edges[2*i];
			node2=edges[2*i+1];
			neighbors_id[cum_neighbors[node1]+(index_now[node1]++)]=node2;
			neighbors_id[cum_neighbors[node2]+(index_now[node2]++)]=node1;
		}
		
		/*==========================================================*/
		/* psi_n and draw node */
		pindex=index+cum_ptcl_frame[t];
		for (i=0; i<ptcl_frame[t]; i++)
		{
			pointer=*(pindex+i);
			re_psi=0; im_psi=0;
			for (j=cum_neighbors[i]; j<cum_neighbors[i+1]; j++)
			{
				pointer1=*(pindex+neighbors_id[j]);
				dx=(double)(*pointer1-*pointer);
				dy=(double)(*(pointer1+1)-*(pointer+1));
				theta=n*atan2(dy, dx);
				re_psi+=cos(theta);
				im_psi+=sin(theta);
			}
			abs_psi_2=(re_psi*re_psi+im_psi*im_psi)/neighbors[i]/neighbors[i];

			if (abs_psi_2<threshold)
			{
				/* liquid-like particle */
				/* red */
				ax=(int)(factor*(*pointer));
				ay=(int)(factor*(*(pointer+1)));
				fprintf(out, "1 3 0 3 4 4 50 -1 20 0.000 1 0.0000 ");
		        fprintf(out, "%d %d %d %d %d %d %d %d\n", 
						ax, ay, r, r, ax, ay, ax+r, ay+r);
			}

		}

		/* draw bonds */
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

	free(neighbors);
	free(cum_neighbors);
	free(neighbors_id);
	free(index_now);
	free(points);
	free(edges);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	free(figfile);

	delete [] size;

}


void psi_n_statistics(colloid_base &p, int n, 
		double threshold, bool tracked)
{
	int i, j, t;
	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();


	const int dim=2;
	int ti=size[0]-1;
	if (tracked)
		--ti;
	if (ti<dim)
		pERROR("dimension is not 2, or not a tracked data!");
	float tt=0;

	float *pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)tt + 1;
	

	int *ptcl_frame=(int *)calloc(total_t, sizeof(int));
	pointer=p_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	int *cum_ptcl_frame=(int *)malloc((total_t+1)*sizeof(int));
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

	float **index=(float **)malloc(size[1]*sizeof(float *));
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	float *points=(float *)malloc(2*maxptcl*sizeof(float));
	float *pointer1;
	int nedge, node1, node2;
	/* max_ nedge: 3*maxptcl
	 * so length of ``edges'': 2*max_nedge = 6*maxptcl 
	 * (two vertices for each edge)
	 * all neighbor ids = 2*nedge
	 */
	int *edges=(int *)malloc(6*maxptcl*sizeof(int));
	int *neighbors=(int *)malloc(maxptcl*sizeof(int));
	int *cum_neighbors=(int *)malloc((maxptcl+1)*sizeof(int));
	index_now=(int *)malloc(maxptcl*sizeof(int));
	int *neighbors_id=(int *)malloc(6*maxptcl*sizeof(int));
	double re_psi, im_psi, dx, dy, theta;
	float **pindex;

	/* |psi|^2 */
	double abs_psi_2;

	/* count: number of non-crystal-like particles */
	unsigned long count=0, sum=0, sum2=0;

	double threshold2=threshold*threshold; 

	printf("# average particle    :\t%6.2f\n", (double)(size[1])/total_t);
	printf("# total frame         :\t%d\n", total_t);
	printf("# threshold           :\t%f\n", threshold);
	
	for (t=0; t<total_t; t++)
	{
		count=0;

		pointer1=points;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			*(pointer1++)=*pointer;
			*(pointer1++)=*(pointer+1);
		}
		if (n==6)
			dct(ptcl_frame[t], points, nedge, edges);
		else if (n==4)
			dct_RNG(ptcl_frame[t], points, nedge, edges);

		for (i=0; i<ptcl_frame[t]; i++)
		{
			neighbors[i]=0;
			index_now[i]=0;
		}
		for (i=0; i<2*nedge; i++)
			++neighbors[edges[i]];

		cum_neighbors[0]=0;
		for (i=0; i<ptcl_frame[t]; i++)
			cum_neighbors[i+1]=cum_neighbors[i]+neighbors[i];
		
		for (i=0; i<nedge; i++)
		{
			node1=edges[2*i];
			node2=edges[2*i+1];
			neighbors_id[cum_neighbors[node1]+(index_now[node1]++)]=node2;
			neighbors_id[cum_neighbors[node2]+(index_now[node2]++)]=node1;
		}
		
		/*==========================================================*/
		/* psi_n and draw node */
		pindex=index+cum_ptcl_frame[t];
		for (i=0; i<ptcl_frame[t]; i++)
		{
			pointer=*(pindex+i);
			re_psi=0; im_psi=0;
			for (j=cum_neighbors[i]; j<cum_neighbors[i+1]; j++)
			{
				pointer1=*(pindex+neighbors_id[j]);
				dx=(double)(*pointer1-*pointer);
				dy=(double)(*(pointer1+1)-*(pointer+1));
				theta=n*atan2(dy, dx);
				re_psi+=cos(theta);
				im_psi+=sin(theta);
			}
			abs_psi_2=(re_psi*re_psi+im_psi*im_psi)/neighbors[i]/neighbors[i];

			if (abs_psi_2<threshold2)
			{
				/* liquid-like particle */
				/* red */
				++count;
			}

		}

		sum+=count;
		sum2+=count*count;

		printf("%d\t %lu\n", t, count);

	}

	free(neighbors);
	free(cum_neighbors);
	free(neighbors_id);
	free(index_now);
	free(points);
	free(edges);
	free(ptcl_frame);
	free(cum_ptcl_frame);

	double ave=(double)(sum)/total_t;
	printf("# average non-crystal particles:\t%6.6f\t+/-\t%6.6f\n", ave,
			sqrt((double)(sum2)/total_t-ave*ave));

	delete [] size;

}
