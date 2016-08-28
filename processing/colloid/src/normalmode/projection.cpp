#include "normalmode.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
void projection(const char *file, double a0, bool Remove_Drift)
{
	int i, k, t, p;
	// read gdf 
	char *gdffilename=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffilename);
	free(gdffilename);

	int *size = ptid.get_size();
	float *ptid_ptr = ptid.get_array_pointer();
	
	const int total_p = (int)ptid_ptr[size[2]-1]+1;

	const int dim = 2;
	if (dim+2>size[0])
		pERROR("Not a 2D data or tracked data!");
	

	const int ti=size[0]-2; //time index in the array
	const int ii=ti+1; // particle index in the array

	float *pointer=ptid_ptr+ti;
	float tt = *pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( tt < *pointer )
			tt = *pointer;
	}
	const int total_t = (int)tt + 1 ;
	
	if ( total_t < 2*total_p )
	{
		fprintf(stderr, 
			"# Error: there should be enough independent configurations!\n");
		fprintf(stderr, "#        frame_number >= 2*particle_number\n");
		exit (1);
	}
	
	if ((total_p*total_t)!=size[1])
		pERROR("bad tracking, there are broken trajectories.");

	// rescale 
	if (a0>0)
	{
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			for (k=0; k<dim; k++)
				*(pointer+k)/=a0;
			pointer+=size[0];
		}
	}
	
	if (Remove_Drift)
		remove_drift(ptid, dim);
	

	//=====================================================================
	// get mean position
	float * mp_ptr = (float *)calloc(dim*total_p, sizeof(float));
	POINTER_NULL(mp_ptr);

	float *pointer1;
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		for (k=0; k<dim; k++)
			*(pointer1+k)+=*(pointer+k);
		pointer+=size[0];
	}

	// mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
		for (k=0; k<dim; k++)
			*(pointer++)/=total_t;

	
	//=====================================================================
	// get displacement
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(pointer1+k);
		pointer+=size[0];
	}

	free(mp_ptr);
	
	//=====================================================================
	// get normal modes from .ev file
	char *evfilename=getfilename(file, ".ev");
	int Gn;
	double *E;
	double *G;
	if ( readev(Gn, &E, &G, evfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", 
				evfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	

	if ( Gn != 2*total_p )
		pERROR("reading files wrong!");
	
	
	double *omega=(double *)malloc( (Gn-2)*sizeof(double) );
	POINTER_NULL(omega);

	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);

	free(E);

	//===================================================
	double *dx=(double *) malloc (Gn*sizeof(double));
	POINTER_NULL(dx);
	double *projection=(double *)calloc((Gn-2), sizeof(double));
	POINTER_NULL(projection);
	// assuming there is  no particle loss in the data,
	// so that i'th particle in t frame is at i*total_t+t
	
	double *G_pointer;
	
	double *dxG=(double *) malloc (Gn*sizeof(double));
	POINTER_NULL(dxG);
	double *dx2=(double *) malloc (Gn*sizeof(double));
	POINTER_NULL(dx2);

	double pj;
	
	const int total_t_size0=total_t*size[0];

	for (t=0; t<total_t; t++)
	{
		pointer=ptid_ptr+t*size[0];
		
		for (p=0; p<total_p; p++)
		{
			// get dx
			*(dx++)=(double)(*pointer); // dx
			*(dx++)=(double)(*(pointer+1)); // dy
			
			pointer+=total_t_size0;
		}
		dx -= Gn;

		G_pointer = G+2*Gn;
		
		for (i=0; i<Gn-2; i++)
		{
			for (p=0; p<total_p; p++)
			{
				*(dxG++) = (*dx) * (*(G_pointer++));  // dx*G
				*(dx2++) = (*dx) * (*dx);         // dx*dx
				++dx;
				
				*(dxG++) = (*dx) * (*(G_pointer++));  // dy*G
				*(dx2++) = (*dx) * (*dx);         // dy*dy
				++dx;
			}
			dx  -= Gn;
			dxG -= Gn;
			dx2 -= Gn;
			pj = pairwise(Gn, dxG)/sqrt(pairwise(Gn, dx2));
			projection[i] += pj*pj/total_t;
		}
	}

	free (dx);
	free (dx2);
	free (dxG);
	free (G);

	//===================================================
	// print out the data
	char *pjfilename = getfilename(file, ".pj");
	FILE *pjfile=fopen(pjfilename, "w");
	FILE_NULL(pjfile, pjfilename);

	fprintf(pjfile, "# Projection of displacement\n");
	fprintf(pjfile, "# omega projection cum_projection\n");

	double cum=0;
	for (i=Gn-3; i>=0; i--)
	{
		cum += projection[i];
		fprintf(pjfile, "%6.6f %6.6f %6.6f\n", omega[i], projection[i], cum);
	}
	fclose (pjfile);
	
	free(omega);
	free(projection);

	//=========================================================
	// plot by GNUplot, using pipe
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set logscale x\n"
				   "set xlabel 'omega'\n"
				   "set ylabel 'cum_projection'\n");

	fprintf(pipe, "plot '%s' using 1:3\n", pjfilename);

	fclose(pipe);

	free(pjfilename);

	delete [] size;
}
