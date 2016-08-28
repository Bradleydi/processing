#include "dynamics.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

void vibration_Gaussian(const char *file, bool Remove_drift)
{
	int i, j, k, t, p;
	
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);

	int *size=ptid.get_size();
	const int ti=size[0]-2;
	const int dim=2;

	float *ptid_ptr=ptid.get_array_pointer();
	
	
	// get total_p & total_t
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	float* pointer=ptid_ptr+ti;
	float T=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if ( T < *pointer )
			T = *pointer;
	}
	const int total_t=(int)(T)+1;

	
	// remove_drift
	if (Remove_drift)
		remove_drift(ptid, dim);


	/*
	//rescale
	if (a0>0)
	{
		float aa0=(float)a0;
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			for (k=0; k<dim; k++)
				*(pointer+k)/=aa0;
			pointer+=size[0];
		}
	}*/

	// get index
	j=total_t*total_p;
	if ( j > size[1]*2 )
		pERROR("bad tracking!");

	int * index=(int *)malloc(j*sizeof(int));
	POINTER_NULL(index);
	
	for (i=0; i<j; i++)
		index[i]=-1;
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		index[(int)(*pointer)+(int)(*(pointer+1))*total_t]=i*size[0];
		pointer+=size[0];
	}

	// get mean position
	
	printf("# get mean position ...\n");

	double *mp=(double *)calloc(dim*total_p, sizeof(double));
	POINTER_NULL(mp);

	int count;
	for (p=0; p<total_p; p++)
	{
		count=0;
		for (t=0; t<total_t; t++)
		{
			i=index[t+p*total_t];
			if ( i != -1 )
			{
				++count;
				pointer=ptid_ptr+i;
				for (k=0; k<dim; k++)
				{
					mp[p*dim+k] += (double) (*(pointer+k));
				}
			}
		}
		if (count==0)
		{
			fprintf(stderr, "# Error: no particle with index %d\n", p);
			exit (1);
		}
		for (k=0; k<dim; k++)
		{
			mp[p*dim+k]/=count;
		}
	}



	double *dx2=(double *) malloc( total_t*sizeof(double) );
	POINTER_NULL(dx2);
	double *dx4=(double *) malloc( total_t*sizeof(double) );
	POINTER_NULL(dx4);

	double *NonGaussian=(double *) malloc (total_p*sizeof(double));
	POINTER_NULL(NonGaussian);

	double dx, dy, s2, s4;

	// only valid for dim=2
	for (p=0; p<total_p; p++)
	{
		count=0;
		for (t=0; t<total_t; t++)
		{
			i=index[t+p*total_t];
			if ( i != -1 )
			{
				pointer=ptid_ptr+i;
				/*
				for (k=0; k<dim; k++)
				{
					dx2[count] = (double) (*(pointer+k)) - mp[p*dim+k];
					dx2[count] *= dx2[count];
					dx4[count++] = dx2[count]*dx2[count];
				}*/
				dx=(double)(*pointer)-mp[p*dim];
				dy=(double)(*(pointer+1))-mp[p*dim+1];
				dx2[count]=dx*dx+dy*dy;
				dx4[count]=dx2[count]*dx2[count];
				++count;
			}
		}
		s2=pairwise(count, dx2)/count;
		s4=pairwise(count, dx4)/count;
		NonGaussian[p]=s4/s2/s2/2.0-1.0;
	}

	free(index);
	free(dx2);
	free(dx4);

	//=============================================================
	printf("# Particle number     :        %d\n", total_p);
	printf("# Frame number        :        %d\n", total_t);
	//============================================================
	//
	// According to practice, there are always a very few particles
	// having a very large DW factor.
	// So first the values will be adjusted.
	// Assuming they are Gaussian distributed.
	/*
	double aveDW = pairwise(total_p, DW)/total_p;
	double *DW2=(double *)malloc(total_p*sizeof(double));
	POINTER_NULL(DW2);
	for (p=0; p<total_p; p++)
	{
		DW2[p]=DW[p]-aveDW;
		DW2[p]*=DW2[p];
	}
	double sdDW = sqrt(pairwise(total_p, DW2)/total_p);
	*/

	const char subfix []=".VG";
	char *dwfile=getfilename(file, subfix);
	FILE * dwf=fopen(dwfile, "w");
	if ( dwf==NULL )
	{
		fprintf(stderr, "# Error: can not open file \"%s\"!\n", dwfile);
		exit (1);
	}
	
	fprintf(dwf, "# Vibrition nonGaussian parameter of \"%s\"\n", file);
	fprintf(dwf, "# Particle number     :        %d\n", total_p);
	fprintf(dwf, "# Frame number        :        %d\n", total_t);

	for (p=0; p<total_p; p++)
	{
		for (k=0; k<dim; k++)
			fprintf(dwf, "%f ", mp[p*dim+k]);
		fprintf(dwf, "%f\n", NonGaussian[p]);

		/*
		if ( DW[p] < aveDW-3*sdDW )
			fprintf(dwf, "%f\n", aveDW-3*sdDW);
		else if ( DW[p] > aveDW+3*sdDW )
			fprintf(dwf, "%f\n", aveDW+3*sdDW);
		else
			fprintf(dwf, "%f\n", DW[p]);
			*/
	}

	fclose(dwf);

	free(mp);
	free(NonGaussian);
	printf("# Vibrition nonGaussian parameter done.\n");

	/*=============================
	 * plot by gnuplot
	 *=============================*/

	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	char NULLsubfix [2];
	NULLsubfix[0]='\0';
	char *NULLfilename = getfilename(file, NULLsubfix);

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set pm3d map\n"
				   "set palette rgbformulae 33,13,10\n"
				   "set size ratio -1\n"
				   "set autoscale\n"
				   "set title \"%s vibrition nonGaussian parameter\"\n"
				   "#set xlabel 'x'\n"
				   "#set ylabel 'y'\n", NULLfilename);

	fprintf(pipe, "splot '%s' using 1:2:3 with points palette pt 7\n", dwfile);
	fprintf(pipe, "replot\n");

	fclose(pipe);

	free(dwfile);
	free(NULLfilename);

	delete [] size;
}
