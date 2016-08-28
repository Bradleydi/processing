#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

void plot_K(const char *file)
{
	// read K
	int Kn;
	double *K;
	readK(Kn, &K, file);

	// read mean position
	int total_p;
	float *mp_ptr;
	readmp(total_p, &mp_ptr, file);
	
	if ( Kn != 2*total_p)
		pERROR("read files error!");

	int i, j;
	const int dim=2;

	float xmin=mp_ptr[0], xmax=mp_ptr[0];
	float ymin=mp_ptr[1], ymax=mp_ptr[1];
	float *pointer=mp_ptr+dim;
	for (i=1; i<total_p; i++)
	{
		if ( xmin > *pointer)
			xmin = *pointer;
		else if ( xmax < *pointer)
			xmax = *pointer;
		++pointer;

		if ( ymin > *pointer)
			ymin = *pointer;
		else if ( ymax < *pointer)
			ymax = *pointer;
		++pointer;
	}

	char *filename=getfilename(file, ".Kf");

	FILE *Kf=fopen(filename, "w");
	FILE_NULL(Kf, filename);
	/*
	double threshold=0.0, keff;
	for (i=0; i<total_p*2; i++)
		threshold += K[i*Kn+i];
	threshold /= total_p*10; // 1/3 of average Kxxand Kyy
	for (i=0; i<total_p; i++)
	{
		for (j=i+1; j<total_p; j++) // not itself
		{
			keff = fabs(K[2*i*Kn+2*j] + K[2*i*Kn+Kn+2*j+1]);
			if ( keff >= threshold )
			{
				fprintf(Kf, "%f %f %f\n%f %f %f\n\n",
						mp_ptr[i*dim], mp_ptr[i*dim+1], keff,
						mp_ptr[j*dim], mp_ptr[j*dim+1], keff);
			}
		}
	}*/
	float dx, dy, d2, dcut=1.5, dcut2=dcut*dcut;
	double keff;
	for (i=0; i<total_p; i++)
	{
		for (j=i+1; j<total_p; j++) // not itself
		{
			dx=(mp_ptr[j*dim]-mp_ptr[i*dim]);
			dy=(mp_ptr[j*dim+1]-mp_ptr[i*dim+1]);
			d2 = dx*dx+dy*dy;
			if ( d2 < dcut2 )
			{
				keff = fabs(K[2*i*Kn+2*j] + K[2*i*Kn+Kn+2*j+1]);
				fprintf(Kf, "%f %f %f\n%f %f %f\n\n",
						mp_ptr[i*dim], mp_ptr[i*dim+1], keff,
						mp_ptr[j*dim], mp_ptr[j*dim+1], keff);
			}
		}
	}
	Fclose(Kf, filename);
	free(mp_ptr);
	free(K);
	
	char *eps = getfilename(file, "_K.eps");
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set terminal postscript eps color enhanced\n"
				   "set output \"%s\"\n"
				   "set size ratio -1\n"
				   //"set autoscale\n"
				   "set xrange [%f:%f]\n"
				   "set yrange [%f:%f]\n"
				   "unset key\n"
				   "unset xtics\n"
				   "unset ytics\n"
				   "set border 0\n"
				   "set palette rgbformulae 33,13,10\n"
				   "#set xlabel 'x/a_0'\n"
				   "#set ylabel 'y/a_0'\n", eps, 
				   xmin-1.0, xmax+1.0, ymin-1.0, ymax+1.0
				   );


	fprintf(pipe, "plot '%s' using 1:2:3 with lines lc palette lw 3\n", filename);

	fclose(pipe);
	free(filename);
}
