#include "structure.h"
#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"

#include <fftw3.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>


void sk(const char *file, bool sk2D, bool tracked)
{
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

	float *pointer=p_ptr+ti;
	float t1=*pointer, t2=t1;
	for (i=0; i<size[1]; i++)
	{
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
		pointer+=size[0];
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin=(int)t1;
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
		pERROR("particle fluctuation is too large!");

	int *index_now=(int *)calloc(total_t, sizeof(int));
	POINTER_NULL(index_now);

	float **index=(float **)malloc(size[1]*sizeof(float *));
	POINTER_NULL(index);
	pointer=p_ptr;
	float xmin=1.0e6, xmax=-1.0e6, ymin=xmin, ymax=xmax;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;

		if ( xmin > *pointer )
			xmin = *pointer;
		else if ( xmax < *pointer )
			xmax = *pointer;

		++pointer;
		if ( ymin > *pointer )
			ymin = *pointer;
		else if ( ymax < *pointer )
			ymax = *pointer;


		pointer+=size[0]-1;
	}
	free(index_now);

	xmin=floor(xmin);
	ymin=floor(ymin);
	
	// cut off offset xmin, ymin
	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		*(pointer++) -= xmin;
		*pointer -=  ymin;
		pointer+=size[0]-1;
	}

	const int NX=(int)(ceil(xmax-xmin));
	const int NY=(int)(ceil(ymax-ymin));
	const int NN=NX*NY;
	const int NK=(NX/2)+1;
	const int NNK=NY*NK;
	const double PI = acos(-1.0);
	const double Deltakx = 2.0*PI/NX; // Delta kx = 2pi/N/Delta x
	const double Deltaky = 2.0*PI/NY;
	// so the unit is 1pixel^-1.
	
	//=================================================================
	
	// begin to calculate sk for each frame
	double *image=Malloc(double, NN);
	POINTER_NULL(image);
	fftw_complex* dft=(fftw_complex *)fftw_malloc(NNK*sizeof(fftw_complex));
	POINTER_NULL(dft);
	double *sk=Calloc(double, NNK);
	POINTER_NULL(sk);
	// Create plans
	//fftw_plan plan=fftw_plan_dft_r2c_2d(NY, NX, image, dft, FFTW_ESTIMATE);
	fftw_plan plan=fftw_plan_dft_r2c_2d(NY, NX, image, dft, FFTW_MEASURE);

	double Re, Im;
	int xi, yi, j;
	float xf, yf, xfyf;
	// begin the big loop
	for (t=0; t<total_t; t++)
	{
		for (i=0; i<NN; i++)
			image[i]=0.0;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			// index_x is column, index_y is row
			// so row-major means index_y*Ncol+index_x
			xf=*pointer;
			xi=(int)xf;
			xf -= xi;
			yf = *(pointer+1);
			yi=(int)yf;
			yf -= yi;
			xfyf=xf*yf;
			//image[((int)(*(pointer+1)))*NX+(int)(*pointer)]=1.0;
			j=yi*NX+xi;

			image[j]      +=  1.0-xf-yf-xfyf;
			image[j+1]    +=  xf-xfyf;
			image[j+NX]   +=  yf-xfyf;
			image[j+NX+1] +=  xfyf;
		}

		// Compute forward DFT
		fftw_execute(plan);

		// compute sk of frame t
		// sk = rho_k*rho_(-k)=rho_k*rho_k^*
		for (i=0; i<NNK; i++)
		{
			Re = dft[i][0];
			//Im = dft[i][0]; // a bug here
			Im = dft[i][1]; 
			sk[i] += Re*Re + Im*Im;
		}
	}
	double skmax=-1.0;
	for (i=0; i<NNK; i++)
	{
		sk[i] /= total_t*NN;
		if ( skmax < sk[i] )
			skmax = sk[i];
	}

	// Free memory
	fftw_destroy_plan(plan);
	fftw_free(dft);
	free(image);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	free(index);

	//==================================================================
	// put each element to the right place.
	// Since in the fftw format, the range is from [0:NY-1]x[0:NX/2],
	// and what we need is [-NY/2:NY/2]x[-NX/2:NX/2]
	// so we need to construct it from the original data
	//==================================================================

	/*
	char *tmpfile=getfilename(file, ".SK");
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);
	
	if ( (int)fwrite(sk, sizeof(double), NNK, infile) != NNK )
	{
		fprintf(stderr, "# Error: write to file '%s' failed!\n", tmpfile);
		exit (1);
	}
	fclose(infile);

	FILE *pipe=popen("gnuplot -p", "w");
	if (pipe==NULL)
		pERROR("Gnuplot pipe can not open.");
	fprintf(pipe, "reset\n"
				  "unset key\n"
				  "set size ratio -1\n"
				  "set pm3d map\n"
				  "set palette gray\n"
				  "unset colorbox\n"
				  "set cbrange [0:]\n"
				  "splot '%s' binary array=%dx%d format=\"%%double\" "
				  "with image palette\n"
				  "exit\n", 
				  tmpfile, NK, NY);
	fclose(pipe);
	printf("reset\n"
				  "unset key\n"
				  "set size ratio -1\n"
				  "set pm3d map\n"
				  "set palette gray\n"
				  "unset colorbox\n"
				  "set cbrange [0:]\n"
				  "splot '%s' binary array=%dx%d format=\"%%double\" "
				  "with image palette\n"
				  "exit\n", 
				  tmpfile, NK, NY);
	*/
	int i1=0, i2=0;
	const int NX_2=NX/2;
	const int NY_2=NY/2;
	int NXY=(NY_2*2+1)*(NX_2*2+1);
	double *SK=Malloc(double, NXY);
	POINTER_NULL(SK);
	--NXY;
	// first quadrant
	// in SK (i,j)  0<=i<= NY_2,  NX_2<=j<=NX/2*2
	// in sk (NY_2-i, j-NX_2)
	// so it's in sk region
	for (i=0; i<=NY_2; i++)
	{
		// i1=i*(NX_2*2+1)+j; j=NX_2
		i1=i*(NX_2*2+1)+NX_2;
		//i2=(NY_2-i)*NK+j-NX_2; j=NX_2
		i2=(NY_2-i)*NK;
		for (j=NX_2; j<=NX/2*2; j++)
		{
			SK[i1]=sk[i2++];
			SK[NXY-i1]=SK[i1];
			++i1;
		}
	}
	// fourth quadrant
	// in SK (i,j)  NY_2<i<= NY_2*2,  NX_2<=j<=NX/2*2
	// in sk (NY_2-i, j-NX_2)
	// so it's not in sk region
	// and it should be shift to (NY+NY_2-i, j-NX_2)
	// hence the first index range is [NY-NY_2, NY)
	for (i=NY_2+1; i<=NY_2*2; i++)
	{
		// i1=i*(NX_2*2+1)+j; j=NX_2
		i1=i*(NX_2*2+1)+NX_2;
		//i2=(NY+NY_2-i)*NK+j-NX_2; j=NX_2
		i2=(NY+NY_2-i)*NK;
		for (j=NX_2; j<=NX/2*2; j++)
		{
			SK[i1]=sk[i2++];
			SK[NXY-i1]=SK[i1];
			++i1;
		}
	}
	++NXY;
	free(sk);
	printf("# saving data...\n");

	if ( sk2D )
	{
		// plot S(k)  // forget to change y from [0, NY-1] to [-NY/2+1, NY/2-1]
		char *tmpfile=getfilename(file, ".SK");
		FILE *infile=fopen(tmpfile, "w");
		FILE_NULL(infile, tmpfile);
		
		if ( (int)fwrite(SK, sizeof(double), NXY, infile) != NXY )
		{
			fprintf(stderr, "# Error: write to file '%s' failed!\n", tmpfile);
			exit (1);
		}
		fclose(infile);

		char *epsfile=getfilename(file, "_SK.eps");
		FILE *pipe=popen("gnuplot -p", "w");
		if (pipe==NULL)
			pERROR("Gnuplot pipe can not open.");
		fprintf(pipe, "reset\n"
					  "set terminal postscript eps color enhanced\n"
					  "set output '%s'\n"
					  "unset key\n"
					  "unset tics\n"
					  "set size ratio -1\n"
					  "set xrange [-pi:pi]\n"
					  "set yrange [-pi:pi]\n"
					  "set pm3d map\n"
					  "set palette gray\n"
					  "unset colorbox\n"
					  "set cbrange [0:%f]\n"
					  "set label 1 '+' at 0,0 center front nopoint textcolor lt 1\n"
					  "splot '%s' binary array=%dx%d format=\"%%double\" "
					  "dx=%f dy=%f origin=(-pi, -pi, 0) with image palette\n"
					  "exit\n", epsfile, 
					  //(NX_2-NX)*Deltakx, NX_2*Deltakx, // they are PI
					  //(NY_2-NY)*Deltaky, NY_2*Deltaky,
					  skmax/20.0, // to make the image brighter
					  tmpfile, 
					  NX_2*2+1, NY_2*2+1, // array axb, row first: a,column no.
					  Deltakx, Deltaky);
		fclose(pipe);
		free(epsfile);
		//if ( remove(tmpfile) != 0 )
		//	fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);
	}
	else
	{
		// get angle-averaged S(k)
		int nbin, window=3;
		//if ( NX > NY ) nbin=(int)sqrt(NY_2);
		//else nbin=(int)sqrt(NX_2);
		// since the resolution can be Deltakx or Deltaky, so we use
		double binsize, k;//=PI/nbin;
		if ( NX > NY ) 
		{ binsize=Deltaky*window; nbin=NY_2/window+1; }
		else 
		{ binsize=Deltakx*window; nbin=NX_2/window+1; }

		unsigned int *count=Calloc(unsigned int, nbin);
		POINTER_NULL(count);
		double *skr=Calloc(double, nbin); POINTER_NULL(skr);
		i2=0;
		for (i=NY_2; i>=0; i--)
		{
			for (j=-NX_2; j<=NX_2; j++)
			{
				k=sqrt(i*i*Deltaky*Deltaky+j*j*Deltakx*Deltakx);
				if ( k < PI)
				{
					i1=(int)(k/binsize);
					skr[i1] += SK[i2++];
					if (i==0)	++count[i1];
					else	count[i1] +=2;
				}
			}
		}


		char *skfile=getfilename(file, ".sk");
		FILE *skf=fopen(skfile, "w");
		FILE_NULL(skf, skfile);
		
		fprintf(skf, "# sk(r)\n"
					 "# r unit: 1 pixel^-1\n");
		for (i=0; i<nbin; i++)
		{
			if ( count[i] != 0 )
			{
				skr[i] /= count[i]/2.0;
				fprintf(skf, "%f %f %u\n", (i+0.5)*binsize, skr[i], count[i]);
			}
			else
				fprintf(skf, "%f 0.0 0\n", (i+0.5)*binsize);
		}
		fclose(skf);
		free(count);
		free(skr);
	}

	free(SK);
	delete [] size;
}
