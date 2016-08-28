#include "data_preprocess.h"
#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

void ptid_drift(const char *filename, int dim)
{
	colloid_base ptid;
	char *gdffile=getfilename(filename,".gdf");
	readgdf(ptid, gdffile);

	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();
	const int ti = size[0]-2;
	const int ii = size[0]-1;

	/* test tracked data */
	if ( dim == 0 )
		dim = ti;
	else
	{
		if ( dim < 0 )
			pERROR("given dimension less than 0!");
		if ( dim+2 > size[0] )
			pERROR("not a tracked data!");
	}
	
	if (dim!=2)
		pERROR("not implimented yet!");

	/* get total_p */
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	int *frame_ptcl = (int *) calloc ( total_p, sizeof(int) );
	POINTER_NULL(frame_ptcl);

	/* get total_t */
	int i;
	float *pointer=ptid_ptr+ti;
	++frame_ptcl[(int)(*(pointer+1))];
	float tt=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( tt < *pointer )
			tt = *pointer;
		++frame_ptcl[(int)(*(pointer+1))];
	}

	const int total_t=(int)tt+1;
	
	/*=============================
	 * find out particles existing in all frames
	 *=============================*/
	unsigned char *good=(unsigned char *)calloc(total_p, sizeof(unsigned char));
	POINTER_NULL(good);
	int all_p=0;
	for (i=0; i<total_p; i++)
	{
		if ( frame_ptcl[i] == total_t )
		{
			++good[i];
			++all_p;
		}
	}
	free(frame_ptcl);
	printf("# number of particles                       : %d\n", total_p);
	printf("# number of particles existing in all frames: %d\n", all_p);
	if ( all_p < 2*total_p/3)
		pERROR("bad tracking, too many broken trajectories!");

	double *xdrift=(double *)calloc(total_t, sizeof(double));
	double *ydrift=(double *)calloc(total_t, sizeof(double));
	POINTER_NULL(xdrift);
	POINTER_NULL(ydrift);

	pointer=ptid_ptr;
	int t;
	for (i=0; i<size[1]; i++)
	{
		if ( good [ (int)(*(pointer+ii)) ] )
		{
			t=(int)(*(pointer+ti));
			xdrift[t] += (double)(*pointer);
			ydrift[t] += (double)(*(pointer+1));
		}
		pointer += size[0];
	}
	free(good);
	delete [] size;
	
	double xc0=(*xdrift)/all_p;
	double yc0=(*ydrift)/all_p;
	
	char subfix [] = ".tmp";
	char *tmp = getfilename(filename, subfix);
	if ( fopen(tmp, "r") != NULL )
	{
		fprintf (stderr, "# Error: temperal file \"%s\" exists!\n", tmp); 
		exit (1);
	}

	FILE *tmpf=fopen(tmp, "w");
	if (tmpf==NULL)
	{
		fprintf (stderr, "# Error: cannot open temperal file \"%s\"!\n", tmp);
		exit (1);
	}

	for (t=0; t<total_t; t++)
	{
		xdrift[t]=xdrift[t]/all_p-xc0;
		ydrift[t]=ydrift[t]/all_p-yc0;
		fprintf(tmpf, "%d %f %f\n", t, xdrift[t], ydrift[t]);
	}
	fclose(tmpf);
	
	double ax, bx, ay, by;

	linearfit(total_t, xdrift, ax, bx);
	linearfit(total_t, ydrift, ay, by);
	free(xdrift);
	free(ydrift);


	//printf("# plot...\n");
	/* plot */
	//=========================================================
	// plot by GNUplot, using pipe

	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe for gnuplot!\n");
		exit (1);
	}
	
	fprintf(pipe,  "reset\n"
				   "set autoscale\n"
				   "set title \"drift of %s\\n"
				   "total particle = %d\\n"
				   "total frame    = %d\"\n"
				   "set xlabel 't'\n"
				   "set ylabel 'drift'\n", 
				   filename, total_p, total_t);

	fprintf(pipe, "plot \"%s\" using 1:2 title \"x\",\\\n"
				       "\"%s\" using 1:3 title \"y\",\\\n"
					   "%f*x+%f title\"x-fit\",\\\n"
					   "%f*x+%f title\"x-fit\"\n"
					   , tmp, tmp, ax, bx, ay, by);

	fclose(pipe);

	if( remove( tmp ) != 0 )
	{
		fprintf(stderr, "# Error: deleting file \"%s\" failed!", tmp);
		exit (1);
	}
	free(tmp);
}


/*
	//=========================================================
	// plot by GNUplot, using pipe
	char epsTLsubfix [17];
	sprintf(epsTLsubfix, "_mode_%d.eps", n);
	char *eps = getfilename(file, epsTLsubfix);


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
				   "set terminal postscript eps color enhanced\n"
				   "set output \"%s\"\n"
				   "set size ratio -1\n"
				   "set autoscale\n"
				   "set title \"%s\\n"
				   " %d'th normal mode\\nomega=%f\"\n"
				   "set xlabel 'x/a_0'\n"
				   "set ylabel 'y/a_0'\n", eps, NULLfilename, n, E);


	fprintf(pipe, "plot '%s' using 1:2:3:4 with vectors\n", modefilename);

	fclose(pipe);

	free(modefilename);
	free(NULLfilename);
	free(eps);
}
}*/
