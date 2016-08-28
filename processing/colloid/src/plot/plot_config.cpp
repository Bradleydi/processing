#include "plot.h"

#include "colloid_base.h"
#include "miscellaneous.h"

#include<cstdio>
#include<cstdlib>
#include<cmath>

using namespace std;


/* using pipes and temp files */
void plot_config(colloid_base &p, const int n, bool tracked)
{
	int i;

	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	int ti=size[0]-1;
	if (tracked)
		--ti;

	float *pointer=p_ptr+ti;
	float t1=*pointer, t2=t1;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if ( t1 > *pointer )
			t1 = *pointer;
		else if ( t2 < *pointer )
			t2 = *pointer;
	}
	int total_t=(int)(t2-t1)+1;
	const int tmin=(int)t1;
	const int tmax=(int)t2;

	if (n<tmin)
	{
		fprintf(stderr, "frame no. should be greater than %d!\n", tmin);
		exit (0);
	}
	else if ( n >= tmax+1 )
	{
		fprintf(stderr, "frame no. exceeds the total frame number %d!\n", 
				tmax+1);
		exit (0);
	}


	/*==================================================================
	 * generate data file
	 *==================================================================*/

	char *tmpfilename = tmpnam (NULL);
	printf("# %s\n", tmpfilename);
	FILE *infile=fopen(tmpfilename, "w");
	if (infile==NULL)
	{
		fprintf(stderr, "# Error: can not create a temporyary file!\n");
		exit (1);
	}
	
	fprintf(infile, "# xy position of frame no. \t %d\n", n);

	for (i=0; i<size[1]; i++)
	{
		if ( (int)(*(p_ptr+ti)) == n )
			fprintf(infile, "%6.6f\t%6.6f\n", *p_ptr, *(p_ptr+1));
		p_ptr+=size[0];
	}

	fclose(infile);
	
	
	/*==================================================================
	 * Using pipes to GNUPLOT
	 *==================================================================*/

	FILE *pipe = popen("gnuplot -persist","w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: pipe can not open!\n");
		exit (1);
	}

	fprintf(pipe, "reset\n");
	fprintf(pipe, "unset key\n");
	fprintf(pipe, "set size ratio -1\n");
	fprintf(pipe, "set xlabel 'x'\n");
	fprintf(pipe, "set ylabel 'y'\n");
	if ( vGNUPLOT > vNEWGNUPLOT )
		fprintf(pipe, "set pointsize %.2f\n", 
				35.0/sqrt((float)(size[1])/total_t));

	else
		fprintf(pipe, "set style line 7 pt 7 lt -1\n");


	fprintf(pipe, "plot '%s' using 1:2 with points ls 7\n", tmpfilename);

	fclose(pipe);

	if ( remove(tmpfilename) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed!\n", tmpfilename);
}


void plot_config2(colloid_base &p, const int n, 
		const char *file, bool tracked)
{
	int i;

	int *size=p.get_size();
	float *p_ptr=p.get_array_pointer();

	int ti=size[0]-1;
	if (tracked)
		--ti;


	/*==================================================================
	 * generate data file
	 *==================================================================*/

	const char subfix1 []=".dat";
	char *filename1=getfilename(file, subfix1);
	FILE *infile=fopen(filename1, "w");
	if (infile==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename1);
		exit (1);
	}
	
	fprintf(infile, "# xy position of frame no. \t %d\n", n);

	for (i=0; i<size[1]; i++)
	{
		if ( (int)(*(p_ptr+ti)) == n )
			fprintf(infile, "%6.6f\t%6.6f\n", *p_ptr, *(p_ptr+1));
		p_ptr+=size[0];
	}

	fclose(infile);
	
	
	/*==================================================================
	 * generate GNUPLOT script file
	 *==================================================================*/

	const char subfix3 []=".eps";
	char *epsfilename=getfilename(file, subfix3);

	const char subfix2 []=".plt";
	char *filename2=getfilename(file, subfix2);
	infile=fopen(filename2, "w");
	if (infile==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename2);
		exit (1);
	}

	fprintf(infile, "reset\n\n");
	fprintf(infile, "set terminal postscript eps color enhanced\n");
	fprintf(infile, "set output '%s'\n\n", epsfilename);
	fprintf(infile, "unset key\n\n");
	fprintf(infile, "set xlabel 'x'\n");
	fprintf(infile, "set ylabel 'y'\n\n");
	
	if ( vGNUPLOT > vNEWGNUPLOT )
		fprintf(infile, "set pointsize 0.7f\n");	
	else
		fprintf(infile, "set style line 7 pt 7 lt -1\n");

	//fprintf(infile, "#set tmarg 2.2\n");
	//fprintf(infile, "#set bmarg 2\n");
	//fprintf(infile, "#set lmarg 4\n");
	//fprintf(infile, "#set rmarg 2\n\n");
	
	//fprintf(infile, "#set size ratio 47/65 #0.723\n\n");
	fprintf(infile, "filename=\"%s\"\n\n", filename1);
	fprintf(infile, "plot filename using 1:2 with points ls 7\n\n");
	
	fprintf(infile, "quit\n");
	

	fclose(infile);
	free(filename1);
	free(filename2);
	free(epsfilename);

}
