#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void showposition(const char *file, const int factor, bool Remove_Drift)
{
	char *gdffile = getfilename(file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);

	int *size = p.get_size();
	float *p_ptr = p.get_array_pointer();

	if (Remove_Drift)
		remove_drift(p, 2);

	float *pointer = p_ptr;
	int i;
	float xmin=*pointer, xmax=xmin, ymin=*(pointer+1), ymax=ymin;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];

		if ( xmin > *pointer )
			xmin = *pointer;
		else if ( xmax < *pointer )
			xmax = *pointer;

		if ( ymin > *(pointer+1) )
			ymin = *(pointer+1);
		else if ( ymax < *(pointer+1) )
			ymax = *(pointer+1);
	}

	const int ixmax=(int)(xmax+1.5);
	const int ixmin=(int)(xmin-0.5);
	const int iymax=(int)(ymax+1.5);
	const int iymin=(int)(ymin-0.5);

	const int col = (ixmax-ixmin)*factor;
	const int row = (iymax-iymin)*factor;

	unsigned char *image=(unsigned char *)calloc(col*row, sizeof(unsigned char));
	POINTER_NULL(image);

	const int offset = ixmin+iymin*col;

	pointer=p_ptr;
	for (i=0; i<size[1]; i++)
	{
		*(image+((int)(*pointer+0.5)+((int)(*(pointer+1)+0.5))*col-offset)*factor)=1;
		pointer+=size[0];
	}

	char *tmpfile = tmpnam (NULL);
	printf("# %s\n", tmpfile);
	FILE *infile=fopen(tmpfile, "w");
	FILE_NULL(infile, tmpfile);

	fprintf(infile, "# xy position of %s\n", gdffile);
	int x,y;
	for (x=0; x<col; x++)
	{
		for (y=0; y<row; y++)
		{
			if (image[x+y*col])
				fprintf(infile, "%d %d 1\n", x+ixmin, y+iymin);
			else
				fprintf(infile, "%d %d 0\n", x+ixmin, y+iymin);
		}
	}
	fclose(infile);
	free(image);

	free(gdffile);

	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"gnuplot\"!\n");
		exit (1);
	}

	fprintf(pipe, "reset\n"
				  "unset key\n"
				  "set size ratio -1\n"
				  "set tmargin at screen 0.95\n"
				  "unset key\n"
				  "set pm3d map\n"
				  "set palette gray\n"
				  "unset colorbox\n"
				  "set xrange [%d:%d]\n"
				  "set yrange [%d:%d]\n"
				  "splot \"%s\" using 1:2:3 with image palette\n"
				  "exit\n", ixmin, ixmax, iymin, iymax, tmpfile);
	fclose(pipe);

	if ( remove(tmpfile) != 0 )
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed.\n", tmpfile);
}
