#include "miscellaneous.h"

#include <stdio.h>
#include <stdlib.h>

void show_matrix(int N, double *A)
{
	char *tmpfilename = tmpnam (NULL);
	FILE *infile=fopen(tmpfilename, "wb");
	FILE_NULL(infile, tmpfilename);
	fwrite(A, sizeof(A[0]), N*N, infile);
	fclose(infile);

	FILE *pipe = popen("gnuplot -persist","w");
	if (pipe==NULL) {
		fprintf(stderr, "# Error: pipe can not open!\n");
		exit (1);
	}
	fprintf(pipe, "reset\n"
			      "unset key\n"
				  "set size ratio -1\n"
				  "set pm3d map\n"
				  "set xrange [0:%d]\n"
				  "set yrange [0:%d]\n"
				  "splot '%s' binary array=%dx%d format=\"%%double\" with image\n",
				  N, N, tmpfilename, N, N);
	fclose(pipe);
	if ( remove(tmpfilename) != 0 )
	{
		fprintf(stderr, "# Error: tmpfile '%s' cannot be removed!\n", tmpfilename);
		exit (1);
	}
}
