#include "colloid_base.h"
#include "io.h"
#include "plot.h"

using namespace std;

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv [])
{
	int i, j=0, k=0;
	int twindow=10, tstep=1, framerate=2;
	char *str_frameindex=NULL;
	bool Tracked=false;

	if (argc==1) phelp();

	for (i=2; i<argc; i++)
	{
		if ( argv[i][0]=='-' ) // control option
		{
			if ( argv[i][1]=='w') // timewindow
				twindow=atoi(argv[++i]);
			else if ( argv[i][1]=='s') // tstep
				tstep=atoi(argv[++i]);
			else if ( argv[i][1]=='f') // framerate
				framerate=atoi(argv[++i]);
			else if ( argv[i][1]=='t') // tracked
				Tracked=true;
		}
		else
		{
			while ( argv[i][j]!='\0' )
				if ( argv[i][j++] == ':' )	// check whether is frameindex
					++k;

			if ( k == 2 || k == 1 )
			{
				str_frameindex=argv[i];
				k=3; // to prevent another frameindex argument
			}
			else
				phelp();
		}
	}
	if (str_frameindex==NULL) 
		traj_plot(argv[1], twindow, tstep, ":", framerate, Tracked);
	else
		traj_plot(argv[1], twindow, tstep, str_frameindex, framerate, Tracked);
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\ttraj_plot gdf_data_file [-w twindow] [-s tstep] [frameindex]"
			" [Tracked] [-f framerate]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\tdata file\n");
	printf("\t-w twindow\n"
		   "\t\ttime window, default = 10\n");
	printf("\t-s tstep\n"
		   "\t\ttime interval between two showed frames, default = 1\n");
	printf("\tframeindex\n"
		   "\t\tframeindex formatted string, to specify\n"
		   "\t\twhich frame should be plot, default =\":\"\n");
	printf("string 'frameindex' format:\n"
		   "\tN\n"
		   "\t\tset chosen frame as N (only one frame)\n"
		   "\tN1:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tThe increment is set to 1\n"
		   "\t\tA ':' means choosing all frames.\n"
		   "\tN1:c:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tc: increment, if not specified, set to 1\n");
	printf("\tTracked\n");
	printf("\t\t-t tracked, defalt = false\n");
	printf("\t-f framerate\n");
	printf("\t\t-f number, defalt = 2\n");
	exit (0);
}
