#include "colloid_base.h"
#include "io.h"
#include "plot.h"

using namespace std;

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int agrc, char *argv [])
{
	if (agrc==1)
		phelp();
	colloid_base p;
	readgdf(p, argv[1]);

	if (argv[1][0]=='t' &&
		argv[1][1]=='r' &&
		argv[1][2]=='k')
		plot_config(p, atoi(argv[2]), true);
	else
		plot_config(p, atoi(argv[2]), false);
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tplot_config gdf_data_file frame_index\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data should be prefixed by 'trk'\n");
	printf("\tframe_index\n");
	printf("\t\tindex of frame to be showed\n");
	exit (0);
}
