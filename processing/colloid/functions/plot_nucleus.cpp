#include "colloid_base.h"
#include "io.h"
#include "structure.h"

using namespace std;

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv [])
{
	if (argc==3)
		plot_nucleus(argv[1], argv[2], 2, false);
	else if (argc==4 && argv[3][0]=='-' && argv[3][1]=='t')
		plot_nucleus(argv[1], argv[2], 2, true);
	else if (argc==5 && argv[3][0]=='-' && argv[3][1]=='f')
		plot_nucleus(argv[1], argv[2], atoi(argv[4]), false);
	else if (argc==6 && argv[3][0]=='-' && argv[3][1]=='t'
			&& argv[4][0]=='-' && argv[4][1]=='f')
		plot_nucleus(argv[1], argv[2], atoi(argv[5]), true);
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tplot_nucleus gdf_data_file frameindex [Tracked] [-f framerate]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\tdata file\n");
	printf("\tframeindex\n"
		   "\t\tframeindex formatted string, to specify\n"
		   "\t\twhich frame should be plot\n");
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
	printf("\t\t-t tracked\n");
	printf("\tframerate\n");
	printf("\t\t-f number\n");
	exit (0);
}
