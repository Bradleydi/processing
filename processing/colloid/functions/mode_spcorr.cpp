#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==3)
	{
		if ( argv[2][0] != '-' )
			mode_spcorr(atoi(argv[2]), argv[1]);
		else if ( argv[2][0] == '-' && argv[2][1] == 'a' )
			mode_spcorr_all(argv[1], -1.); // all modes
		else
			phelp();

	}
	else if ( argc == 4 && argv[2][0] == '-' && argv[2][1] == 'a' )
		mode_spcorr_all(argv[1], atof(argv[3])); // parts of modes
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tmode_spcorr gdf_data_file n [ -a [ omegaRange ] ]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\tn\n");
	printf("\t\tchoose to plot the n'th mode\n");
	printf("\t\twe can also replace this argument n by '-a'\n"
		   "\t\t to get correlations of all modes\n"
		   "\t\ta omegaRange can be given to specify the upper bound of\n"
		   "\t\tthe plotted frequency axis\n");
	exit (0);
}
