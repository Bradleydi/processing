#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv [] )
{
	if ( argc == 2 )
		traj_length(argv[1]);
	else phelp();
	return 0;
}

void phelp()
{
	printf("\tplot the trajectories length in space, to see where the track\n"
		   "\tis good. The length is renormalized with respect to the total\n"
		   "\tframe number.\n"
		   "\tTo get statistics, use 'ptid_info gdffile -PLOT'\n");
	printf("Usage:\n");
	printf("\ttraj_length gdffile\n");
	printf("Parameters:\n");
	printf("\tgdffile\n");
	printf("\t\ttracked data file\n");
	exit (0);
}
