#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv [] )
{
	if ( argc == 2 )
		IR_SSD(argv[1], /* dim = */ 2, /*window= */ 0);
	else if ( argc == 3 )
		IR_SSD(argv[1], /* dim = */ 2, /*window= */ atoi(argv[2]));
	else phelp();
	return 0;
}

void phelp()
{
	printf("\tTo check whether there is particle irreversible rearrangement (IR)\n");
	printf("Usage:\n");
	printf("\tIR gdffile [dim] [window]\n");
	printf("Parameters:\n");
	printf("\tgdffile\n");
	printf("\t\ttracked data file\n");
	printf("\twindow\n"
		   "\t\twindow for steady state detection (SSD)\n"
		   "\t\tdefalut is maxframe/10\n"
		   "\t\tusually this is fairly good\n");
	exit (0);
}
