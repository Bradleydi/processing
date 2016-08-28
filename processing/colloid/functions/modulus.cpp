#include "elasticity.h"

#include <cstdio>
#include <cstdlib>

void phelp ();

int main (int argc, char *argv[])
{
	if (argc==5)
		el_modulus (argv[1], atoi(argv[2]), atof(argv[3]), atoi(argv[4]));
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tmodulus gdfdata b a0 dim\n");
	printf("Parameters:\n");
	printf("\tgdfdata\n");
	printf("\t\tpretracked data gdf file\n");
	printf("\tb\n");
	printf("\t\tcoarse grain box number bxb\n");
	printf("\ta0\n");
	printf("\t\tlattice constant\n");
	printf("\tdim\n");
	printf("\t\tdimension\n");
	exit (0);
}
