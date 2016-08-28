#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==4)
	{
		softspot(argv[1], atoi(argv[2]), atoi(argv[3]));
	}
	else 
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tnormal_mode gdf_data_file [a0]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\ta0\n");
	printf("\t\tpixel per unit length, usually chosen as lattice constant\n");
	printf("\t\tobtained by radial distribution function g(r)\n");
	exit (0);
}
