#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==1)
		phelp();
	colloid_base ptid;
	readgdf(ptid, argv[1]);

	if (argc==2)
	{
		normal_mode(ptid, argv[1], 1, true);
	}
	else if (argc==3)
	{
		if (argv[2][0]=='-' && argv[2][1]=='r')
			normal_mode(ptid, argv[1], false);
		else
			normal_mode(ptid, argv[1], atof(argv[2]), true);
	}
	else if (argc==4)
		normal_mode(ptid, argv[1], atof(argv[2]), false);
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tnormal_mode gdf_data_file [a0] [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\ta0\n");
	printf("\t\tpixel per unit length, usually chosen as lattice constant\n");
	printf("\t\tobtained by radial distribution function g(r)\n");
	printf("\tremove_drift\n");
	printf("\t\t-r remove the drift\n");
	exit (0);
}
