#include "dynamics.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	//void DWfactor(const char *file, double a0, const int & dim, bool Remove_drift)
	if (argc==1)
		phelp();
	else if (argc==2)
		DWfactor(argv[1], -1, 2, true);
	else if (argc==3)
	{
		if (argv[2][0]=='-' && argv[2][1]=='r')
			DWfactor(argv[1], -1, 2, false);
		else
			DWfactor(argv[1], atof(argv[2]), 2, true);
	}
	else if (argc==4)
		DWfactor(argv[1], atof(argv[2]), 2, false);
	else
		phelp();
	
	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tDWfactor gdf_data_file [a0] [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\ta0\n");
	printf("\t\tpixel per unit length, usually chosen as lattice constant\n");
	printf("\t\tobtained by radial distribution function g(r)\n");
	printf("\tremove_drift\n");
	printf("\t\t-r not remove the drift\n");
	exit (0);
}
