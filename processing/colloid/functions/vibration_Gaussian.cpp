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
		vibration_Gaussian(argv[1], true);
	else if (argc==3 && argv[2][0]=='-' && argv[2][1]=='r')
		vibration_Gaussian(argv[1], false);
	else
		phelp();
	
	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tvibration_Gaussian gdf_data_file [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\tremove_drift\n");
	printf("\t\t-r not remove the drift\n");
	exit (0);
}
