#include "colloid_base.h"
#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==1)
		phelp();

	if (argc==2)
		Kr(argv[1]);
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tKr gdf_data_file\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	exit (0);
}
