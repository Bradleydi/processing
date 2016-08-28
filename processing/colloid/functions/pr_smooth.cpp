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
		participation_ratio_smooth(argv[1], /*h=*/0.0);
	else if (argc==3)
		participation_ratio_smooth(argv[1], /*h=*/(double)atof(argv[2]));
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Kernel density estimation of DOS\n");
	printf("Usage:\n");
	printf("\tDOS_density gdf_data_file [h] [nbin]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\th\n");
	printf("\t\tbandwidth\n");
	printf("\tnbin\n");
	printf("\t\tnumber of bins\n");
	exit (0);
}
