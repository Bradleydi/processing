#include "statistics.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void phelp();

int main (int argc, char *argv[])
{
	//char filename []="simulation8000.gdf";
	if (argc==1)
		phelp();
	if (argc==2)
		pt_info(argv[1], /*dim=*/ 0, /*show_histogram=*/ false);
	else if (argc==3)
		if (strcmp(argv[2], "-PLOT")==0)
			pt_info(argv[1], /*dim=*/ 0, /*show_histogram=*/ true);
		else
			pt_info(argv[1], /*dim=*/ atoi(argv[2]), /*show_histogram=*/ false);
	else if (argc==4 ||  (strcmp(argv[3], "-PLOT")==0) )
		pt_info(argv[1], /*dim=*/ atoi(argv[2]), /*show_histogram=*/ true);
	else phelp();
	return 0;
}

void phelp()
{
	printf("\tPrint and plot the information of a given untracked data\n");
	printf("Usage:\n");
	printf("\tpt_info gdf_data_file [dim] [-PLOT]\n");
	printf("Parameters:\n");
	printf("\tdim\n");
	printf("\t\tspecify dimension of the data\n");
	printf("\t\tdefault value is column_number-2\n");
	printf("\t-PLOT\n");
	printf("\t\tplot the histogram of particles in each frame and\n"); 
	printf("\t\tsome other details\n");
	exit (0);
}
