#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <iostream>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc!=3)
		phelp();
	//char filename []="simulation8000.gdf";
	char *file=getfilename(argv[1], ".gdf");
	colloid_base ptid;
	readgdf(ptid, file);
	free(file);

	if (argc==3)
		mode_hist(ptid, atoi(argv[2]), argv[1], true);
	else
		phelp();

	return 0;
}

void phelp ()
{
	printf("\tCalculate the amplitude distribution of a given normal mode\n"
		   "Usage:\n"
		   "\tmode_hist gdffile mode_id\n"
		   "Parameters:\n"
		   "\tgdffile\n"
		   "\t\tgdfdata, normalmode should be called before this program\n"
		   "\tmode_id\n"
		   "\t\tID of the mode to plot the mode histogram\n");
	exit (1);
}
