#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv[])
{
	if (argc==5) 
		e_softspot(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), 0.9);
	else if (argc==6)
		e_softspot(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]),
				atof(argv[5]));
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("\tPlot soft spots of ellipsoids crystals/glasses\n"
		   "\tTo calculate normal modes, use 'e_normalmode'\n"
		   "Usage:\n"
		   "\te_softspot gdfdata a b a0 [threshold]\n"
		   "Parameters:\n"
		   "\tgdfdata\n"
		   "\t\ttracked gdf data\n"
		   "\ta b\n"
		   "\t\tsemi-long and semi-short axis length\n"
		   "\t\tunit: pixels\n"
		   "\ta0\n"
		   "\t\tpixel per unit length, for ellipse, use b as new unit\n"
		   "\tthreshold\n"
		   "\t\tthreshold of the contribution of displacement field.\n"
		   "\t\tDefault is 0.9\n"
		   "Output:\n"
		   "\ttwo eps figures with subfix '_SPT' and '_SPR'\n");
	exit (0);
}
