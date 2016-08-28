#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <iostream>
#include <cstdlib>

using namespace std;

void phelp ();

int main (int argc, char *argv[])
{
	//char filename []="simulation8000.gdf";
	char *file=getfilename(argv[1], ".gdf");
	colloid_base ptid;
	readgdf(ptid, file);
	free(file);

	if (argc==3)
		uur(ptid, argv[1], atof(argv[2]), true);
	else if (argc==4 && argv[3][0]=='-' && argv[3][1]=='r')
		uur(ptid, argv[1], atof(argv[2]), false);
	else
		phelp();

	return 0;
}

void phelp ()
{
	printf("\tCalculate the displacement correlation\n"
		   "Usage:\n"
		   "\tuur gdffile a0 [ Remove_Drift ]\n"
		   "Parameters:\n"
		   "\tgdffile\n"
		   "\t\tgdfdata, normalmode should be called before this program\n"
		   "\ta0\n"
		   "\t\tlattice constant, can be obtained from g(r)\n"
		   "\tRemove_Drift\n"
		   "\t\t-r       not remove drift\n"
		   "\t\tdefalt   remove drift\n");
	exit (1);
}
