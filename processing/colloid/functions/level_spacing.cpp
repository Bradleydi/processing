#include "colloid_base.h"
#include "normalmode.h"
#include "statistics.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==2)
		level_spacing(argv[1], 0, KDE);
	else if (argc==3)
	{
		if ( argv[2][0]=='-' && argv[2][1]=='l' )
			level_spacing(argv[1], 0, LLR);
		else
			level_spacing(argv[1], atof(argv[2]), KDE);
	}
	else if (argc==4 && argv[3][0]=='-' && argv[3][1]=='l' )
		level_spacing(argv[1], atof(argv[2]), LLR);
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n"
		   "\tlevel spacing distribution of normal modes\n");
	printf("Usage:\n");
	printf("\tlevel_spacing gdf_data_file\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	exit (0);
}
