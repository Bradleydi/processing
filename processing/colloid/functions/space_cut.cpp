#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char * argv[])
{
	if ( argc == 3 )
		space_cut(/*filename*/ argv[1], argv[2]);
	else
		phelp();

	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tspace_cut gdf_data_file region\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked/untracked data\n");
	printf("\tregion\n");
	printf("\t\tthe remaining region to set, \n"
		   "\t\tformated as \"xmin:xmax:ymin:ymax\"\n");
	exit (0);
}
