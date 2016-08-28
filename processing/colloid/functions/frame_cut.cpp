#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char * argv[])
{
	if ( argc == 3 )
		frame_cut(/*filename*/ argv[1], atoi(argv[2]));
	else if ( argc == 4 && argv[3][0]=='-' && argv[3][1]=='t' )
		frame_cut(/*filename*/ argv[1], atoi(argv[2]), true);
	else
		phelp();

	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tframe_cut gdf_data_file maxt [Tracked]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked/untracked data\n");
	printf("\tmaxt\n");
	printf("\t\tmaximal frame number to set\n");
	printf("\tTracked\n");
	printf("\t\t-t tracked\n");
	exit (0);
}
