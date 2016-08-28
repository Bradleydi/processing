#include "structure.h"

#include <cstdio>
#include <cstdlib>


using namespace std;

void phelp();


int main (int argc, char *argv[])
{
	if (argc == 2)
		s2(argv[1], 2, false);
	else if (argc == 3 && argv[2][0]=='-' && argv[2][1]=='t' )
		s2(argv[1], 2, true);
	else if (argc == 3 && argv[2][0]=='-' && argv[2][1]=='d' 
			 && argv[2][2]=='3' )
		s2(argv[1], 3, false);
	else if (argc == 4 && argv[2][0]=='-' && argv[2][1]=='t' 
			 && argv[3][0]=='-' && argv[3][1]=='d' 
			 && argv[3][2]=='3' )
		s2(argv[1], 3, true);
	else if (argc == 4 && argv[2][0]=='-' && argv[2][1]=='d' 
			 && argv[2][2]=='3'
			 && argv[3][0]=='-' && argv[3][1]=='t' )
		s2(argv[1], 3, true);
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\ts2 gdf \n");
	printf("Parameters:\n");
	printf("\ttrk_ptid_gdf\n");
	printf("\t\ttracked data gdf file with prefix \"trk\"\n");
	printf("\tpt_gdf\n");
	printf("\t\tfilename of the resulted pt data.\n"); 
	printf("\t\tif not specified, filename will be the one\n");
	printf("\t\twithout the prefix \"trk\"\n");
	exit (0);
}
