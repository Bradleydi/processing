#include "structure.h"

#include <cstdio>
#include <cstdlib>


using namespace std;

void phelp();


int main (int argc, char *argv[])
{
	if (argc == 2)
		triangle(argv[1], 0.1, false);
	else if (argc == 3 && argv[2][0]=='-' && argv[2][0]=='t' )
		triangle(argv[1], 0.1, true);
	else if (argc == 4 && argv[2][0]=='-' && argv[2][1]=='c' )
		triangle_config(argv[1], atof(argv[3]), false);
	else if (argc == 5 && argv[2][0]=='-' && argv[2][1]=='c' 
			 && argv[4][0]=='-' && argv[4][1]=='t' )
		triangle_config(argv[1], atof(argv[3]), true);
	else if (argc == 4 && argv[2][0]=='-' && argv[2][1]=='l' )
		triangle_cluster(argv[1], atof(argv[3]), false);
	else if (argc == 5 && argv[2][0]=='-' && argv[2][1]=='l' 
			 && argv[4][0]=='-' && argv[4][1]=='t' )
		triangle_cluster(argv[1], atof(argv[3]), true);
	else if (argc == 4 && argv[2][0]=='-' && argv[2][1]=='p' )
		triangle_cluster_plot(argv[1], atoi(argv[3]), false);
	else if (argc == 5 && argv[2][0]=='-' && argv[2][1]=='p' 
			 && argv[4][0]=='-' && argv[4][1]=='t' )
		triangle_cluster_plot(argv[1], atoi(argv[3]), true);
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tptid2pt trk_ptid_gdf [pt_gdf] \n");
	printf("Parameters:\n");
	printf("\ttrk_ptid_gdf\n");
	printf("\t\ttracked data gdf file with prefix \"trk\"\n");
	printf("\tpt_gdf\n");
	printf("\t\tfilename of the resulted pt data.\n"); 
	printf("\t\tif not specified, filename will be the one\n");
	printf("\t\twithout the prefix \"trk\"\n");
	exit (0);
}
