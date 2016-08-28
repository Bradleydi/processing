#include "colloid_base.h"
#include "io.h"
#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	if (argc==3)
	{
		int n=atoi(argv[2]);
		plot_mode(n, argv[1]);
	}
	else if (argc==4 && argv[3][0]=='-' && 
			argv[3][1]=='e' && argv[3][2]=='p' && argv[3][3]=='s') 
	{
		int n=atoi(argv[2]);
		plot_mode_eps(n, argv[1]);
	}
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tplot_mode gdf_data_file n [-eps]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\tn\n");
	printf("\t\tchoose to plot the n'th mode\n");
	printf("\t-eps\n"
		   "\t\tsave the plot as eps file\n");
	exit (0);
}
