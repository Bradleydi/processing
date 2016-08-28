#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv[])
{
	if (argc==2)
		plot_K(argv[1]);
	/*
	else if (argc==4 && argv[3][0]=='-' && 
			argv[3][1]=='e' && argv[3][2]=='p' && argv[3][3]=='s') 
	{
		int n=atoi(argv[2]);
		plot_mode_eps(n, argv[1]);
	}*/
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("\tPlot normal modes of ellipsoids crystals/glasses\n"
		   "\tTo calculate normal modes, use 'e_normalmode'\n"
		   "Usage:\n"
		   "\te_plot_mode gdf_data_file n\n"
		   "Parameters:\n"
		   "\tgdf_data_file\n"
		   "\t\ttracked gdf data\n"
		   "\tn\n"
		   "\t\tchoose to plot the n'th mode\n");
	/*
	printf("Usage:\n");
	printf("\tplot_mode gdf_data_file n [-eps]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\tn\n");
	printf("\t\tchoose to plot the n'th mode\n");
	printf("\t-eps\n"
		   "\t\tsave the plot as eps file\n");
	*/
	exit (0);
}
