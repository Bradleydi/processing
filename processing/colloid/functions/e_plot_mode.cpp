#include "normalmode.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void phelp();

int main (int argc, char *argv[])
{
	if (argc==3)
	{
		int n=atoi(argv[2]);
		e_plot_mode(argv[1], n);
	}
	else if (argc==6) 
	{
		int n=atoi(argv[5]);
		e_plot_mode2(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), n);
	}
	else if (argc==7) 
	{
		int n=atoi(argv[5]);
		e_plot_mode3(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]), n, atof(argv[6]));
	}
	else if (argc==8 && strcmp(argv[7], "-T")==0 ) 
	{
		int n=atoi(argv[5]);
		e_plot_mode_trans(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]),
				n, atof(argv[6]));
	}
	else if (argc==8 && strcmp(argv[7], "-R")==0 ) 
	{
		int n=atoi(argv[5]);
		e_plot_mode_rotat(argv[1], atof(argv[2]), atof(argv[3]), atof(argv[4]),
				n, atof(argv[6]));
	}
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
