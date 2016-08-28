#include "normalmode.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv[])
{
	if (argc==5)
	{
		e_plot_config(argv[1], atoi(argv[4]), 
				atof(argv[2]), atof(argv[3]), false);
	}
	else if (argc==6 && argv[5][0]=='-' && argv[5][1]=='t')
	{
		e_plot_config(argv[1], atoi(argv[4]), 
				atof(argv[2]), atof(argv[3]), true);
	}
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
	printf("\tPlot configuration of ellipsoids crystals/glasses\n"
		   "\tTo calculate normal modes, use 'e_normalmode'\n"
		   "Usage:\n"
		   "\te_plot_config gdf_data_file a b n [Tracked]\n"
		   "Parameters:\n"
		   "\tgdf_data_file\n"
		   "\t\tgdf data\n"
		   "\ta, b\n"
		   "\t\tsemi-long and semi-short axis length\n"
		   "\t\tunit: pixels\n"
		   "\tn\n"
		   "\t\tchoose to plot the n'th frame\n"
		   "\tTracked\n"
		   "\t\t-t\ttracked\n"
		   "\t\tDefault\tuntracked\n");
	exit (0);
}
