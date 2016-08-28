#include "normalmode.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void phelp();

int main (int argc, char *argv[])
{
	if (argc < 3)
		phelp();
	else if (argc == 3)
	{
		if (strcmp(argv[2], "-DOS")==0)
			e_DOS(argv[1]);
		if (strcmp(argv[2], "-DOSLOG")==0)
			e_DOS_log(argv[1]);
		else if (strcmp(argv[2], "-CDOS")==0)
			e_CDOS(argv[1]);
		else if (strcmp(argv[2], "-PR")==0)
			e_PR(argv[1]);
		else if (strcmp(argv[2], "-P")==0)
			e_P(argv[1]);
		else if (strcmp(argv[2], "-CORR")==0)
			e_corr(argv[1]);
		else
			phelp();
	}
	else if (argc == 4)
	{
		// a0 = 1.0, not remove drift
		e_normalmode(argv[1], -1.0, atof(argv[2]), atof(argv[3]), false);
		e_DOS(argv[1]);
		e_CDOS(argv[1]);
		e_PR(argv[1]);
	}
	else if (argc==5)
	{
		if (argv[4][0]=='-' && argv[4][1]=='r')
			e_normalmode(argv[1], -1.0, atof(argv[2]), atof(argv[3]), true);
		else if (argv[4][0] != '-')
			e_normalmode(argv[1], atof(argv[4]), \
					atof(argv[2]), atof(argv[3]), false);
		else
			phelp();
		e_DOS(argv[1]);
		e_CDOS(argv[1]);
		e_PR(argv[1]);
	}
	else if (argc==6 && argv[5][0]=='-' && argv[5][1]=='r' )
	{
			e_normalmode(argv[1], atof(argv[4]), \
					atof(argv[2]), atof(argv[3]), true);
			e_DOS(argv[1]);
			e_CDOS(argv[1]);
			e_PR(argv[1]);
	}
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("\tTo calculate the normal modes for ellipsoidal crystals/glasses\n"
		   "\tAlso calculate the CDOS and PR\n"
		   "\tTo plot each mode, use 'e_plot_mode'\n");
	printf("Usage:\n");
	printf("\tnormal_mode gdf_data_file a b [a0] [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\ta b\n"
		   "\t\tsemi-long and semi-short axis length\n"
		   "\t\tunit: pixels\n");
	printf("\ta0\n");
	printf("\t\tpixel per unit length, for ellipse, use b as new unit\n");
	printf("\tremove_drift\n");
	printf("\t\t-r remove the drift\n");
	printf("Further usage:\n");
	printf("\tnormal_mode gdfdata -options\n");
	printf("Avaible options:\n"
		   "\t-DOS\tdensity of states\n"
		   "\t-CDOS\tcumulative density of states\n"
		   "\t-PR\tparticipation ratio\n");
	exit (0);
}
