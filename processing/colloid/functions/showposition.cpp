#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

void phelp();

int main (int argc, char *argv [] )
{
	if ( argc == 2 )
		showposition(argv[1], 1, false);
	else if (argc == 3)
		if ( argv[2][0]=='-' && argv[2][1]=='r' ) 
			showposition(argv[1], 1, true);
		else
			showposition(argv[1], atoi(argv[2]), false);
	else if ( argc == 4 && argv[3][0]=='-' && argv[3][1]=='r' ) 
		showposition(argv[1], atoi(argv[2]), true);
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("\tShow particles positions of a given data, all frames are\n"
		   "\toverlaid together.\n"
		   "\tTo plot each frame as an aimation, use 'config_plot'\n");
	printf("Usage:\n");
	printf("\tshowposition gdf_data_file [factor] [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\tdata file\n");
	printf("\tfactor\n"
		   "\t\tto scale the output figures to zoom in and zoom out,\n"
		   "\t\tusually an experimental data always comes from featuring an\n"
		   "\t\timage, hence this parameter needs not to set.\n"
		   "\t\tHowever, for a simulation data, this factor should be set\n"
		   "\t\tto an appropriate value, e.g. 100.\n"
		   "\t\tDefulat value is 1.0\n");
	printf("\tremove_drift\n");
	printf("\t\t-r remove the drift\n");
	exit (0);
}
