#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void phelp();

int main (int argc, char *argv[])
{
	int i;
	for ( i=0; i<argc; i++)
		if ( strcmp(argv[i], "-l")==0 || strcmp(argv[i], "-L")==0 ||
			strcmp(argv[i], "-llr")==0 || strcmp(argv[i], "-LLR")==0 )
		{
			//drift_llr(file, dim, smooth);
			// local linear regression
			int dim=2;
			double smooth=1.0;
			for ( i=0; i<argc; i++)
				if ( strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-D")==0 ||
				   strcmp(argv[i], "-dim")==0 || strcmp(argv[i], "-DIM")==0 )
					dim=atoi(argv[i+1]);
				else if ( strcmp(argv[i], "-s")==0  || 
						strcmp(argv[i], "-smooth")==0)
					smooth=atof(argv[i+1]);

			drift_llr(argv[1], dim, smooth);
			return 0;
		}
	// just linear fit
	if (argc == 2)
		ptid_drift(argv[1], 2);
	else if (argc == 3)
		ptid_drift(argv[1], atoi(argv[2]));
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("\tShow the drift of the data, and check whether need to remove\n"
		   "\t the drift.\n");
	printf("Usage:\n");
	printf("\tptid_drift gdfdata [dim] \n");
	printf("Parameters:\n");
	printf("\tgdfdata\n");
	printf("\t\tTRACKED data gdf file\n");
	printf("\tdim\n");
	printf("\t\tdimension of the data, if not specified, set to 2.\n"); 
	printf("\t\tRecently, dim!=2 has not been implemented yet.\n");
	exit (0);
}
