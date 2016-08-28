#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>

void phelp ();

int main (int argc, char *argv[])
{
	if (argc==2)
		is_pt (argv[1]);
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tis_pt gdfdata\n");
	printf("Parameters:\n");
	printf("\tgdfdata\n");
	printf("\t\tpretracked data gdf file\n");
	exit (0);
}
