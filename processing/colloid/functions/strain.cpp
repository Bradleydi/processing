#include "elasticity.h"

#include <cstdio>
#include <cstdlib>

void phelp ();

int main (int argc, char *argv[])
{
	/*
void strain( const char *file, 
				  const float & Lambda,
				  bool Remove_Drift);
				  */
	if (argc==3)
		strain(argv[1], atof(argv[2]), false);
	else if (argc==4)
	{
		if (argv[3][0]=='-' && argv[3][1]=='r')
			strain(argv[1], atof(argv[2]), true);
		else
			strain(argv[1], atof(argv[2]), false);
	}
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tstrain gdfdata Lambda [Remove_Drift]\n");
	printf("Parameters:\n");
	printf("\tgdfdata\n");
	printf("\t\tpretracked data gdf file\n");
	printf("\tLambda\n");
	printf("\t\tif 2 particles' distance is less than Lambda, then they are\n");
	printf("\t\tviewed as neighbors. Usually set it to be 1.2~1.5a0\n");
	printf("\tRemove_Drift\n");
	printf("\t\t-r        remove drift\n");
	printf("\t\tdefault   notremove drift\n");
	exit (0);
}
