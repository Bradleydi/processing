#include "elasticity.h"

#include <cstdio>
#include <cstdlib>

void phelp ();

int main (int argc, char *argv[])
{
	/* void strain_plot(const char *file, int frame); */
	if (argc==3)
		plot_strain(argv[1], argv[2], "-b", 2);
	else if (argc==4)
		plot_strain(argv[1], argv[2], argv[3], 2);
	else if (argc==5 && argv[3][0]=='-' && argv[3][1]=='f')
		plot_strain(argv[1], argv[2], "-b", atoi(argv[4]));
	else if (argc==6 && argv[4][0]=='-' && argv[4][1]=='f')
		plot_strain(argv[1], argv[2], argv[3], atoi(argv[5]));
	/*
	   else if (argc==4)
	{
		if (argv[3][0]=='-' && argv[3][1]=='r')
			strain(argv[1], atof(argv[2]), true);
		else
			strain(argv[1], atof(argv[2]), false);
	}*/
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tstrain_plot gdfdata frameindex [which_strain]\n");
	printf("Parameters:\n");
	printf("\tgdfdata\n");
	printf("\t\tpretracked data gdf file\n");
	printf("\tframeindex\n"
		   "\t\tframeindex formatted string, to specify\n"
		   "\t\twhich frame should be plot\n");
	printf("\twhich_strain\n");
	printf("\t\t-b\tbulk strain (default)\n");
	printf("\t\t-s\tshear strain\n");
	printf("\t\t-o\tlocal rotation\n");
	printf("string 'frameindex' format:\n"
		   "\tN\n"
		   "\t\tset chosen frame as N (only one frame)\n"
		   "\tN1:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tThe increment is set to 1\n"
		   "\t\tA ':' means choosing all frames.\n"
		   "\tN1:c:N2\n"
		   "\t\tN1: lower bound, if not specified, set to lower\n"
		   "\t\tN2: upper bound, if not specified, set to upper\n"
		   "\t\tc: increment, if not specified, set to 1\n");
	exit (0);
}
