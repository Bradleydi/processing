#include "structure.h"

#include <cstdio>
#include <cstdlib>


using namespace std;

void phelp();


int main (int argc, char *argv[])
{
	if ( argc == 1 ) phelp();
	bool sk2D=false, tracked=false;
	int i;
	for (i=2; i<argc; i++)
	{
		if (argv[i][0] == '-')
			if ( argv[i][1] == '2' )
				sk2D = true;
			else if ( argv[i][1] == 't' )
				tracked = true;
			else
				phelp();
		else
			phelp();
	}
	sk(argv[1], sk2D, tracked);
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tsk gdffile [-2] [-t]\n");
	printf("Purpose:\n"
		   "\t\tstructure factor\n");
	printf("Parameters:\n");
	printf("\tgdffile\n");
	printf("\t\tdata gdf file\n");
	printf("\t-2\n");
	printf("\t\tshow 2D structure factor, defalt is angular averaged\n"); 
	printf("\t-t\n");
	printf("\t\tracked data\n");
	exit (0);
}
