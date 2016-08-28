#include "colloid_base.h"
#include "data_preprocess.h"

#include <cstdio>
#include <cstdlib>


using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	//char filename []="simulation8000.gdf";
	if (argc==1)
		phelp();
	if (argc==2)
		goodtrack(argv[1], /*threshold=*/ 1.0);
	else if (argc==3)
		goodtrack(argv[1], /*threshold=*/ atof(argv[2]));
	else
		phelp();
	return 0;
}

void phelp()
{
	printf("Usage:\n");
	printf("\tgoodtrack tracked_file [threshold]\n");
	printf("Parameters:\n");
	printf("\tthreshold\n");
	printf("\t\tIn a given time window, particle with trajectory length\n");
	printf("\t\tlonger than threshold*window_width, it's a good particle\n"); 
	exit (0);
}
