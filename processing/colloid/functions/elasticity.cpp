#include "colloid_base.h"
#include "io.h"
#include "elasticity.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char *argv[])
{
	//char filename []="simulation8000.gdf";
	colloid_base ptid;
	readgdf(ptid, argv[1]);

	if (argc==3)
	{
		char outfile [] = "tmp_elasticity.dat";
		elasticity(ptid, /* b= */ 10, /* Lambda= */ atof(argv[2]) /*0.15*/ /*076*/, /* nbin= */ 100, outfile);
	}
	else if (argc==4)
	{
		char outfile [] = "nbin.dat";
		elasticity(ptid, /* b= */ atoi(argv[3]), /* Lambda= */ atof(argv[2]) /*0.15*/ /*076*/, /* nbin= */ 100, outfile);
	}
	else if (argc==5)
		elasticity(ptid, /* b= */ atoi(argv[3]), /* Lambda= */ atof(argv[2]) /*0.15*/ /*076*/, /* nbin= */ 100, /*outfile*/ argv[4]);

	return 0;
}
