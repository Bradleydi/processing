#include "colloid_base.h"
#include "data_preprocess.h"

#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char * argv[])
{
	//char filename []="trkxyT23d6.gdf";
	trajectory_filter(/*filename*/ argv[1], atoi(argv[2]));
	return 0;
}
