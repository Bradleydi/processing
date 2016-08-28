#include "colloid_base.h"
#include "io.h"
#include "structure.h"

int main(int argc, char * argv[])
{
	//char filename []="trkfeng8000.gdf";
	colloid_base cb;
	readgdf(cb,argv[1]);
	gr2d(cb, true);
	return 0;
}
