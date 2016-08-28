#include "colloid_base.h"
#include "io.h"
#include "dynamics.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void phelp();

int main (int argc, char *argv[])
{
	//char filename []="trkfeng8000.gdf";
	colloid_base ptid;
	readgdf(ptid, argv[1]);

	//cout << argc << endl;

	//void dxpdf(colloid_base& ptid, const int & dt, const int & nbin)
	if (argc==2)
		msd(ptid);
	else if (argc==3)
	{
		if (argv[2][0]=='-' && argv[2][1]=='r')
			msd(ptid, 2, false);
		else
			msd(ptid, atoi(argv[2]));
	}
	else if (argc==4)
		msd(ptid, atoi(argv[2]), false);
	else
		phelp();

	return 0;
}


void phelp()
{
	printf("Usage:\n");
	printf("\tmsd gdf_data_file [dim] [remove_drift]\n");
	printf("Parameters:\n");
	printf("\tgdf_data_file\n");
	printf("\t\ttracked data\n");
	printf("\tdim\n");
	printf("\t\tdimension\n");
	printf("\tremove_drift\n");
	printf("\t\t-r remove the drift\n");
	exit (0);
}
