#include "colloid_base.h"
#include "io.h"
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

	// test "trk" prefix
	if ( argv[1][0] == 't' &&
		 argv[1][1] == 'r' &&
		 argv[1][2] == 'k' )
	{
		if (argc==2)
		{
			// get pt filename
			int i, n=3;
			while ( argv[1][n] != '\0' )
				++n;

			n-=2;
			char *filename=(char *)malloc(n*sizeof(char));
			
			for (i=0; i<n; i++)
				filename[i]=argv[1][i+3];

			// whether file exist?
			FILE *fp = fopen(filename, "r");
			if ( fp != NULL )
			{
				fclose(fp);
				fprintf(stderr, "# Error: %s exists!\n", filename);
				exit (1);
			}

			colloid_base ptid, pt;
			readgdf(ptid, argv[1]);
			ptid2pt(ptid, pt);
			writegdf(pt, filename);
		}
		else if (argc==3)
		{
			// whether file exist?
			FILE *fp1 = fopen(argv[2], "r");
			if ( fp1 != NULL )
			{
				fclose(fp1);
				fprintf(stderr, "# Error: %s exists!\n", argv[2]);
				exit (1);
			}

			colloid_base ptid, pt;
			readgdf(ptid, argv[1]);
			ptid2pt(ptid, pt);
			writegdf(pt, argv[2]);
		}
		else
		{
			fprintf(stderr, "# Error: wrong argument(s)!\n");
			exit (1);
		}
		return 0;
	}
	else
	{
		fprintf(stderr, "# Error: ptid_gdf file should be with prefix \"trk\"\n");
		return 1;
	}
}

void phelp()
{
	printf("Usage:\n");
	printf("\tptid2pt trk_ptid_gdf [pt_gdf] \n");
	printf("Parameters:\n");
	printf("\ttrk_ptid_gdf\n");
	printf("\t\ttracked data gdf file with prefix \"trk\"\n");
	printf("\tpt_gdf\n");
	printf("\t\tfilename of the resulted pt data.\n"); 
	printf("\t\tif not specified, filename will be the one\n");
	printf("\t\twithout the prefix \"trk\"\n");
	exit (0);
}
