#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>


bool is_pt(const char *file)
{
	colloid_base pt;
	char *gdffile = getfilename (file, ".gdf");
	readgdf(pt, gdffile);
	int *size=pt.get_size();
	float *pt_ptr=pt.get_array_pointer();

	// check the last column
	pt_ptr+=size[0]-1;
	int T = (int)(*pt_ptr);
	int tmp, i;
	for (i=1; i<size[1]; i++)
	{
		tmp=(int)(*pt_ptr);
		if ( T > tmp ) // no decrease
		{
			printf("# '%s' is NOT valid position data.\n"
				   "# frame no. should not decrease.\n", gdffile);
			return false;
		}
		if ( (T+1) < tmp ) // no gap
		{
			printf("# '%s' is NOT valid position data.\n"
				   "# should be NO frame no. gap.\n", gdffile);
			return false;
		}
		if ( T!=tmp )
			T=tmp;
		pt_ptr += size[0];
	}
	printf("# '%s' is valid position data.\n", gdffile);
	delete [] size;
	free(gdffile);
	return true;
}


bool is_ptid(const char *file)
{
	colloid_base pt;
	char *gdffile = getfilename (file, ".gdf");
	readgdf(pt, gdffile);
	int *size=pt.get_size();
	float *pt_ptr=pt.get_array_pointer();

	// check the time and index
	pt_ptr+=size[0]-2;
	int T = (int)(*pt_ptr);
	int I = (int)(*(pt_ptr+1));
	int tmp, i;
	for (i=1; i<size[1]; i++)
	{
		tmp=(int)(*(pt_ptr+1));
		if ( I > tmp ) // no decrease
		{
			printf("# '%s' is NOT valid tracked data.\n"
				   "# particle no. should not decrease.\n", gdffile);
			return false;
		}
		if ( (I+1) < tmp ) // no gap
		{
			printf("# '%s' is NOT valid tracked data.\n"
				   "# should be NO particle no. gap.\n", gdffile);
			return false;
		}
		if ( I==tmp )  // same particle: no decrease, there can be gap
		{
			tmp=(int)(*pt_ptr);
			if ( T > tmp ) // no decrease
			{
				printf("# '%s' is NOT valid tracked data.\n"
					   "# frame no. should not decrease.\n", gdffile);
				return false;
			}
			else if ( T < tmp )
				T=tmp;
		}
		else
		{
			I=tmp; T=(int)(*pt_ptr);
		}
		pt_ptr += size[0];
	}
	printf("# '%s' is valid tracked data.\n", gdffile);
	delete [] size;
	free(gdffile);
	return true;
}
