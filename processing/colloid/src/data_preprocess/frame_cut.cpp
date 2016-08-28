#include "data_preprocess.h"
#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

void frame_cut(const char *file, const int maxt, bool Tracked)
{
	char *gdffile = getfilename(file, ".gdf");
	colloid_base cb;
	readgdf(cb, gdffile);

	int *size=cb.get_size();
	float *cb_ptr=cb.get_array_pointer();

	int ti=size[0]-1;
	if (Tracked)
		--ti;

	int i;
	float *pointer=cb_ptr+ti;
	float tt = *pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( tt < *pointer )
			tt = *pointer;
	}
	const int total_t = (int)tt + 1;

	if ( maxt > total_t )
	{
		fprintf(stderr, "# ERROR: total frame number is %d\n", total_t);
		exit (1);
	}

	int good=0;
	tt = (float)maxt;
	pointer=cb_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( *pointer < tt )
			++good;
		pointer += size[0];
	}

	colloid_base newcb;
	newcb.reserve_memory(size[0], good);
	float *newcb_ptr=newcb.get_array_pointer();
	
	int k;
	pointer=cb_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		if ( *pointer < tt )
			for (k=0; k<size[0]; k++)
				*(newcb_ptr++) = *(cb_ptr+k);
		pointer += size[0];
		cb_ptr += size[0];
	}
	delete [] size;

	char suffix [20];
	sprintf(suffix, "_fc%d.gdf", maxt);
	char *newgdffile = getfilename (file, suffix);
	writegdf(newcb, newgdffile);

	char *infofile = getfilename (newgdffile, ".info");
	FILE *infof = fopen (infofile, "a");
	fprintf(infof, "get from\n\tframe_cut %s %d\n\n", gdffile, maxt);
	fclose(infof);

	free(gdffile);
	free(newgdffile);
	free(infofile);
}
