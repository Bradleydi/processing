#include "data_preprocess.h"
#include "colloid_base.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>

using namespace std;


void space_cut(const char *file, const char *region)
{
	/*============================================================*/
	/*============================================================*/
	/*============================================================*/
	// test string "region"
	// region should be formated by xmin:xmax:ymin:ymax
	// anyone can be omitted
	int i=0, j=0; 
	int *count=(int *)calloc(4, sizeof(int));
	POINTER_NULL(count);
	while ( region[i] != '\0' )
	{
		if ( region[i] != ':' ) 
		{
			++count[j];
			if ( region[i] < '0') // check is a number
			{
				if ( region[i] != '+' && region[i] !='-' && region[i] != '.')
					pERROR("input is not a number.");
			}
			else if (region[i] > '9' ) 
				pERROR("input is not a number.");
		}
		else
			++j;
		++i;
	}
	if ( j != 3 )
		pERROR("region should be formated by xmin:xmax:ymin:ymax.");
	
	printf("# %d %d %d %d\n", count[0], count[1], count[2],count[3]);
	
	float  xmin = -1.0e6;
	float  xmax = 1.0e6;
	float  ymin = -1.0e6;
	float  ymax = 1.0e6;

	int maxcount=count[0];
	for (i=1; i<4; i++)
		if ( maxcount < count[i] )
			maxcount = count[i];

	char *tmp = (char *)malloc((maxcount+1)*sizeof(char));
	POINTER_NULL(tmp);

	float **m=(float **)malloc(4*sizeof(float*));
	POINTER_NULL(m);
	m[0]=&xmin; m[1]=&xmax; m[2]=&ymin; m[3]=&ymax;
	i=0; j=0; int k;
	for (i=0; i<4; i++)
	{
		if ( count[i]!=0 )
		{
			for (k=0; k<count[i]; k++)
				tmp[k]=region[j++];
			tmp[k]='\0';
			*(m[i])=(float)atof(tmp);
		}
		++j;
	}

	free(tmp);
	free(m);

	if ( xmin >= xmax || ymin >= ymax )
		pERROR("region set is wrong, please check ...");

	printf("# new region is set to:\n"
		   "#\t[%g, %g] x [%g, %g]\n", 
		   xmin, xmax, ymin, ymax);

	/*============================================================*/
	/*============================================================*/
	/*============================================================*/

	char *gdffile = getfilename (file, ".gdf");
	colloid_base p;
	readgdf(p,gdffile);

	float * p_ptr=p.get_array_pointer();
	int * size=p.get_size();
	float * pointer;
	
	char * inregion = (char *)calloc(size[1], sizeof(char));
	POINTER_NULL(inregion);
	pointer=p_ptr;
	k=0;
	for (i=0; i<size[1]; i++)
	{
		if (*pointer >= xmin && *pointer <= xmax &&
			*(pointer+1) >= ymin && *(pointer+1) <= ymax)
		{
			inregion[i]=1;
			++k;
		}
		pointer+=size[0];
	}

	colloid_base pnewdata;
	pnewdata.reserve_memory(size[0], k);
	//float *newdata = (float *)malloc(size[0]*k*sizeof(float));
	float *newdata=pnewdata.get_array_pointer();
	POINTER_NULL(newdata);
	pointer=p_ptr;
	k=0;
	for (i=0; i<size[1]; i++)
	{
		if (inregion[i])
			for (j=0; j<size[0]; j++)
				newdata[k++]=*(pointer++);
		else
			pointer+=size[0];
	}
	free(inregion);
	
	/*=========================================================*/
	/*=========================================================*/
	/*=========================================================*/
	// write to a new gdf file
	char str [20];
	sprintf(str, "_sc.gdf");
	char *newgdffile = getfilename(file, str);
	i=0;
	while ( FILE_EXIST(newgdffile) )
	{
		free(newgdffile);
		sprintf(str, "_sc%d.gdf", i++);
		newgdffile = getfilename(file, str);
	}
	printf("# new gdffile name:   \"%s\"\n", newgdffile);

	writegdf(pnewdata, newgdffile);

	//free(newdata);

	char *infofile = getfilename(newgdffile, ".info");
	FILE *infof = fopen(infofile, "w");
	FILE_NULL(infof, infofile);
	fprintf(infof, "%s get from\n\t%s\n"
			       "by\n\tspace_cut %s %s\n"
				   "new region is set as\n"
		   		   "\t[%g, %g] x [%g, %g]\n\n", 
		   		   newgdffile, gdffile, file, region, xmin, xmax, ymin, ymax);
	fclose(infof);
	
	free(gdffile);
	free(newgdffile);
	free(infofile);

	delete [] size;
}
