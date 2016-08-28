#include "data_preprocess.h"

#include "colloid_base.h"

#include <cstdlib>

using namespace std;

void ptid2pt(colloid_base &ptid, colloid_base &pt)
{
	int i, j;
	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();
	

	const int ti=size[0]-2;
	float *pointer=ptid_ptr+ti;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer += size[0];
	}

	const int total_t=(int)tt+1;
	
	int *ptcl_frame = (int *)calloc(total_t, sizeof(int));
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer += size[0];
	}

	int *cum_ptcl_frame = (int *)malloc((total_t+1)*sizeof(int));
	cum_ptcl_frame[0]=0;
	for ( i=0; i<total_t; i++)
	{
		cum_ptcl_frame[i+1] = cum_ptcl_frame[i] + ptcl_frame[i];
		ptcl_frame[i]=0;
	}

	float **index=(float **)malloc(size[1]*sizeof(float*));
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=(int)(*(pointer+ti));
		index[cum_ptcl_frame[j]+ptcl_frame[j]]=pointer;
		++ptcl_frame[j];
		pointer+=size[0];
	}

	free(ptcl_frame);
	free(cum_ptcl_frame);

	int npt=size[0]-1;
	pt.reserve_memory(npt, size[1]);
	float *pt_ptr=pt.get_array_pointer();
	for (i=0; i<size[1]; i++)
	{
		pointer=index[i];
		for (j=0; j<npt; j++)
			*(pt_ptr++)=*(pointer++);
	}
	free(index);
	delete [] size;
}
