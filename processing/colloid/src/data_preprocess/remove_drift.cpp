#include "data_preprocess.h"
#include "colloid_base.h"
#include "miscellaneous.h"
#include "io.h"

#include <cstdio>
#include <cstdlib>

using namespace std;

// only valid for linear drift
void remove_drift(colloid_base &ptid, int dim)
{
	printf("# Removing drift ...\n");
	int *size=ptid.get_size();
	float *ptid_ptr=ptid.get_array_pointer();
	const int ti = size[0]-2;
	const int ii = size[0]-1;

	/* test tracked data */
	if ( dim == 0 )
		dim = ti;
	else
	{
		if ( dim < 0 )
			pERROR("given dimension less than 0!");
		if ( dim+2 > size[0] )
			pERROR("not a tracked data!");
	}
	
	if (dim!=2)
		pERROR("not implimented yet!");
	
	/* get total_p */
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	int *frame_ptcl = (int *) calloc ( total_p, sizeof(int) );
	POINTER_NULL(frame_ptcl);


	/* get total_t */
	int i;
	float *pointer=ptid_ptr+ti;
	++frame_ptcl[(int)(*(pointer+1))];
	float tt=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer += size[0];
		if ( tt < *pointer )
			tt = *pointer;
		++frame_ptcl[(int)(*(pointer+1))];
	}

	const int total_t=(int)tt+1;

	/*=============================
	 * find out particles existing in all frames
	 *=============================*/
	unsigned char *good=(unsigned char *)calloc(total_p, sizeof(unsigned char));
	POINTER_NULL(good);
	int all_p=0;
	for (i=0; i<total_p; i++)
	{
		if ( frame_ptcl[i] == total_t )
		{
			++good[i];
			++all_p;
		}
	}
	free(frame_ptcl);
	printf("# number of particles                       : %d\n", total_p);
	printf("# number of particles existing in all frames: %d\n", all_p);
	if ( all_p < 2*total_p/3)
		pERROR("bad tracking, too many broken trajectories!");

	double *xdrift=(double *)calloc(total_t, sizeof(double));
	double *ydrift=(double *)calloc(total_t, sizeof(double));
	POINTER_NULL(xdrift);
	POINTER_NULL(ydrift);

	pointer=ptid_ptr;
	int t;
	for (i=0; i<size[1]; i++)
	{
		if ( good [ (int)(*(pointer+ii)) ] )
		{
			t=(int)(*(pointer+ti));
			xdrift[t] += (double)(*pointer);
			ydrift[t] += (double)(*(pointer+1));
		}
		pointer += size[0];
	}
	free(good);
	
	//double xc0=(*xdrift)/all_p;
	//double yc0=(*ydrift)/all_p;
	for (t=0; t<total_t; t++)
	{
		//xdrift[t]=xdrift[t]/all_p-xc0;
		//ydrift[t]=ydrift[t]/all_p-yc0;
		xdrift[t]=xdrift[t]/all_p;
		ydrift[t]=ydrift[t]/all_p;
	}
	/*============================================================*/
	/* linearize the drift */

	double ax, bx, ay, by;
	linearfit(total_t, xdrift, ax, bx);
	linearfit(total_t, ydrift, ay, by);

	for (t=0; t<total_t; t++)
	{
		//xdrift[t]=ax*t+bx+xc0;
		//ydrift[t]=ay*t+by+yc0;
		xdrift[t]=ax*t+bx;
		ydrift[t]=ay*t+by;
	}
	
	/*============================================================*/
	/* remove the drift */
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		*pointer     -= (float)(xdrift[t]);
		*(pointer+1) -= (float)(ydrift[t]);
		pointer += size[0];
	}
	
	free(xdrift);
	free(ydrift);
	delete [] size;
}
