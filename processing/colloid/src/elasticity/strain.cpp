#include "elasticity.h"
#include "data_preprocess.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

void strain( const char *file, 
				  const float Lambda,
				  bool Remove_Drift)
{
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);
	
	int i, j, k, t, p;
	int *size=ptid.get_size();
	int ncol=size[0];
	int ti=ncol-2; //time index in the array (=dim)
	int ii=ncol-1; //index index in the array
	int dim=ti; // dimension
	float *ptid_ptr=ptid.get_array_pointer();
	float *pointer, *pointer1;

	// total particle number
	// = last element of ptid+1
	const int total_p=(int)ptid_ptr[size[2]-1]+1;

	// frame number per particle
	int *frame_ptcl =Calloc(int, total_p); POINTER_NULL(frame_ptcl);

	pointer=ptid_ptr+ti;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		++frame_ptcl[(int)(*(pointer+1))];
		pointer+=ncol;
	}
	// total time
	const int total_t=(int)(tt)+1;
	
	if (total_p*total_t > 2*size[1])
		pERROR("too many broken trajectories.");
	
	if (Remove_Drift)
		remove_drift(ptid, dim);

	//=====================================================================
	// get mean position
	
	float * mp_ptr = (float *)calloc(dim*total_p, sizeof(float));
	POINTER_NULL(mp_ptr);

	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		for (k=0; k<dim; k++)
			*(pointer1+k)+=*(pointer+k);
		pointer+=size[0];
	}

	// mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
		for (k=0; k<dim; k++)
			*(pointer++)/=total_t;  //frame_ptcl[p];

	
	//=====================================================================
	// get size of the system
	float xmin=*mp_ptr;
	float xmax=xmin;
	float ymin=*(mp_ptr+1);
	float ymax=ymin;
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
	{
		if ( xmin > *pointer )
			xmin = *pointer;
		else if ( xmax < *pointer )
			xmax = *pointer;

		++pointer;

		if ( ymin > *pointer )
			ymin = *pointer;
		else if ( ymax < *pointer)
			ymax = *pointer;

		++pointer;
	}
	
	//======================================================================
	// get derivative
	//float Lambda=1.1;  // coarse_graining_length, I'll write a code to search it
					// automatically with initial input to make average neighbor
					// number as 6.
	float Lambda2=Lambda*Lambda;
	float dd; // distance^2
	float dx, dy;
	
	int p1;
	int *neighbor_count=Malloc(int, total_p); POINTER_NULL(neighbor_count);
	// find neighbor, based on the mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
	{
		neighbor_count[p]=0;
		
		pointer1=mp_ptr;
		for(i=0; i<total_p; i++)
		{
			dd=0;
			for (j=0; j<dim; j++)
			{
				dx=*(pointer+j)-*(pointer1+j);
				dd+=dx*dx;
				if (dd>Lambda2)
					break;
			}
			if (dd<=Lambda2)
				++neighbor_count[p];
			pointer1+=dim;
		}
		--neighbor_count[p]; // minus itself
		pointer+=dim;
	}

	// average neighbor
	unsigned long sum2=0L;
	int *csum=Malloc(int, total_p+1); POINTER_NULL(csum);
	csum[0]=0;
	unsigned int NANcount=0;
	for (i=0; i<total_p; i++)
	{
		j=neighbor_count[i];
		if (j==0)
		{
			++NANcount;
			printf("# Error: no neighbor for particle id\t=\t%d\n", i);
		}
		else
		{
			csum[i+1]=csum[i]+j;
			sum2+=j*j;
		}
	}
	if (NANcount!=0)
		printf("# Total particle with no neighbor:\t%d\n", NANcount);
	double ave=(double)(csum[total_p])/total_p;
	printf("# max distance setted  :\t%.6f\n"
		   "# average neighbor     :\t%.6f\n"
		   "# standard derivative  :\t%.6f\n",
		   Lambda, ave, sqrt((double)(sum2)/total_p-ave*ave));
	
	// get neighbor id
	int *neighbor_id=Malloc(int, csum[total_p]); POINTER_NULL(neighbor_id);
	pointer=mp_ptr;
	i=0;
	for (p=0; p<total_p; p++)
	{
		pointer1=mp_ptr;
		for(k=0; k<total_p; k++)
		{
			dd=0;
			for (j=0; j<dim; j++)
			{
				dx=*(pointer+j)-*(pointer1+j);
				dd+=dx*dx;
				if (dd>Lambda2)
					break;
			}
			if (dd<=Lambda2 && k!=p)
				neighbor_id[i++]=k;
			pointer1+=dim;
		}
		pointer+=dim;
	}

	free(neighbor_count);

	
	//==============================================================
	int ncg=dim*(dim+1);
	colloid_base cg;
	cg.reserve_memory(ncg, size[1]);
	float* cg_ptr=cg.get_array_pointer();
	for (i=0; i<ncg*size[1]; i++)
		cg_ptr[i]=0.0;
	
	//================================================================
	printf("# Least squares\n");
	
	// least square
	// df=(\p f / \p x)dx + (\p f / \p y) dy
	float Xxx, Xxy, Xyy, Yxx, Yxy, Yyx, Yyy;
	float dx0, dy0;
	float delta;

	int * index=Malloc(int, total_t*total_p);
	for (p=0; p<total_t*total_p; p++)
		index[p]=-1;

	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		index[(int)(*(pointer+1)*total_t+*pointer)]=i;
		pointer+=ncol;
	}


	float exx, exy, eyx, eyy; // partial derivative
	float chix, chiy; // nonaffineness
	//======================================================
	
	pointer=cg_ptr;
	pointer1=ptid_ptr;
	float *pointer2;
	for (i=0; i<size[1]; i++)
	{
		p=(int)(*(pointer1+ii));

		if ((csum[p]+1) >= csum[p+1]) // no enough equations
		{
			*pointer=1.0e6;
			pointer+=ncg;
			pointer1+=ncol; // without this there will be a bug
			continue;
		}

		t=(int)(*(pointer1+ti));
		//pointer=cg_ptr+i*ncg;
		
		Xxx=0; Xxy=0; Xyy=0;
		Yxx=0; Yxy=0; Yyx=0; Yyy=0;
		for (k=csum[p]; k<csum[p+1]; k++)
		{
			p1=neighbor_id[k];
			j=index[p1*total_t+t];	
			if (j!=-1)
			{
				pointer2=ptid_ptr+j*ncol;
				// for 2D only
				dx0=mp_ptr[p*dim]-mp_ptr[p1*dim];
				dy0=mp_ptr[p*dim+1]-mp_ptr[p1*dim+1];
				Xxx+=dx0*dx0;
				Xxy+=dx0*dy0;
				Xyy+=dy0*dy0;
				dx=*pointer1-*pointer2;
				dy=*(pointer1+1)-*(pointer2+1);
				Yxx+=dx0*dx;
				Yxy+=dy0*dx;
				Yyx+=dx0*dy;
				Yyy+=dy0*dy;
			}
		}
		// solve equation
		delta=Xxx*Xyy-Xxy*Xxy;
		if (delta==0)
			printf("# Error: divided by zero.\n"
                   "# Check neighbors for each particle!\n");
		exx=(Yxx*Xyy-Yxy*Xxy)/delta;
		exy=(Yxy*Xxx-Yxx*Xxy)/delta;
		eyx=(Yyx*Xyy-Yyy*Xxy)/delta;
		eyy=(Yyy*Xxx-Yyx*Xxy)/delta;
			
		*(pointer++)=exx-1;
		*(pointer++)=exy;
		*(pointer++)=eyx;
		*(pointer++)=eyy-1;
			
		// nonaffineness
		chix=0; chiy=0;
		for (k=csum[p]; k<csum[p+1]; k++)
		{
			p1=neighbor_id[k];
			j=index[p1*total_t+t];
			if (j!=-1)
			{
				pointer2=ptid_ptr+j*ncol;
				dx0=mp_ptr[p*dim]-mp_ptr[p1*dim];
				dy0=mp_ptr[p*dim+1]-mp_ptr[p1*dim+1];
				dx=*pointer1-*pointer2;
				dy=*(pointer1+1)-*(pointer2+1);

				delta=dx-exx*dx0-exy*dy0;
				chix+=delta*delta;
				delta=dy-eyx*dx0-eyy*dy0;
				chiy+=delta*delta;
			}
		}
		*(pointer++)=chix;
		*(pointer++)=chiy;
		//*(pointer++)+=chix+chiy;

		pointer1+=ncol;
	}

	//======================================================
	char *cgfile=getfilename(file, "_stn.gdf");
	writegdf(cg, cgfile);
	free(cgfile);

	delete [] size;
	free(frame_ptcl);
	free(mp_ptr);
	free(csum);
	free(neighbor_id);
	free(index);
}
