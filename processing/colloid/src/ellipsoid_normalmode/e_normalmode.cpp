#include "normalmode.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "io.h"
#include "lapack.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

/* Here a0 is just a characteristic length, such as the length of short axis
 * , or first peak of g(r).
 * here a,b are the long and short axis with unit pixel
 * it's used for calculating the ratio of rotational "mass" v.s. translational
 * mass. The form is referred as moment of inertia.
 * gdfdata should be x,y,ang,t,id
 *
 * To adapt to the experimental data, one need to consider the missing data.
 * Now I'm trying to learn the statistics involving missing data, however,
 * Prof. Han wants to see the result asap, hence here I only implement a
 * simple method for randomly missing data. Actually, in experiments, data are
 * always missing around the defects, hence as a result, the missing data 
 * are highly correlated both in space and time.
 * Whatever, I need just some method to work now.
 *
 * I don't know whether remove drift works, but I know, ptid_drift doesn't work.
 */
void e_normalmode(const char *file, double a0, double a, double b,
		bool Remove_Drift)
{
	// read gdf data
	char *gdffile=getfilename(file, ".gdf");
	colloid_base ptid;
	readgdf(ptid, gdffile);
	free(gdffile);
	
	int i, j, k, t, p;
	int *size=ptid.get_size();
	const int dim=2; // dimension
	if (dim+3>size[0])
		pERROR("Not a 2D ellipsoid data or tracked data!");
	
	const int ti=size[0]-2; //time index in the array
	const int ii=ti+1; // particle index in the array
	
	float *ptid_ptr=ptid.get_array_pointer();
	float *pointer, *pointer1;

	// rescale 
	if (a0>0)
	{
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			*pointer /= a0;
			*(pointer+1) /= a0;
			pointer+=size[0];
		}
		a /= a0;
		b /= a0;
	}

	// total particle number
	// = last element of ptid+1
	const int total_p=(int)ptid_ptr[size[2]-1]+1;

	pointer=ptid_ptr+ti;
	float tt=0;
	for (i=0; i<size[1]; i++)
	{
		if ( tt < *pointer )
			tt = *pointer;
		pointer+=size[0];
	}
	const int total_t=(int)tt + 1;
	
	printf("# total particle\t\t:\t%d\n", total_p);
	printf("# total frame   \t\t:\t%d\n", total_t);

	int *ptcl_frame=Calloc(int, total_t); POINTER_NULL(ptcl_frame);
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[(int)(*pointer)];
		pointer+=size[0];
	}

	int *cum_ptcl_frame=(int *)malloc((total_t+1)*sizeof(int));
	POINTER_NULL(cum_ptcl_frame);
	cum_ptcl_frame[0]=0;
	cum_ptcl_frame[1]=ptcl_frame[0];
	int maxptcl=ptcl_frame[0];
	for (i=1; i<total_t; i++)
	{
		cum_ptcl_frame[i+1]=cum_ptcl_frame[i]+ptcl_frame[i];
		if ( maxptcl < ptcl_frame[i] )
			maxptcl = ptcl_frame[i];
	}

	if (maxptcl*total_t > 2*size[1])
		pERROR("too many broken trajectories!");

	int *index_now=(int *)calloc(total_t, sizeof(int)); POINTER_NULL(index_now);
	float **index=(float **)malloc(size[1]*sizeof(float *));
	POINTER_NULL(index);
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		t=(int)(*(pointer+ti));
		index[cum_ptcl_frame[t]+(index_now[t]++)]=pointer;
		pointer+=size[0];
	}
	free(index_now);

	// whether goodframe, remove bad frame.
	// cannot missing too many
	int total_gt=0; // total_goodframe
	const int goodframe_limit = (int)(0.9*size[1]/total_t); 
	for (i=0; i<total_t; i++)
		if (ptcl_frame[i]>goodframe_limit) ptcl_frame[total_gt++]=i;
	// here ptcl_frame saves frame_id which is a goodframe
	int *goodframe=ptcl_frame;

	// check whether goodptcl
	int total_gp=0; // total_goodptcl
	const int goodptcl_limit = (int)(0.90*total_gt); 
	int *frame_ptcl=Calloc(int, total_p); POINTER_NULL(frame_ptcl);
	for (j=0; j<total_gt; j++)
	{
		t=goodframe[j];
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
			++frame_ptcl[(int)(*(index[i]+ii))];
	}
	int *goodptcl=Malloc(int, total_p); POINTER_NULL(goodptcl);
	double delta=0.0;
	for (p=0; p<total_p; p++)
		if (frame_ptcl[p]>goodptcl_limit) 
			{goodptcl[p]=(total_gp++); delta+=frame_ptcl[p];}
		else goodptcl[p]=-1;
	// goodptcl !=-1 if the ptcl is a good one
	delta /= total_gp*total_gt;
	printf("# total good particle\t\t:\t%d\n", total_gp);
	printf("# total good frame   \t\t:\t%d\n", total_gt);
	

	if ( total_gt < 3*total_gp )
	{
		fprintf(stderr, "# Error: there should be enough independent configurations!\n");
		fprintf(stderr, "#        frame_number >= 3*particle_number\n");
		exit (1);
	}
	if (delta <=0.95)
		pERROR("too small delta. Track better");
	
	if (Remove_Drift)
		remove_drift(ptid, dim);
	//=====================================================================
	// get mean position
	// I rearranged that data structure
	// [x1, x2, x3,..., y1, y2, y3,..., ang1, ang2, ang3,...]
	const int total_gp2 = 2*total_gp;
	int Gn=total_gp*3;
	float * mp_ptr = Calloc(float, Gn); POINTER_NULL(mp_ptr);

	for (j=0; j<total_gt; j++)
	{
		t=goodframe[j];
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			p=goodptcl[(int)(*(pointer+ii))]; // reindex the good ptcl
			if ( p>=0 )
			{
				pointer1=mp_ptr+p;
				*pointer1 += *pointer;                    // x
				*(pointer1+total_gp) += *(pointer+1);      // y
				*(pointer1+total_gp2) += *(pointer+2);     // ang
			}
		}
	}

	// mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
		if (goodptcl[p]>=0)
		{
			*pointer /= frame_ptcl[p];
			*(pointer+total_gp) /= frame_ptcl[p];
			*(pointer+total_gp2) /= frame_ptcl[p];
			++pointer;
		}
	e_writemp(total_gp, mp_ptr, file);
	//=====================================================================
	// get displacement
	double *u = Malloc(double, Gn); POINTER_NULL(u);
	// get covariance matrix
	printf("# getting covariance matrix...\n");
	double *G=Calloc(double, Gn*Gn); POINTER_NULL(G);
	printf("%d\n", Gn);
	char UPLO='U'; double alpha=1.0; int INCX=1;
	for (j=0; j<total_gt; j++)
	{
		t=goodframe[j];
		for (i=0; i<Gn; i++) u[i]=0.0;
		for (i=cum_ptcl_frame[t]; i<cum_ptcl_frame[t+1]; i++)
		{
			pointer=index[i];
			p=goodptcl[(int)(*(pointer+ii))]; // reindex the good ptcl
			// get displacement
			if ( p>=0 )
			{
				pointer1=mp_ptr+p;
				u[p] = (double)(*pointer - *pointer1);                 // x
				u[p+total_gp] = (double)(*(pointer+1) - *(pointer1+total_gp));   // y
				u[p+total_gp2] = (double)(*(pointer+2) - *(pointer1+total_gp2));  // ang
			}
			// if the particle is not in the frame, u[p]=0.0
		}
		// calculate uu^T
		dsyr_(&UPLO, &Gn, &alpha, u, &INCX, G, &Gn);
	}
	free(mp_ptr);
	free(u);
	free(ptcl_frame);
	free(cum_ptcl_frame);
	free(frame_ptcl);
	free(goodptcl);
	free(index);
	ptid.free_memory(); // free the memory
	delete [] size;
	//show_matrix(Gn, G);

	printf("%f %f %f %f %f\n", G[1], G[Gn], G[Gn+2], G[0], G[2*Gn+1]);
	//============================================================
	// to make the covariance matrix correct...
	const int total_pGn = total_gp*Gn;
	double T=0., O=0.;
	for (i=0; i<total_gp2; i++)
		T += G[i*Gn+i];
	for (i=total_gp2; i<Gn; i++)
		O += G[i*Gn+i];
	double I = T/2./O;
	//double I = (a*a+b*b)/2.0;
	double sqrtI = sqrt(I);
	for (i=0; i<total_gp; i++)
	{
		//for (j=0; j<=i; j++) // i should be larger than j, hence choose
			// XX YX YY AX AY AA, the former index should be larger
			// here should be j from 0 to total_gp, since xiyj != xjyi
		for (j=0; j<total_gp; j++)
		{
			k=i*Gn+j;  // XX
			G[k] /= total_gt; // XX 
				G[k+total_gp] /= total_gt; // XY
				G[k+total_gp2] *= sqrtI/total_gt;  // XA
			k += total_pGn;  // YX
			G[k] /= total_gt; // YX 
			G[k+total_gp] /= total_gt; // YY
				G[k+total_gp2] *= sqrtI/total_gt;  // YA
			k += total_pGn;  // AX
			G[k] *= sqrtI/total_gt; // AX
			G[k+total_gp] *= sqrtI/total_gt; // AY
			G[k+total_gp2] *= I/total_gt;  // AA
			/*if (j!=i)
			{
				k=i*Gn+i+total_gp; // XY
				G[k] /= total_gt; // XY 
				G[k+total_gp] *= sqrtI/total_gt; // XA
				G[k+total_gp+total_pGn] *= sqrtI/total_gt; // YA
			}*/
		}
		/*
		for (j=0; j<=i; j++) // i should be larger than j, hence choose
			// XX YX YY AX AY AA, the former index should be larger
		{
			k=i*Gn+j;  // XX
			G[k] /= total_gt;
			G[k+total_gp] /= total_gt;
			G[k+total_gp2] *= sqrtI/total_gt;
			k += total_pGn;  // YX
			
			G[k+total_gp] /= total_gt;
			G[k+total_gp2] *= sqrtI/total_gt;
			k += total_pGn;  // AX
			G[k+total_gp2] *= I/total_gt; 
			*//*
			if(j!=i)
			{
				G[k] *= sqrtI/total_gt;
				G[k+total_gp] *= sqrtI/total_gt;
				G[k - total_pGn] /= total_gt;
			}*//*
		}
		//k = i*Gn+j+2*total_pGn;  // AX
		G[k] *= sqrtI/total_gt;
		G[k+total_gp] *= sqrtI/total_gt;
		G[k - total_pGn] /= total_gt; */
	}
	//show_matrix(Gn, G);

	delta=1.0; // this is actually wrong, but I don't know how to make it better
	// I don't know how to deal with missing data
	double delta_1=1.0/delta, delta_2=delta_1*delta_1;
	for(i=0; i<Gn; i++)
	{
		G[i*Gn+i] *= delta_1;
		for (j=0; j<i; j++)
		{
			k=i*Gn+j;
			G[k] *=delta_2;
			G[j*Gn+i]=G[k];
		}
	}
	printf("%f %f %f %f %f\n", G[1], G[Gn], G[Gn+2], G[0], G[2*Gn+1]);
	//printf("# DONE!\n");
	
	// write G to file.dcm
	e_writedcm(Gn, G, file);
	//show_matrix(Gn, G);

	double *E=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(E);
	printf("# Calculating eigenvalues...\n");
	dsyev(G, Gn, E);

	if (E[0]<0)
	{
		printf("# First eigenvalue: %6.10f\n", E[0]);
		pERROR("negative eigenvalue(s)");
	}

	// write E and G to file.ev
	e_writeev(Gn, E, G, file);

	free(G);
	free(E);
}
