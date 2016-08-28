#include "elasticity.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

void coarse_grain(colloid_base& cg, 
		          const colloid_base& ptid, 
				  const int& b, 
				  const float & Lambda,
				  bool Remove_Drift)
{
	
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
	int total_p=(int)ptid_ptr[size[2]-1]+1;

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
	int total_t=(int)(tt)+1;
	
	int *ptcl_frame;
	float *drift=NULL;
	if (Remove_Drift)
	{
		//particle number per frame
		ptcl_frame=Calloc(int, total_t); POINTER_NULL(ptcl_frame);
		//particle number per frame
		pointer=ptid_ptr+ti;
		for (i=0; i<size[1]; i++)
		{
			++ptcl_frame[(int)(*pointer)];
			pointer+=ncol;
		}
		// get drift sum
		drift=Calloc(float, dim*total_t); POINTER_NULL(drift);
		pointer=ptid_ptr;
		for (i=0; i<size[1]; i++)
		{
			pointer1=drift+(int)(*(pointer+ti))*dim;
			for (k=0; k<dim; k++)
				*(pointer1+k)+=*(pointer+k);
			pointer+=ncol;
		}

		// average
		pointer=drift;
		for (t=0; t<total_t; t++)
		{
			for (k=0; k<dim; k++)
				*(pointer+k)/=ptcl_frame[t];
			pointer+=dim;
		}
		free(ptcl_frame);
	}

	//=====================================================================
	// get mean position
	float * mp_ptr = Calloc(float, dim*total_p); POINTER_NULL(mp_ptr);
	float * pointer2;
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		pointer1=mp_ptr+(int)(*(pointer+ii))*dim;
		// remove drift
		if (Remove_Drift)
		{
			pointer2=drift+(int)(*(pointer+ti))*dim;
			for (k=0; k<dim; k++)
				*(pointer1+k)+=*(pointer+k)-*(pointer2+k);
		}
		else
		{
			for (k=0; k<dim; k++)
				*(pointer1+k)+=*(pointer+k);
		}
		pointer+=ncol;
	}
	if (Remove_Drift)
		free(drift);

	// mean position
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
		for (k=0; k<dim; k++)
			*(pointer++)/=frame_ptcl[p];

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
			ymax = * pointer;

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
		csum[i+1]=csum[i]+j;
		sum2+=j*j;
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

	//======================================================
	// coarse-grain box center and particle ids in them
	float xboxsize=(1+1e-6)*(xmax-xmin)/b;
	float yboxsize=(1+1e-6)*(ymax-ymin)/b;

	int * box_index = Malloc(int, total_p); POINTER_NULL(box_index);
	// particle number in boxes
	int * ptcl_box = Calloc(int, b*b); POINTER_NULL(ptcl_box); 
	pointer=mp_ptr;
	for (p=0; p<total_p; p++)
	{
		i=(int)((*pointer-xmin)/xboxsize);
		j=(int)((*(pointer+1)-ymin)/yboxsize);
		k=i*b+j; // located in k'th box
		box_index[p]=k;
		++ptcl_box[k];
		if (box_index[p]>=b*b)
			pERROR("calculating box index wrong.");
		pointer+=dim;
	}

	printf("\n# statistics of particles in a box:\n");
	int bmin=ptcl_box[0], bmax=ptcl_box[0];
	int * Hb= Calloc(int, 100); POINTER_NULL(Hb);
	++Hb[ptcl_box[0]];
	for (i=1; i<b*b; i++)
	{
		p=ptcl_box[i];
		if (bmin > p)
			bmin=p;
		else if (bmax < p)
			bmax=p;
		++Hb[p];
	}
	printf("# min                :\t%d\n"
		   "# max                :\t%d\n", bmin, bmax);
	printf("# histogram:\n");
	for (i=0; i<100; i++)
		if (Hb[i]!=0)
			printf("# \t%d\t%d\n", i, Hb[i]);
	printf("\n");
	free(Hb);
	
	float approxp=total_p/b/b;
	printf("# approximate particle number in a coarse-grained box:\t%.6f\n",
			approxp);

	unsigned char * mask_box = Malloc(unsigned char, b*b);
	POINTER_NULL(mask_box);
	int threshold=(int)(approxp*0.8);
	for (i=0; i<b*b; i++)
	{
		if (ptcl_box[i] < threshold)
			mask_box[i]=0;
		else
			mask_box[i]=1;
	}
	// set mask to false if contain ptcl_no_neighbor, causing nan error.
	//mp_ptr=mp.get_array_pointer();
	for (p=0; p<total_p; p++)
	{
		//p=index_ptcl_no_neighbor[k];
		if (neighbor_count[p]<3)
			mask_box[box_index[p]]=0;
	}
	free(neighbor_count);

	int boxcount=0;
	int * new_box_index = Malloc(int, b*b); POINTER_NULL(new_box_index);
	for (i=0; i<b*b; i++)
		if (mask_box[i])
		{
			new_box_index[i]=boxcount;
			++boxcount;
		}
	
	//==============================================================
	int ncg=dim*(dim+2)+2;
	cg.reserve_memory(ncg, boxcount*total_t);
	float* cg_ptr=cg.get_array_pointer();
	pointer=cg_ptr;
	for (t=0; t<total_t; t++)
	{
		for (i=0; i<b; i++)
		{
			for (j=0; j<b; j++)
			{
				if (mask_box[i*b+j])
				{
					*(pointer++)=xmin+(i+0.5)*xboxsize;
					*(pointer++)=ymin+(j+0.5)*yboxsize;
					for (k=2; k<ncg-1; k++)
						*(pointer++)=0;
					*(pointer++)=t;
				}
			}
		}
	}
	
	//================================================================
	printf("# Least squares\n");
	
	// least square
	// df=(\p f / \p x)dx + (\p f / \p y) dy
	float Xxx, Xxy, Xyy, Yxx, Yxy, Yyx, Yyy;
	float dx0, dy0;
	float delta;

	if (total_t*total_p>2*size[1])
	{
		cout << "# Error: track is bad, too many broken trajectories!\n";
		exit (1);
	}
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
	// coarse-grain

	printf("# coarse-graining\n");
	// approximate particle number in a box
	//float approxp=xboxsize*yboxsize*total_p/(xmax-xmin)/(ymax-ymin);
	
	//cout << "ptid: " << ptid_ptr << ' ' << ptid_ptr+size[2] << endl;
	
	pointer1=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		p=(int)(*(pointer1+ii));
		t=(int)(*(pointer1+ti));
		if (mask_box[box_index[p]])
		{
			//pointer=cg_ptr+(t*b*b+box_index[p])*ncg+2;
			pointer=cg_ptr+(t*boxcount+new_box_index[box_index[p]])*ncg+2;
			
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
			
			*(pointer++)+=exx-1;
			*(pointer++)+=exy;
			*(pointer++)+=eyx;
			*(pointer++)+=eyy-1;
			
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
			*(pointer++)+=chix;
			*(pointer++)+=chiy;
			*(pointer++)+=chix+chiy;
		}

		pointer1+=ncol;
	}

	pointer=cg_ptr;
	//for (p=0;p<b*b*total_p; p++)  important bug
	for (t=0; t<total_t; t++)
		for (p=0;p<b*b; p++)
		{
			if (mask_box[p])
			{
				for (k=2;k<ncg-1; k++)
					*(pointer+k)/=ptcl_box[p];//approxp;
				pointer+=ncg;
			}
		}

	//======================================================

	delete [] size;
	free(frame_ptcl);
	free(mp_ptr);
	free(csum);
	free(neighbor_id);
	free(box_index);
	free(ptcl_box);
	free(mask_box);
	free(new_box_index);
	free(index);
}


void elasticity(const colloid_base& ptid, 
		        const int & b,
				const float & Lambda,
				const int & nbin,
				const char * filename,
				const float & constrain,
				bool Remove_Drift)
{
	colloid_base cg;
	coarse_grain(cg, ptid, b, Lambda, Remove_Drift);
	char *cgfile=getfilename(filename, "_cg.gdf");
	writegdf(cg, cgfile);
	free(cgfile);
	int *size=cg.get_size();
	float * cg_ptr=cg.get_array_pointer();

	printf("# check NaN...\n");
	int i, cc=1;
	for (i=0; i<size[2]; i++)
		if (cg_ptr[i]!=cg_ptr[i])
			printf("# NaN encountered for %d time(s)\n", cc++); 
	if (cc!=1)
		exit (1);
	else
		printf("# No NaN found.\n");

	printf("# statistics\n");
	
	float * dV = Malloc(float, size[1]); POINTER_NULL(dV);
	float * theta = Malloc(float, size[1]); POINTER_NULL(theta);

	cg_ptr+=2;
	for (i=0; i<size[1]; i++)
	{
		// e_xx+e_yy=u_xx+u_yy
		dV[i]=*cg_ptr+*(cg_ptr+3);
		// theta=(u_yx-u_xy)/2
		theta[i]=(*(cg_ptr+2)-*(cg_ptr+1))/2;
		cg_ptr+=size[0];
	}
	cg_ptr=cg.get_array_pointer();

	double MINdV=dV[0], MAXdV=dV[0];
	double MINtheta=theta[0], MAXtheta=theta[0];

	double meandV=dV[0], sddV=dV[0]*dV[0];
	double meantheta=theta[0], sdtheta=theta[0]*theta[0];

	double dVi, thetai;
	for (i=1; i<size[1]; i++)
	{
		dVi=dV[i];
		if ( MINdV > dVi )
			MINdV=dVi;
		else if ( MAXdV < dVi )
			MAXdV=dVi;
		
		thetai=theta[i];
		if ( MINtheta > thetai )
			MINtheta=thetai;
		else if ( MAXtheta < thetai )
			MAXtheta=thetai;

		meandV+=dVi;
		sddV+=dVi*dVi;
		
		meantheta+=thetai;
		sdtheta+=thetai*thetai;

		if (thetai > MAXtheta)
			printf("%f %f\n", thetai, MAXtheta);
	}

	meandV/=size[1];
	sddV=sqrt(sddV/size[1]-meandV*meandV);

	meantheta/=size[1];
	sdtheta=sqrt(sdtheta/size[1]-meantheta*meantheta);

	double binsizedV=(1+1e-6)*(MAXdV-MINdV)/nbin;
	double binsizetheta=(1+1e-6)*(MAXtheta-MINtheta)/nbin;

	// histogram
	printf("histogram\n");
	unsigned long * HdV = Calloc(unsigned long, nbin);
	POINTER_NULL(HdV);
	unsigned long * Htheta = Calloc(unsigned long, nbin);
	POINTER_NULL(Htheta);
	
	for (i=0; i<size[1]; i++)
	{
		//if (dV[i]>=-constrain && dV[i]<constrain)
		++HdV[(int)((dV[i]-MINdV)/binsizedV)];
		if ((int)((theta[i]-MINtheta)/binsizetheta)>=nbin)
			printf("# Error: %d %f %f %d  wrong\n", i, theta[i], MAXtheta,
					(int)((theta[i]-MINtheta)/binsizetheta));
		//if (theta[i]>=-constrain && theta[i] <constrain)
		++Htheta[(int)((theta[i]-MINtheta)/binsizetheta)];
	}

	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo=localtime(	&rawtime );

	FILE * out=fopen(filename, "w");
	FILE_NULL(out, filename);

	fprintf(out, "# COLLOID/ELASTICITY                                   %s\n",
			asctime(timeinfo));

	fprintf(out, "# Independent statistics: \t%d\n", size[1]);

	fprintf(out, 
			"#             Volumn fluctuation\t\tOrientation fluctuation\n#\n");
	fprintf(out, "# mean    :\t%.6f\t\t%.6f\n", meandV, meantheta);
	fprintf(out, "# sd      :\t%.6f\t\t%.6f\n", sddV, sdtheta);
	fprintf(out, "#\n"
				 "# min     :\t%.6f\t\t%.6f\n", MINdV, MINtheta);
	fprintf(out, "# max     :\t%.6f\t\t%.6f\n", MAXdV, MAXtheta);
	fprintf(out, "#\n"
				 "# nbin    :\t%d\t\t\t%d\n", nbin, nbin);
	fprintf(out, "# binsize :\t%.6f\t\t%.6f\n", binsizedV, binsizetheta);
	fprintf(out, "#\n"
				 "# bin_center : probability        "
				 "bin_center : probability\n");
	for (i=0; i<nbin; i++)
		fprintf(out, "%.6f\t%.6f\t%.6f\t%.6f\n",
					MINdV+(i+0.5)*binsizedV,
					(float)(HdV[i])/size[1]/binsizedV, 
					MINtheta+(i+0.5)*binsizetheta,
					(float)(Htheta[i])/size[1]/binsizetheta);
	fclose(out);
	printf("# statistical data is stored in '%s'\n", filename);

	delete [] size;
	free(dV);
	free(theta);
	free(HdV);
	free(Htheta);
}
