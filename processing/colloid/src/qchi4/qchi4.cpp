#include "colloid_base.h"
#include "qchi4.h"

#include <cstdio>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

void pt_qchi4(colloid_base& pt, double a, int dim)
{
	int i, j, k, t;
	int * size=pt.get_size();
	if (dim!=2)
		dim=size[0]-1;
	else if (dim>size[0]-1)
	{
		cout << "# Error: dim is larger than the given data.\n";
		exit (0);
	}
	else if (dim <0)
	{
		cout << "# Error: dim given is negative.\n";
		exit (0);
	}

	float * pt_ptr=pt.get_array_pointer();

	int total_t=int(pt_ptr[size[2]-1])+1;
	
	
	// remove_drift
	int * ptcl_frame=new int [total_t];
	for (t=0; t<total_t; t++)
		ptcl_frame[t]=0;

	float * pointer=pt_ptr+(size[0]-1); // t
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[int(*pointer)];
		pointer+=size[0];
	}

	int * cum_ptcl_frame=new int [total_t+1];
	cum_ptcl_frame[0]=0;
	for (t=0; t<total_t; t++)
		cum_ptcl_frame[t+1]=cum_ptcl_frame[t]+ptcl_frame[t];

	t=dim*total_t;
	float * drift = new float [t];
	for (i=0; i<t; i++)
		drift[i]=0;

	pointer=pt_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+(size[0]-1)))*dim;
		for (k=0; k<dim; k++)
			*(drift+(j+k))+=*(pointer+k);
		pointer+=size[0];
	}
	// average
	pointer=drift;
	for (t=0; t<total_t; t++)
	{
		for (k=0; k<dim; k++)
			*(pointer++)/=ptcl_frame[t];
	}


	// remove drift
	//cout << "remove drift\n";
	pointer=pt_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+(size[0]-1)))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(drift+(j+k));
		pointer+=size[0];
	}

	delete [] drift;

	// get V
	float * min =new float [dim];
	float * max =new float [dim];
	pointer=pt_ptr;
	for (i=0; i<size[1]; i++)
	{
		for (k=0; k<dim; k++)
		{
			if (min[k] > *(pointer+k))
				min[k]=*(pointer+k);
			else if (max[k] < *(pointer+k))
				max[k]=*(pointer+k);
		}
		pointer+=size[0];
	}
	float V=1;
	for (k=0; k<dim; k++)
		V*=(max[k]-min[k]);

	//=============================================================
	// divided into box and get box id
	unsigned int nxbox=(int)((max[0]-min[0])/a)+1;
	unsigned int *xbox_id=NULL;
	unsigned int *ybox_id=NULL;
	unsigned int *count_xbox=NULL;
	unsigned int *cum_xbox=NULL;
	int *xbox_points=NULL;
	int ti=size[0]-1;
	
	if (dim==2)
	{
		xbox_id=new unsigned int [size[1]];
		ybox_id=new unsigned int [size[1]];
		count_xbox=(unsigned int*)calloc(nxbox*total_t,sizeof(unsigned int));
		pointer=pt_ptr;
		for (i=0; i<size[1]; i++)
		{
			xbox_id[i]=(int)((*pointer-min[0])/a);
			t=(int)(*(pointer+ti));
			++count_xbox[nxbox*t+xbox_id[i]];
			ybox_id[i]=(int)((*(pointer+1)-min[1])/a);
			pointer+=size[0];
		}
		cum_xbox=new unsigned int [nxbox*total_t+1];
		cum_xbox[0]=0;
		for (i=0; i<nxbox*total_t; i++)
		{
			cum_xbox[i+1]=cum_xbox[i]+count_xbox[i];
			count_xbox[i]=0;
		}
		if (cum_xbox[nxbox*total_t] != size[1])
		{
			printf("# Error: cum_xbox!\n");
			exit(1);
		}

		xbox_points=new int [size[1]];
		pointer=pt_ptr;
		for (i=0; i<size[1]; i++)
		{
			j=(int)((*pointer-min[0])/a);
			t=(int)(*(pointer+ti));
			j+=t*nxbox;
			xbox_points[cum_xbox[j]+(count_xbox[j]++)]=i;
			pointer+=size[0];
		}

		free(count_xbox);
	}

	delete [] min;
	delete [] max;


	//=============================================================
	//=============================================================
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo=localtime( &rawtime );

	// output to stdout
	//cout.setf(ios::floatfield);
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(6);
	//cout.width(12);

	cout << "# COLLOID/DYNAMICS/QCHI4                                "
		 << asctime (timeinfo) << endl;

	double ave_p=double(size[1])/total_t;
	unsigned long long sd=0;
	for (i=0; i<total_t; i++)
		sd+=ptcl_frame[i]*ptcl_frame[i];
	cout << "# Average particles           :\t\t" << ave_p 
		 << "    +/-    " << sqrt(double(sd)/total_t-ave_p*ave_p) << endl;
	cout << "# Frame number                :\t\t" << total_t << endl;

	cout << "# Overlap cutoff a            :\t\t" << a << endl;


	cout << "\n# dt\tq\tchi4\tindependent_statistics";
	cout << "\n# ==\t=\t====\t======================\n";
	//============================================================


	// get dts
	int dts [150];
	dts[0]=1;
	int cdt=1;
	for (i=0; i<150; i++)
	{
		j=int(pow(1.15,i)+0.5);
		if (j>total_t-100) // should more than 1000 independent statistics
			break;
		if (j==dts[cdt-1])
			continue;
		dts[cdt++]=j;
	}


	// get qchi4
	//cout << "get dx\n";
	float * pointer1, *pointer2;
	int count, tt;
	//double q, chi4;
	unsigned long long q, q2, qt;
	double d, dx;
	int t1, t2;
	double aa=a*a;
	double ave;
	unsigned int ixbox, iybox, jbox, low, up, dij;

	if (dim==2)
	{
		for (tt=0; tt<cdt; tt++) // dts
		{
			q=0ULL; q2=0ULL;
			count=total_t-dts[tt];
			for (t1=0; t1<count; t1++)
			{
				qt=0ULL;
				t2=t1+dts[tt];
				pointer1=pt_ptr+(cum_ptcl_frame[t1]*size[0]);
				for (i=cum_ptcl_frame[t1]; i< cum_ptcl_frame[t1+1]; i++)
				{
					ixbox=xbox_id[i];
					iybox=ybox_id[i];
					
					low=cum_xbox[ixbox-1+t2*nxbox];
					up=cum_xbox[ixbox+2+t2*nxbox];

					for (k=low; k<up; k++)
					{	
						j=xbox_points[k];
						jbox=ybox_id[j];
						dij = (iybox>=jbox) ? (iybox-jbox) : (jbox-iybox);
						if (dij<=1)
						{
							pointer2=pt_ptr+(j*size[0]);
							dx=*pointer1-*pointer2;
							d=dx*dx;
							if (d<=aa)
							{
								dx=*(pointer1+1)-*(pointer2+1);
								d+=dx*dx;
								if (d<=aa)
									++qt;
							}
						}
					}
					pointer1+=size[0];
				}
				q+=qt;
				q2+=qt*qt;
			}
			
			// output the msd
			ave=double(q)/count;
			//cout << dts[tt] << '\t' << ave << 
			//	 '\t' << double(q2)/count-ave*ave << endl;
			cout << dts[tt] << '\t' << ave/ave_p << '\t'
				 << (double(q2)/count-ave*ave)*V/ave_p/ave_p 
				 << '\t' << count << endl;
		}

		delete [] xbox_id;
		delete [] ybox_id;
		delete [] cum_xbox;
		delete [] xbox_points;
	}
	else 
	{
		for (tt=0; tt<cdt; tt++) // dts
		{
			q=0ULL; q2=0ULL;
			count=total_t-dts[tt];
			for (t1=0; t1<count; t1++)
			{
				qt=0ULL;
				t2=t1+dts[tt];
				pointer1=pt_ptr+(cum_ptcl_frame[t1]*size[0]);
				for (i=cum_ptcl_frame[t1]; i< cum_ptcl_frame[t1+1]; i++)
				{
					pointer2=pt_ptr+(cum_ptcl_frame[t2]*size[0]);
					
					for (j=cum_ptcl_frame[t2]; j<cum_ptcl_frame[t2+1]; j++)
					{	
						dx=*pointer1-*pointer2;
						d=dx*dx;
						for (k=1; k<dim; k++)
						{
							dx=*(pointer1+k)-*(pointer2+k);
							d+=dx*dx;
							//if (d>aa)
							//	break;
						}
						if (d<=aa)
							++qt;
						pointer2+=size[0];
					}
					pointer1+=size[0];
				}
				q+=qt;
				q2+=qt*qt;
			}
			
			// output the msd
			ave=double(q)/count;
			cout << dts[tt] << '\t' << ave/ave_p << '\t'
				 << (double(q2)/count-ave*ave)*V/ave_p/ave_p 
				 << '\t' << count << endl;
		}
	}

	delete [] ptcl_frame;
	delete [] cum_ptcl_frame;
	delete [] size;
}
