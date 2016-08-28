#include "dynamics.h"

#include "colloid_base.h"
#include "data_preprocess.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

void msd(colloid_base& ptid, const int& dim, bool Remove_drift)
{
	int i, j, k, t, p;
	int * size=ptid.get_size();
	const int ti=size[0]-2;
	if (dim<1 || dim>ti)
	{
		fprintf(stderr, "# Error: wrong dimension!\n"
	         			"# \tDimension can only be 1 to %d\n", ti);
		exit (1);
	}
	float * ptid_ptr=ptid.get_array_pointer();

	
	// get total_p & total_t
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	float* pointer=ptid_ptr+ti;
	float T=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if ( T < *pointer )
			T = *pointer;
	}
	const int total_t=(int)(T)+1;

	
	// remove_drift
	if (Remove_drift)
		remove_drift(ptid, dim);

	// get index
	j=total_t*total_p;
	int * index=new int [j];
	POINTER_NULL(index);
	for (i=0; i<j; i++)
		index[i]=-1;
	pointer=ptid_ptr+ti;
	for (i=0; i<size[1]; i++)
	{
		index[int(*pointer)+int(*(pointer+1))*total_t]=i*size[0];
		pointer+=size[0];
	}


	//=============================================================
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo=localtime( &rawtime );


	printf("# COLLOID/DYNAMICS/MSD                                  %s\n",
		 asctime (timeinfo) );

	printf("# Particle number     :        %d\n", total_p);
	printf("# Frame number        :        %d\n", total_t);


	cout << "\n# MSD (mean squared displacements)";
	cout << "\n#                    dim";
	cout << "\n#       ";
	for (k=0; k<dim; k++)
		cout <<	"===============";
	cout << "\n# dt    ";
	for (k=0; k<dim; k++)
		cout << k << "               ";
	cout << "total        independent statistics";
	cout << "\n# ===   ";
	for (k=0; k<dim; k++)
		cout << "===             ";
	cout << "=====        =====================\n";
	//============================================================


	// get dts
	int dts [150];
	dts[0]=1;
	int cdt=1;
	for (i=0; i<150; i++)
	{
		j=int(pow(1.15,i)+0.5);
		if (j>total_t) 
			break;
		if (j==dts[cdt-1])
			continue;
		dts[cdt++]=j;
	}


	// get dx
	//cout << "get dx\n";
	float * pointer2;
	int count, tt;
	double * mean = new double [dim];
	double * dx2 = new double [dim];
	POINTER_NULL(dx2);
	double d;

	for (tt=0; tt<cdt; tt++) // dts
	{
		count=0;
		j=total_t-dts[tt];
		for (i=0; i<dim; i++) // initialized
		{
			mean[i]=0;
			dx2[i]=0;
		}
		for (p=0; p<total_p; p++) // each particle
			for (t=0; t<j; t++)   // each time
			{
				i=index[t+p*total_t];
				if (i!=-1)
				{
					pointer=ptid_ptr+i;
					i=index[t+dts[tt]+p*total_t];
					if (i!=-1)
					{
						pointer2=ptid_ptr+i;
						for (k=0; k<dim; k++)
						{
							d=*(pointer2++)-*(pointer++);
							mean[k]+=d;
							dx2[k]+=d*d;
						}
						++count;
					}
				}
			}
		
		// output the msd
		cout << dts[tt] << '\t'; 
		
		for (k=0; k<dim; k++)
		{
			mean[k]/=count;
			printf("%6.6f\t", mean[k]);
		}
		
		d=0;
		for (k=0; k<dim; k++)
		{
			dx2[k]=dx2[k]/count-mean[k]*mean[k];
			d+=dx2[k];
			printf("%6.6f\t", dx2[k]);
		}
		printf("%6.6f\t%d\n", d, count);

	}


	delete [] index;
	delete [] mean;
	delete [] dx2;
	delete [] size;
}
