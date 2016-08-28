#include "colloid_base.h"
#include "io.h"
#include "data_preprocess.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

void trajectory_filter(colloid_base & cb, const colloid_base & ptid, const int & n)
{
	int i,j,k,c;
	float * ptid_ptr=ptid.get_array_pointer();
	int * size=ptid.get_size();
	
	int total_p=int(ptid_ptr[size[2]-1])+1;
	bool * newindex=new bool [total_p];

	int * frame_ptcl=new int [total_p];
	for (i=0; i<total_p; i++)
		frame_ptcl[i]=0;

	float * pointer1, * pointer2;
	pointer1=ptid_ptr+size[0]-1;
	for (i=0; i<size[1]; i++)
	{
		frame_ptcl[int(*pointer1)]++;
		pointer1+=size[0];
	}

	j=0; c=0;
	for (i=0; i<total_p; i++)
		if (frame_ptcl[i]>=n)
		{
			newindex[i]=true;
			c+=frame_ptcl[i];
			j++;
		}
		else
			newindex[i]=false;
	
	cout << "Particle number with trajectory length larger than " << n
		 << ":\n\t" << j << "\t( total tracked particle:\t"
		 << total_p << " )\n";

	cb.reserve_memory(size[0], c);

	pointer1=ptid_ptr;
	pointer2=cb.get_array_pointer();
	int p, pnew=0;
	for (p=0; p<total_p; p++)
	{
		//p=int(*pointer1+size[0]-1);	
		if (newindex[p]) // is a new index
		{
			for (j=0; j<frame_ptcl[p]; j++)
			{
				for (k=0; k<size[0]-1; k++)
					*(pointer2++)=*(pointer1++);
				*(pointer2++)=pnew;
				pointer1++;
			}
			pnew++;
		}
		else
			pointer1+=frame_ptcl[p]*size[0];
	}


	delete [] size;
	delete [] frame_ptcl;
	delete [] newindex;
}


void trajectory_filter(const char * outfile, const char * infile, const int & n)
{
	colloid_base cb, ptid;
	readgdf(ptid, infile);
	trajectory_filter(cb, ptid, n);
	writegdf(cb, outfile);
}


void trajectory_filter(const char * infile, const int & n)
{
	int c=0;
	while (infile[c]!='\0' && c<100)
		c++;
	if (infile[c-4]!='.' &&
		infile[c-3]!='g' &&
		infile[c-2]!='d' &&
		infile[c-1]!='f')
	{
		cout << "Error: Input filename error, should have suffix .gdf\n";
		exit (0);
	}
	int i=0;
	char *outfile=new char [c+10];
	for (i=0; i<c-4; i++)
		outfile[i]=infile[i];
	
	outfile[i++]='_';
	char tmp [15];
	sprintf(tmp, "%i.gdf", n);
	for (int j=0; j<15; j++)
		if (tmp[j]!='\0')
			outfile[i++]=tmp[j];
		else
			break;
	outfile[i]='\0';
	
	colloid_base cb, ptid;
	readgdf(ptid, infile);
	trajectory_filter(cb, ptid, n);
	writegdf(cb, outfile);
	delete [] outfile;
}


void space_filter(colloid_base & cb, 
		          const colloid_base & ptid,
				  const float & xmin_new,
				  const float & xmax_new,
				  const float & ymin_new,
				  const float & ymax_new,
				  const int & n)
{
	int i,j,k,c;
	float * ptid_ptr=ptid.get_array_pointer();
	int * size=ptid.get_size();
	float * pointer1, * pointer2;
	
	bool * inregion = new bool [size[1]];
	pointer1=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		if (*pointer1 >= xmin_new &&
			*pointer1 <= xmax_new &&
			*(pointer1+1) >= ymin_new &&
			*(pointer1+1) <= ymax_new)
		{
			inregion[i]=true;
		}

		pointer1+=size[0];
	}
	
	int total_p=int(ptid_ptr[size[2]-1])+1;
	bool * newindex=new bool [total_p];

	int * frame_ptcl=new int [total_p];
	for (i=0; i<total_p; i++)
		frame_ptcl[i]=0;

	pointer1=ptid_ptr+size[0]-1;
	for (i=0; i<size[1]; i++)
	{
		if (inregion[i])
			frame_ptcl[int(*pointer1)]++;
		pointer1+=size[0];
	}


	j=0; c=0;
	for (i=0; i<total_p; i++)
		if (frame_ptcl[i]>=n)
		{
			newindex[i]=true;
			c+=frame_ptcl[i];
			j++;
		}
		else
			newindex[i]=false;
	
	cout << "Particle number with trajectory length larger than " << n
		 << "\n and in region\n"
		 << "x: ( " << xmin_new << " , " << xmax_new << " )\n"
		 << "x: ( " << ymin_new << " , " << ymax_new << " )\n"
		 << ":\n\t" << j << "\t( total tracked particle:\t"
		 << total_p << " )\n";

	cb.reserve_memory(size[0], c);

	pointer1=ptid_ptr;
	pointer2=cb.get_array_pointer();
	int p, pnew=0;
	i=0;
	for (p=0; p<total_p; p++)
	{
		//p=int(*pointer1+size[0]-1);	
		if (newindex[p]) // is a new index
		{
			for (j=0; j<frame_ptcl[p]; j++)
			{
				if (inregion[i])
				{
					for (k=0; k<size[0]-1; k++)
						*(pointer2++)=*(pointer1++);
					*(pointer2++)=pnew;
					pointer1++;
					i++;
				}
				else
				{
					i++;
					pointer1+=size[0];
				}
			}
			pnew++;
		}
		else
		{
			i+=frame_ptcl[p];
			pointer1+=frame_ptcl[p]*size[0];
		}
	}


	delete [] size;
	delete [] frame_ptcl;
	delete [] newindex;
}


void pt_remove_drift(colloid_base& pt)
{
	int i, j, k, t;
	int *size=pt.get_size();
	int ncol=size[0];
	int ti=ncol-1; //time index in the array (=dim)
	int dim=ti; // dimension
	float *arr=pt.get_array_pointer();


	int total_time=(int)arr[size[2]-1]+1; // total time
	// = last element of pt+1

	//particle number per frame
	int * ptcl_frame=new int [total_time];
	for (t=0; t<total_time; t++)
		ptcl_frame[t]=0;
	//particle number per frame
	float* pointer=arr+ti;
	for (i=0; i<size[1]; i++)
	{
		ptcl_frame[int(*pointer)]++;
		pointer+=ncol;
	}

	// get drift sum
	t=dim*total_time;
	float *drift_pointer=new float [t];

	for (i=0; i<t; i++)
		drift_pointer[i]=0;
	pointer=arr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+ti))*dim;
		for (k=0; k<dim; k++)
			*(drift_pointer+j+k)+=*(pointer+k);
		pointer+=ncol;
	}
	// average
	pointer=drift_pointer;
	for (t=0; t<total_time; t++)
	{
		for (k=0; k<dim; k++)
			*(pointer+k)/=ptcl_frame[t];
		pointer+=dim;
	}

	// remove drift
	pointer=arr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+ti))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(drift_pointer+j+k);
		pointer+=ncol;
	}
	delete [] size;
	delete [] ptcl_frame;
	delete [] drift_pointer;
}
