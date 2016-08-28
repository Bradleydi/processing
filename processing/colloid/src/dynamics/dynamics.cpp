#include "dynamics.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

void dxpdf(colloid_base& ptid, 
		   const int & dt, const int & dim, const int & nbin)
{
	int i, j, k, t, p;
	int * size=ptid.get_size();
	int tindex=size[0]-2;
	if (dim<1 || dim>tindex)
	{
		cout << "# Error: wrong dimension!\n"
	         << "# \tDimension can only be 1 to " << tindex << endl;
		exit (1);
	}
	float * ptid_ptr=ptid.get_array_pointer();


	// get total_t & total_t
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	float* pointer=ptid_ptr+tindex;
	float T=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if ( T < *pointer )
			T = *pointer;
	}
	const int total_t=(int)(T)+1;


	// remove_drift
	int * ptcl_frame=new int [total_t];
	for (t=0; t<total_t; t++)
		ptcl_frame[t]=0;

	pointer=ptid_ptr+tindex;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[int(*pointer)];
		pointer+=size[0];
	}

	t=dim*total_t;
	float * drift = new float [t];
	for (i=0; i<t; i++)
		drift[i]=0;

	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+tindex))*dim;
		for (k=0; k<dim; k++)
			*(drift+j+k)+=*(pointer+k);
		pointer+=size[0];
	}
	// average
	pointer=drift;
	for (t=0; t<total_t; t++)
	{
		for (k=0; k<dim; k++)
			*(pointer++)/=ptcl_frame[t];
	}

	delete [] ptcl_frame;

	// remove drift
	//cout << "remove drift\n";
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+tindex))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(drift+j+k);
		pointer+=size[0];
	}

	delete [] drift;

	// get index
	j=total_t*total_p;
	int * index=new int [j];
	for (i=0; i<j; i++)
		index[i]=-1;
	pointer=ptid_ptr+tindex;
	for (i=0; i<size[1]; i++)
	{
		index[int(*pointer)+int(*(pointer+1))*total_t]=i*size[0];
		pointer+=size[0];
	}

	// get dx
	//cout << "get dx\n";
	j=total_t-dt;
	float * dx =new float [dim*j*total_p];
	float * pointer2, * dxptr=dx;
	int count=0;
	float * min = new float [dim];
	float * max = new float [dim];
	double * mean = new double [dim];
	double * dx2 = new double [dim];
	double * dx4 = new double [dim];
	float d;
	for (i=0; i<dim; i++)
	{
		min[i]=0; max[i]=0;
		mean[i]=0;
		dx2[i]=0; dx4[i]=0;
	}
	for (p=0; p<total_p; p++)
		for (t=0; t<j; t++)
		{
			i=index[t+p*total_t];
			if (i!=-1)
			{
				pointer=ptid_ptr+i;
				i=index[t+dt+p*total_t];
				if (i!=-1)
				{
					pointer2=ptid_ptr+i;
					for (k=0; k<dim; k++)
					{
						d=*(pointer2++)-*(pointer++);
						if ( min[k] > d )
							min[k] = d;
						else if ( max[k] < d)
							max[k] =d;
						*(dxptr++)=d;
						mean[k]+=d;
						d=d*d;
						dx2[k]+=d;
						dx4[k]+=d*d;
					}
					++count;
				}
			}
		}

	// statistics
	//cout << "statistics\n";
	float * binsize = new float [dim];
	for (k=0; k<dim; k++)
		binsize[k]=(max[k]-min[k])*(1+1e-6)/nbin;

	j=dim*nbin;
	int * hist= new int [j];
	for (i=0; i<j; i++)
		hist[i]=0;

	dxptr=dx;
	for (i=0; i<count; i++)
	{
		for (k=0; k<dim; k++)
		{
			++hist[dim*int((*(dxptr++)-min[k])/binsize[k])+k];
		}
	}

	// statsitcs for total dx
	unsigned int * thist = new unsigned int [nbin];
	for (i=0; i<nbin; i++)
		thist[i]=0;
	float tmin=min[0], tmax=max[0];
	double tmean=0, tdx2=0, tdx4=0;
	for (k=1; k<dim; k++)
	{
		tmin=( tmin <= min[k] ) ? tmin : min[k];
		tmax=( tmax >= max[k] ) ? tmax : max[k];
	}
	float tbinsize=(tmax-tmin)*(1+1e-6)/nbin;
	dxptr=dx;
	j=count*dim;
	for (i=0; i<j; i++)
	{
		d=*(dxptr++);
		thist[int((d-tmin)/tbinsize)]++;
		tmean+=d;
		d=d*d;
		tdx2+=d;
		tdx4+=d*d;
	}

	delete [] dx;

	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo=localtime( &rawtime );

	// output to stdout
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(6);

	cout << "# COLLOID/DYNAMICS/DXPDF                                "
		 << asctime (timeinfo) << endl;

	cout << "# Particle number     :\t\t" << total_p << endl;
	cout << "# Frame number        :\t\t" << total_t << endl;

	cout << "\n# dt                  :\t\t" << dt;
	//cout << "\n#      \tmean   \tdx2   \tdx4   \talpha2\n" ;
	//for (k=0; k<dim; k++)
	//{
	//	cout << "# dim " << k << '\t'
	//		 << mean[k]/count << '\t' 
	//		 << dx2[k]/count << '\t' 
	//		 << dx4[k]/count << '\t'
	//		 << dx4[k]*count/dx2[k]/dx2[k]/3-1 << endl;
	//}

	cout << "\n# statistics of displacement";
	cout << "\n# dim   ";
	for (k=0; k<dim; k++)
		cout << "      \t" << k << "      ";
	cout << "      \ttotal";
	cout << "\n# ======";
	for (k=0; k<dim; k++)
		cout << "======\t=======";
	cout << "======\t=====";
	cout << "\n# min   ";
	for (k=0; k<dim; k++)
		cout << "      \t" << min[k];
	cout << "      \t" << tmin;
	cout << "\n# max   ";
	for (k=0; k<dim; k++)
		cout << "      \t" << max[k];
	cout << "      \t" << tmax;
	cout << "\n# mean  ";
	for (k=0; k<dim; k++)
		cout << "      \t" << mean[k]/count;
	cout << "      \t" << tmean/count/dim;
	cout << "\n# dx2   ";
	for (k=0; k<dim; k++)
		cout << "      \t" << dx2[k]/count;
	cout << "      \t" << tdx2/count/dim;
	cout << "\n# dx4   ";
	for (k=0; k<dim; k++)
		cout << "      \t" << dx4[k]/count;
	cout << "      \t" << tdx4/count/dim;
	cout << "\n# alpha2";
	for (k=0; k<dim; k++)
		cout << "      \t" << dx4[k]*count/dx2[k]/dx2[k]/3-1;
	cout << "      \t" << tdx4*count*dim/tdx2/tdx2/3-1;


	cout << "\n\n# nbin   ";
	for (k=0; k<=dim; k++)
		cout << "      \t" << nbin;
	cout << "\n# binsize";
	for (k=0; k<dim; k++)
		cout << "      \t" << binsize[k];
	cout << "      \t" << tbinsize;
	

	cout << "\n\n# ";
	for (k=0; k<dim; k++)
		cout << "dx(dim" << k << ") probability ";
	cout << "dx(total) probabilit\n";
	for (i=0; i<nbin; i++)
	{
		for (k=0; k<dim; k++)
			cout << (i+0.5)*binsize[k]+min[k]
				 << ' ' << float(hist[i*dim+k])/count/binsize[k] << ' ';
		cout << (i+0.5)*tbinsize+tmin << ' '
			 << float(thist[i])/count/dim/tbinsize << ' ';
		cout << endl;
	}

	delete [] index;
	delete [] min;
	delete [] max;
	delete [] mean;
	delete [] dx2;
	delete [] dx4;
	delete [] hist;
	delete [] thist;
	delete [] size;
	delete [] binsize;
}


void msd_i(colloid_base& ptid, const int& dim)
{
	int i, j, k, t, p;
	int * size=ptid.get_size();
	int tindex=size[0]-2;
	if (dim<1 || dim>tindex)
	{
		cout << "# Error: wrong dimension!\n"
	         << "# \tDimension can only be 1 to " << tindex << endl;
		exit (1);
	}
	float * ptid_ptr=ptid.get_array_pointer();

	
	// get total_t & total_t
	const int total_p=(int)(ptid_ptr[size[2]-1])+1;
	
	float* pointer=ptid_ptr+tindex;
	float T=*pointer;
	for (i=1; i<size[1]; i++)
	{
		pointer+=size[0];
		if ( T < *pointer )
			T = *pointer;
	}
	const int total_t=(int)(T)+1;

	
	// remove_drift
	int * ptcl_frame=new int [total_t];
	for (t=0; t<total_t; t++)
		ptcl_frame[t]=0;

	pointer=ptid_ptr+tindex;
	for (i=0; i<size[1]; i++)
	{
		++ptcl_frame[int(*pointer)];
		pointer+=size[0];
	}

	t=dim*total_t;
	float * drift = new float [t];
	for (i=0; i<t; i++)
		drift[i]=0;

	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+tindex))*dim;
		for (k=0; k<dim; k++)
			*(drift+j+k)+=*(pointer+k);
		pointer+=size[0];
	}
	// average
	pointer=drift;
	for (t=0; t<total_t; t++)
	{
		for (k=0; k<dim; k++)
			*(pointer++)/=ptcl_frame[t];
	}

	delete [] ptcl_frame;

	// remove drift
	//cout << "remove drift\n";
	pointer=ptid_ptr;
	for (i=0; i<size[1]; i++)
	{
		j=int(*(pointer+tindex))*dim;
		for (k=0; k<dim; k++)
			*(pointer+k)-=*(drift+j+k);
		pointer+=size[0];
	}

	delete [] drift;

	// get index
	j=total_t*total_p;
	int * index=new int [j];
	for (i=0; i<j; i++)
		index[i]=-1;
	pointer=ptid_ptr+tindex;
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

	// output to stdout
	cout.setf(ios::fixed, ios::floatfield);
	cout.precision(6);

	cout << "# COLLOID/DYNAMICS/MSD                                  "
		 << asctime (timeinfo) << endl;

	cout << "# Particle number     :\t\t" << total_p << endl;
	cout << "# Frame number        :\t\t" << total_t << endl;


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
	//double * mean = new double [dim];
	double * dx2 = new double [dim];
	double d;

	for (tt=0; tt<cdt; tt++) // dts
	{
		count=0;
		j=total_t-dts[tt];
		for (i=0; i<dim; i++) // initialized
		{
			//mean[i]=0;
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
							//mean[k]+=d;
							dx2[k]+=d*d;
						}
						++count;
					}
				}
			}
		
		// output the msd
		cout << dts[tt] << '\t'; 
		/*
		for (k=0; k<dim; k++)
		{
			mean[k]/=count;
			cout << mean[k] << '\t';
		}
		*/
		d=0;
		for (k=0; k<dim; k++)
		{
			dx2[k]/=count;
			d+=dx2[k];
			cout << dx2[k] << '\t';
		}
		cout << d << '\t' <<  count << endl;

	}


	delete [] index;
	//delete [] mean;
	delete [] dx2;
	delete [] size;
}


void qchi4(colloid_base & pt, float& a)
{
	
}
