#include "colloid_base.h"
#include "io.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

void show_file_size(unsigned int filesize)
{
	char str[10];
	if (filesize<1024)
		sprintf(str, "%dB", filesize);
	else if (filesize>=1024 && filesize<1024*1024)
		sprintf(str, "%.2fKB", (float)filesize/1024);
	else if (filesize>=1024*1024 && filesize<1024*1024*1024)
		sprintf(str, "%.2fMB", (float)filesize/1024/1024);
	else
		sprintf(str, "%.2fGB", (float)filesize/1024/1024/1024);

	cout << str;
}

void readgdf(colloid_base & cb, const char * filename)
{
	ifstream in(filename, ios::in | ios::binary | ios::ate);
	if ( ! in.is_open() )
	{
		cout << "# Error: File can not open!\n";
		exit (1);
	}

	ifstream::pos_type filesize;
	filesize=in.tellg();
	cout << "# Input data file name:\t\t" << filename << endl;
	cout << "# Input data file size:\t\t";
	show_file_size((unsigned int)filesize);
	cout << endl;
	
	int header;
	in.seekg(0,ios::beg);
	in.read((char*)&header, 4);
	if (header !=82991)
	{
		cout << "# Error: Input file is not a .gdf file.\n";
		in.close();
		exit (1);
	}

	int size [5];
	in.read((char*)size, 20);
	int s1=size[1], s2=size[2], s3=size[4];

	cb.reserve_memory(size[1], size[2]);
	float *arr=cb.get_array_pointer();
	in.read((char*)arr, sizeof(float)*s3);
	if (in.tellg()!=filesize)
	{
		cout << "# Error: Input .gdf file may wrong during reading the array.\n";
		exit (1);
	}
	in.close();
	cout << "# Array: " << s1 << "x" << s2 << " (="
		 << s3 << ")\n";

}


void writegdf(colloid_base & cb, const char * filename)
{
	ofstream out (filename, ios::out | ios::binary | ios::trunc);
	if ( ! out.is_open() )
	{
		cout << "# Error: File can not open!\n";
		exit (1);
	}
	
	int headersize [6];
	headersize[0]=82991;
	int *size=cb.get_size();
	headersize[1]=2;
	headersize[2]=size[0];
	headersize[3]=size[1];
	headersize[4]=4; // IDL float type code
	headersize[5]=size[2];
	out.write((char*)headersize, 24); // write header and size
	out.write((char*)cb.get_array_pointer(), sizeof(float)*size[2]);
	cout << "# Totally ";
	show_file_size(sizeof(float)*size[2]+24);
	cout << " disk space is used to store the gdf data.\n";
	cout << "# Array: " << size[0] << "x" << size[1] << " (="
		 << size[2] << ")\n";
	out.close();
	delete [] size;
}
