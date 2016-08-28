#include "colloid_base.h"

#include <iostream>
#include <cstdlib>

using namespace std;

colloid_base::colloid_base() {}

colloid_base::colloid_base(int col, int row)
{
	ncol=col; row=nrow; total=col*row;
	array=new float [total];
}


colloid_base::colloid_base(int col, int row, float value)
{
	ncol=col; nrow=row; total=col*row;
	array=new float [total];
	for (int i=0; i<total; i++)
		array[i]=value;
}

colloid_base::colloid_base(int col, int row, float* values)
{
	ncol=col; nrow=row; total=col*row;
	array=new float [total];
	for (int i=0; i<total; i++)
		array[i]=values[i];
}

colloid_base::colloid_base(const colloid_base& cp)
{
	ncol=cp.ncol;
	nrow=cp.nrow;
	total=ncol*nrow;
	array=new float [total];
	float * cparray=cp.array;
	for (int i=0; i<total; i++)
		array[i]=cparray[i];
}

colloid_base::~colloid_base()
{
	//cout << array << "destroy\n";
	delete [] array;
}

void colloid_base::operator=(const colloid_base& cp)
{
	//colloid_base(cp); will shadow parameter
	ncol=cp.ncol;
	nrow=cp.nrow;
	total=ncol*nrow;
	array=new float [total];
	float * cparray=cp.array;
	for (int i=0; i<total; i++)
		array[i]=cparray[i];
}

colloid_base colloid_base::operator+(const colloid_base& add)
{
	if (ncol!=add.ncol || nrow!=add.nrow)
	{
		cout << "SizeError: two array have different sizes.";
		exit (1);
	}

	float * a=new float [total];
	for (int i=0; i<total; i++)
		a[i]=array[i]+add.array[i];
	return colloid_base(ncol, nrow, a);
}

colloid_base colloid_base::operator-(const colloid_base& minus)
{
	if (ncol!=minus.ncol || nrow!=minus.nrow)
	{
		cout << "SizeError: two array have different sizes.";
		exit (1);
	}

	float * a=new float [total];
	for (int i=0; i<total; i++)
		a[i]=array[i]-minus.array[i];
	return colloid_base(ncol, nrow, a);
}

void colloid_base::set_value(int nc, int nr, float value)
{
	ncol=nc; nrow=nr; total=nc*nr;
	array=new float [total];
	for (int i=0; i<total; i++)
		array[i]=value;
}

void colloid_base::set_value(int nc, int nr, float* arr)
{
	ncol=nc; nrow=nr; total=nc*nr;
	array=new float [total];
	for (int i=0; i<total; i++)
		array[i]=arr[i];
}

void colloid_base::reserve_memory(int nc, int nr)
{
	ncol=nc; nrow=nr; total=nc*nr;
	array=new float [total];
}

float colloid_base::v(int icol, int irow)
{
	return array[icol+irow*ncol];
}

void colloid_base::get_column(float * store, int icol)
{
	store=new float [ncol];
	for (int i=0; i<ncol; i++)
		store[i]=array[icol+i*ncol];
}

void colloid_base::get_row(float * store, int irow)
{
	store=new float [nrow];
	for (int i=0; i<nrow; i++)
		store[i]=array[i+irow*ncol];
}

//void colloid_base::get_size(int** size)
int* colloid_base::get_size() const
{
	int *size=new int [3];
	size[0]=ncol;
	size[1]=nrow;
	size[2]=total;
	return size;
}

//void colloid_base::show_size()
//{
//	cout << "Array: " << ncol << 'x' << nrow
//		 << " (=" << total << ")\n";
//}

// put inline here, not in prototype
// inline function should be put in header file
//inline float * colloid_base::get_array_pointer() 
//{
//	return array;
//}
