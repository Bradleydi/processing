#ifndef MISCELLANEOUS_H
#define MISCELLANEOUS_H

#ifndef vGNUPLOT 
#define vGNUPLOT 4.4
#endif

#ifndef vNEWGNUPLOT 
#define vNEWGNUPLOT 4.0
#endif


#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#define POW3(x) ((x)*(x)*(x))

#define Malloc(type, size) (type*)malloc((size)*sizeof(type))
#define Calloc(type, size) (type*)calloc((size), sizeof(type))

/* test whether the given pointer is NULL */
#define POINTER_NULL(pointer) if (pointer==NULL) { fprintf(stderr, "# Error: \"%s\" line %d\n#        pointer \"%s\" is NULL!\n", __FILE__, __LINE__, #pointer); exit (1);}

/* test whether the given file is open */
#define FILE_NULL(file, filename) if (file==NULL) { fprintf(stderr, "# Error: \"%s\" line %d\n#        cannot open file \"%s\"!\n", __FILE__, __LINE__, filename); exit (1);}

/* wrapper for binary reading and writing */
#define Fread(x, file, filename) if ( (int)fread(&x, sizeof(x), 1, file) != 1 ) { fprintf(stderr, "# Error: \"%s\" line %d\n#        reading '%s' from file '%s' error!\n", __FILE__, __LINE__, #x, filename); exit (1);}

#define Fwrite(x, file, filename) if ( (int)fwrite(&x, sizeof(x), 1, file) != 1 ) { fprintf(stderr, "# Error: \"%s\" line %d\n#        writing '%s' to file '%s' error!\n", __FILE__, __LINE__, #x, filename); exit (1);}

#define FreadN(x, n, file, filename) if ( (int)fread(x, sizeof(x[0]), n, file) != n ) { fprintf(stderr, "# Error: \"%s\" line %d\n#        reading '%s' from file '%s' error!\n", __FILE__, __LINE__, #x, filename); exit (1);}

#define FwriteN(x, n, file, filename) if ( (int)fwrite(x, sizeof(x[0]), n, file) != n ) { fprintf(stderr, "# Error: \"%s\" line %d\n#        writing '%s' to file '%s' error!\n", __FILE__, __LINE__, #x, filename); exit (1);}

#define Fclose(file, filename) if (fclose(file)!=0) { fprintf(stderr, "# Error: \"%s\" line %d\n#        closing file %s failed!\n",  __FILE__, __LINE__, filename); exit (1);}

/*! generate file name by given string and subfix.
 *  @param	file
 *  	type:	const char*
 *  	filename string used to generate a file name
 *  @param	subfix
 *  	type:	const char*
 *  	subfix of the generated file
 *  @return	filename
 *  	type:	char*
 *  	generated filename
 */
char* getfilename(const char*, const char*);


/*! pairwise summation of an array.
 *  This function could be used to avoid the summation of too many 
 *  tiny numbers.
 *  @param	n
 *  	type:	int
 *  	number of elements of the array
 *  @param	array
 *  	type:	double*
 *  	the array
 *  @return sum
 *  	type:	double
 *  	summation of the array
 */
double pairwise(int, double *);


inline void pERROR(const char * message)
{
    (void)fprintf(stderr, "# Error: %s\n", message);
    exit (1);
}


inline unsigned char FILE_EXIST(const char *file)
{
	FILE *f = fopen(file, "r");
	if ( f != NULL )
	{
		fclose (f);
		return 1;
	}
	else
		return 0;
}


/*! Average a two dimensional data set over angles about its center.
 * 	data(x,y) sits at radius r = sqrt( (x-xc)^2 + (y-yc)^2 ) 
 * 	from the center, (xc,yc).  Let R be the integer part
 * 	of r, and dR the fractional part.  Then this point is
 * 	averaged into result(R) with a weight 1-dR and into
 * 	result(R+1) with a weight dR.
 *  @param	nx
 *  	type:	int&
 *  	size in x direction
 *  @param	ny
 *  	type:	int&
 *  	size in y direction
 *  @param	data
 *  	type:	double*
 *  	two dimensional array of any type except string or complex.
 *  	data is saved in the form:
 *  		data[i, j] = data[i+j*nx]
 *  @param	nr
 *  	type:	int&
 *  	size of result array
 *  @return	avg
 *  	type:	double*
 *  	result array averaged over angles as a function of radius from
 *  	the center point, measured in pixels. 
 *  written as a C translation of program written by
 *  	David G. Grier, The University of Chicago, 7/30/92
 */
double * angleavg(int &, int &, double *, int &);

double * points2image(int &, double *, int &, int &, 
					  double * =NULL, bool =true);

void linearfit(int, double *, double &, double &);

void sort(int, double *);
void Merge_sort(double *, double *, int, int);

void sort_int(int, int *);
void Merge_sort_int(int *, int *, int, int);

void getcluster(int Nedge, int *edges, 
		int *pNnode, int *pNcluster, int *cluster);

int* getframeindex(const char *string, int lower, int upper);

void show_matrix(int N, double *A);

/* Select */
double select_NR(int n, double *arr, int k);
double select(int n, double *arr, int k);
double select_wcp(int n, double *arr, double *copy, int k);

/* steady state detection */
void SSD(int N, double *y, int window, int *pNSS, int *SS);
/* batch job for SSD, with workspace allocated before calling
 * workspace = Malloc(char, N)
 */
void SSD_batch(int N, double *y, int window, int *pNSS, int *SS, 
	char *workspace);
/* with statistics, i.e. mean & sigma*/
void SSD_statistics(int N, double *y, int window, int *pNSS, int *SS, 
		double *mean, double *sigma);
void SSD_statistics_batch(int N, double *y, int window, int *pNSS, int *SS, 
		double *mean, double *sigma, char *workspace);
#endif /* MISCELLANEOUS_H */
