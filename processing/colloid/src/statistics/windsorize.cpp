#include <statistics.h>
#include <miscellaneous.h>

/* Thu Nov 22 16:09:03 HKT 2012
 * Here we use Nth_element, (implemented by myself, actuall, C++ also has such
 * a similar program. But to make the program more like C, I use that of myself)
 * I believe that a fast select algorithm should be faster than sort,
 * hence here I don't use sort
 *
 * This program is quite similar to trim. See "statistics/trim.cpp"
 */

/* N: size of the array
 * y: the array to windsorize
 * margin: should be in (0:0.5), to specify the extend of windsorize.
 * 		y[i] < (margin*N)th smallest y should be assigned to the 
 * 			(margin*N)th smallest y
 * 		y[i] > ((1-margin)*N)th smallest y should be assigned to the
 * 			((1-margin)*N)th smallest y
 * Btmargin, Upmargin:
 * 		set the bottom and upper margin separately.
 * return is the resulted array, and also the Bt_id and Up_id for future use.
 * To return id is from the idea of nth_element in c++.
 */

void windsorize(int N, double *y, double margin,
		double *pBt, double *pUp)
{
	if (margin<0.5)
		windsorize_kernel(N, y, margin, 1.0-margin, pBt, pUp);
	else
		windsorize_kernel(N, y, 1.0-margin, margin, pBt, pUp);
}

void windsorize_kernel(int N, double *y, double Btmargin, double Upmargin,
		double *pBt, double *pUp)
{
	if (Btmargin<=0.0 || Upmargin>=1.0 || Btmargin >= Upmargin)
	{
		fprintf(stderr, 
				"# Error: margin should be 0.0 < Btmargin < Upmargin < 1.0\n");
		exit (1);
	}

	/*
	int i, Q1_id, median_id, Q3_id;
	double Q1, median, Q3;
	Nth_element(N, y, N/4, &Q1_id, &Q1);
	Nth_element(N, y, N/2, &median_id, &median);
	Nth_element(N, y, 3*N/4, &Q3_id, &Q3);
	*/
	int i;
	double Bt=select(N, y, (int)((double)N*Btmargin));
	double Up=select(N, y, (int)((double)N*Upmargin));

	for (i=0; i<N; i++)
		if ( y[i] < Bt )
			y[i] = Bt;
		else if ( y[i] > Up )
			y[i] = Up;
	*pBt = Bt; *pUp = Up;
}


/* If to windsorize a lot of arrays, a better way is to use select_wcp
 * to avoid allocating copy arrays in select program
 */
void windsorize_wcp(int N, double *y, double *copy, double margin,
		double *pBt, double *pUp)
{
	if (margin<0.5)
		windsorize_wcp_kernel(N, y, copy, margin, 1.0-margin, pBt, pUp);
	else
		windsorize_wcp_kernel(N, y, copy, 1.0-margin, margin, pBt, pUp);
}

void windsorize_wcp_kernel(int N, double *y, double *copy, 
		double Btmargin, double Upmargin,
		double *pBt, double *pUp)
{
	if (Btmargin<=0.0 || Upmargin>=1.0 || Btmargin >= Upmargin)
	{
		fprintf(stderr, 
				"# Error: margin should be 0.0 < Btmargin < Upmargin < 1.0\n");
		exit (1);
	}

	int i;
	double Bt=select_wcp(N, y, copy, (int)((double)N*Btmargin));
	double Up=select_wcp(N, y, copy, (int)((double)N*Upmargin));

	for (i=0; i<N; i++)
		if ( y[i] < Bt )
			y[i] = Bt;
		else if ( y[i] > Up )
			y[i] = Up;
	*pBt = Bt; *pUp = Up;
	printf("%f %f\n", Bt, Up);
}
