#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

/*
 * N: data size
 * y: data (no x specified)
 * w: window width
 * pNSS: pointer to Number of Steady States
 * SS: pointer to the beginning and ending of Steady States.
 * so length of the i'th steady state is
 * 		SS[2*i+1]-SS[2*i]+1
 * 	SS should be allocated before calling the function.
 * 	One can estimate the size, but N/2*sizeof(int) is always more than enough
 */
void SSD(int N, double *y, int window, int *pNSS, int *SS)
{
	if ( window <= 0 )
	{
		window = N/10;
		fprintf(stderr, "# set window width to %d\n", window);
	}
	else if ( window >= N )
	{
		fprintf(stderr, "# Error: Window width >= sample number\n");
		exit (1);
	}

	int i, j, lastSteady=0;
	double sum=y[0], sum2=y[0]*y[0], sum2ad=0.0;
	double dbW=(double)window;
	// moving window
	for (i=1; i<window; i++)
	{
		//mean += y[i];
		sum += y[i];
		sum2 += y[i]*y[i];
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1]);
	}
	char *steady=Calloc(char, N); POINTER_NULL(steady);
	if (sum2ad > (sum2-sum*sum/dbW))
	{
		for (i=0; i<window; i++)
			steady[i]=1;
		lastSteady = window-1;
	}

	// if I use moving overlaped window, then the resolution will be 1
	// for steady state intervals.
	// If partially overlaped, the resolution will be lower, and for
	// nonoverlapped window, the resolution is the window width, very low.
	for (i=window; i<N; i++) // i: end of the window
	{
		sum += (y[i]-y[i-window]);	
		sum2 += (y[i]*y[i] - y[i-window]*y[i-window]);
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1])
			-(y[i-window+1]-y[i-window])*(y[i-window+1]-y[i-window]);
		if (sum2ad > (sum2-sum*sum/dbW))
		{
			if ( lastSteady + window < i) // lastSteady is small
				j=i-window;
			else
				j=lastSteady;
			while (j<=i)
				steady[j++]=1;
			lastSteady = i;
		}
	}
	/*for (i=0; i<N; i++)
		printf("%d %d\n", i, (int)steady[i]);
	fprintf(stderr, "OK\n");*/

	/*===================================================================
	 * steady state dectection
	 *===================================================================*/
	int SSs=0, SSe, imin=0, NSS=0, SSlength=0;
	// SSs : starting point of a steady state
	// SSe : end point of a steady state
	// imin: from imin to detect the next steady state
	const int SSlengthMin=10; 
	// at least 10 samples can be viewed as steady state 
	char recording;
	while (1)
	{
		// detect steady state
		recording=0;
		for (i=imin; i<N; i++)
		{
			if (steady[i])
			{
				if (recording) // SS is recording
					++SSlength;
				else
				{	SSs=i; SSlength=1; recording=1; }
			}
			else if (recording) // SS is recording
			{
				if (SSlength>=SSlengthMin) // if length is enough
					break;
				else // length is too short to be a steady state
					recording=0;
			}
		}
		++NSS; // a SS is successfully found
		SSe=SSs+SSlength-1;
		//fprintf(stderr, "# %dth steady state %d:%d\n", NSS, SSs, SSe);

		/*
		// find mean and var
		sum=y[SSs]; sum2=y[SSs]*y[SSs];
		for (j=SSs+1; j<=SSe; j++)  
			// cannot use i, since i is used to detect end of the while loop
		{
			sum += y[j];
			sum2 += y[j]*y[j];
		}
		mean = sum/SSlength;
		sigma = sqrt(sum2/SSlength - mean*mean);
		//fprintf(stderr, "# mean=%f simga=%f\n", mean, sigma);
		//fprintf(stderr, "# max possible trend = %f\n", 
		//		2.0*sqrt(3.0)*sigma/dbW);
		*/
		SS[2*(NSS-1)]=SSs;
		SS[2*NSS-1]=SSe;
		imin=SSe+1;
		
		if (i==N) 
		{//fprintf(stderr, "# Only %d steady state(s)\n", NSS); 
			break;}
	}
	*pNSS=NSS;

	free(steady);
}

/* batch job for SSD, with workspace allocated before calling
 * workspace = Malloc(char, N)
 */
void SSD_batch(int N, double *y, int window, int *pNSS, int *SS, 
	char *workspace)
{
	if ( window <= 0 )
	{
		window = N/10;
		fprintf(stderr, "# set window width to %d\n", window);
	}
	else if ( window >= N )
	{
		fprintf(stderr, "# Error: Window width >= sample number\n");
		exit (1);
	}

	int i, j, lastSteady=0;
	double sum=y[0], sum2=y[0]*y[0], sum2ad=0.0, r=7.0/6.0;
	double dbW=(double)window;
	// moving window
	for (i=1; i<window; i++)
	{
		//mean += y[i];
		sum += y[i];
		sum2 += y[i]*y[i];
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1]);
	}
	char *steady=workspace;
	for (i=0; i < N; i++) steady[i]=0;
	if ( r*sum2ad > (sum2-sum*sum/dbW))
	{
		for (i=0; i<window; i++)
			steady[i]=1;
		lastSteady = window-1;
	}

	// if I use moving overlaped window, then the resolution will be 1
	// for steady state intervals.
	// If partially overlaped, the resolution will be lower, and for
	// nonoverlapped window, the resolution is the window width, very low.
	for (i=window; i<N; i++) // i: end of the window
	{
		sum += (y[i]-y[i-window]);	
		sum2 += (y[i]*y[i] - y[i-window]*y[i-window]);
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1])
			-(y[i-window+1]-y[i-window])*(y[i-window+1]-y[i-window]);
		if ( r*sum2ad > (sum2-sum*sum/dbW))
		{
			if ( lastSteady + window < i) // lastSteady is small
				j=i-window;
			else
				j=lastSteady;
			while (j<=i)
				steady[j++]=1;
			lastSteady = i;
		}
	}
	/*for (i=0; i<N; i++)
		printf("%d %d\n", i, (int)steady[i]);
	fprintf(stderr, "OK\n");*/

	/*===================================================================
	 * steady state dectection
	 *===================================================================*/
	int SSs=0, SSe, imin=0, NSS=0, SSlength=0;
	// SSs : starting point of a steady state
	// SSe : end point of a steady state
	// imin: from imin to detect the next steady state
	const int SSlengthMin=10; 
	// at least 10 samples can be viewed as steady state 
	char recording;
	while (1)
	{
		// detect steady state
		recording=0;
		for (i=imin; i<N; i++)
		{
			if (steady[i])
			{
				if (recording) // SS is recording
					++SSlength;
				else
				{	SSs=i; SSlength=1; recording=1; }
			}
			else if (recording) // SS is recording
			{
				if (SSlength>=SSlengthMin) // if length is enough
					break;
				else // length is too short to be a steady state
					recording=0;
			}
		}
		if (recording)
		{
			++NSS; // a SS is successfully found
			SSe=SSs+SSlength-1;
			//fprintf(stderr, "# %dth steady state %d:%d\n", NSS, SSs, SSe);

			/*
			// find mean and var
			sum=y[SSs]; sum2=y[SSs]*y[SSs];
			for (j=SSs+1; j<=SSe; j++)  
				// cannot use i, since i is used to detect end of the while loop
			{
				sum += y[j];
				sum2 += y[j]*y[j];
			}
			mean = sum/SSlength;
			sigma = sqrt(sum2/SSlength - mean*mean);
			//fprintf(stderr, "# mean=%f simga=%f\n", mean, sigma);
			//fprintf(stderr, "# max possible trend = %f\n", 
			//		2.0*sqrt(3.0)*sigma/dbW);
			*/
			SS[2*(NSS-1)]=SSs;
			SS[2*NSS-1]=SSe;
			imin=SSe+1;
		}
		
		if (i==N) 
		{//fprintf(stderr, "# Only %d steady state(s)\n", NSS); 
			break;}
	}
	*pNSS=NSS;
}

/* with statistics, i.e. mean & sigma*/
void SSD_statistics(int N, double *y, int window, int *pNSS, int *SS, 
		double *mean, double *sigma)
{
	if ( window <= 0 )
	{
		window = N/10;
		fprintf(stderr, "# set window width to %d\n", window);
	}
	else if ( window >= N )
	{
		fprintf(stderr, "# Error: Window width >= sample number\n");
		exit (1);
	}

	int i, j, lastSteady=0;
	double sum=y[0], sum2=y[0]*y[0], sum2ad=0.0;
	double dbW=(double)window;
	// moving window
	for (i=1; i<window; i++)
	{
		//mean += y[i];
		sum += y[i];
		sum2 += y[i]*y[i];
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1]);
	}
	char *steady=Calloc(char, N); POINTER_NULL(steady);
	if (sum2ad > (sum2-sum*sum/dbW))
	{
		for (i=0; i<window; i++)
			steady[i]=1;
		lastSteady = window-1;
	}

	// if I use moving overlaped window, then the resolution will be 1
	// for steady state intervals.
	// If partially overlaped, the resolution will be lower, and for
	// nonoverlapped window, the resolution is the window width, very low.
	for (i=window; i<N; i++) // i: end of the window
	{
		sum += (y[i]-y[i-window]);	
		sum2 += (y[i]*y[i] - y[i-window]*y[i-window]);
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1])
			-(y[i-window+1]-y[i-window])*(y[i-window+1]-y[i-window]);
		if (sum2ad > (sum2-sum*sum/dbW))
		{
			if ( lastSteady + window < i) // lastSteady is small
				j=i-window;
			else
				j=lastSteady;
			while (j<=i)
				steady[j++]=1;
			lastSteady = i;
		}
	}
	/*for (i=0; i<N; i++)
		printf("%d %d\n", i, (int)steady[i]);
	fprintf(stderr, "OK\n");*/

	/*===================================================================
	 * steady state dectection
	 *===================================================================*/
	int SSs=0, SSe, imin=0, NSS=0, SSlength=0;
	// SSs : starting point of a steady state
	// SSe : end point of a steady state
	// imin: from imin to detect the next steady state
	const int SSlengthMin=10; 
	// at least 10 samples can be viewed as steady state 
	char recording;
	while (1)
	{
		// detect steady state
		recording=0;
		for (i=imin; i<N; i++)
		{
			if (steady[i])
			{
				if (recording) // SS is recording
					++SSlength;
				else
				{	SSs=i; SSlength=1; recording=1; }
			}
			else if (recording) // SS is recording
			{
				if (SSlength>=SSlengthMin) // if length is enough
					break;
				else // length is too short to be a steady state
					recording=0;
			}
		}
		SSe=SSs+SSlength-1;
		//++NSS; // a SS is successfully found
		//fprintf(stderr, "# %dth steady state %d:%d\n", NSS, SSs, SSe);

		// find mean and var
		sum=y[SSs]; sum2=y[SSs]*y[SSs];
		for (j=SSs+1; j<=SSe; j++)  
			// cannot use i, since i is used to detect end of the while loop
		{
			sum += y[j];
			sum2 += y[j]*y[j];
		}
		mean[NSS] = sum/SSlength;
		sigma[NSS] = sqrt(sum2/SSlength - mean[NSS]*mean[NSS]);
		//fprintf(stderr, "# mean=%f simga=%f\n", mean, sigma);
		//fprintf(stderr, "# max possible trend = %f\n", 
		//		2.0*sqrt(3.0)*sigma/dbW);
		SS[2*NSS]=SSs;
		SS[2*(NSS++)+1]=SSe;
		imin=SSe+1;
		
		if (i==N) 
		{//fprintf(stderr, "# Only %d steady state(s)\n", NSS); 
			break;}
	}
	*pNSS=NSS;

	free(steady);
}


/* with statistics, i.e. mean & sigma*/
void SSD_statistics_batch(int N, double *y, int window, int *pNSS, int *SS, 
		double *mean, double *sigma, char *workspace)
{
	if ( window <= 0 )
	{
		window = N/10;
		fprintf(stderr, "# set window width to %d\n", window);
	}
	else if ( window >= N )
	{
		fprintf(stderr, "# Error: Window width >= sample number\n");
		exit (1);
	}

	int i, j, lastSteady=0;
	double sum=y[0], sum2=y[0]*y[0], sum2ad=0.0;
	double dbW=(double)window;
	// moving window
	for (i=1; i<window; i++)
	{
		//mean += y[i];
		sum += y[i];
		sum2 += y[i]*y[i];
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1]);
	}
	char *steady=workspace;
	for (i=0; i < N; i++) steady[i]=0;
	if (sum2ad > (sum2-sum*sum/dbW))
	{
		for (i=0; i<window; i++)
			steady[i]=1;
		lastSteady = window-1;
	}

	// if I use moving overlaped window, then the resolution will be 1
	// for steady state intervals.
	// If partially overlaped, the resolution will be lower, and for
	// nonoverlapped window, the resolution is the window width, very low.
	for (i=window; i<N; i++) // i: end of the window
	{
		sum += (y[i]-y[i-window]);	
		sum2 += (y[i]*y[i] - y[i-window]*y[i-window]);
		sum2ad += (y[i]-y[i-1])*(y[i]-y[i-1])
			-(y[i-window+1]-y[i-window])*(y[i-window+1]-y[i-window]);
		if (sum2ad > (sum2-sum*sum/dbW))
		{
			if ( lastSteady + window < i) // lastSteady is small
				j=i-window;
			else
				j=lastSteady;
			while (j<=i)
				steady[j++]=1;
			lastSteady = i;
		}
	}
	/*for (i=0; i<N; i++)
		printf("%d %d\n", i, (int)steady[i]);
	fprintf(stderr, "OK\n");*/

	/*===================================================================
	 * steady state dectection
	 *===================================================================*/
	int SSs=0, SSe, imin=0, NSS=0, SSlength=0;
	// SSs : starting point of a steady state
	// SSe : end point of a steady state
	// imin: from imin to detect the next steady state
	const int SSlengthMin=10; 
	// at least 10 samples can be viewed as steady state 
	char recording;
	while (1)
	{
		// detect steady state
		recording=0;
		for (i=imin; i<N; i++)
		{
			if (steady[i])
			{
				if (recording) // SS is recording
					++SSlength;
				else
				{	SSs=i; SSlength=1; recording=1; }
			}
			else if (recording) // SS is recording
			{
				if (SSlength>=SSlengthMin) // if length is enough
					break;
				else // length is too short to be a steady state
					recording=0;
			}
		}
		SSe=SSs+SSlength-1;
		//++NSS; // a SS is successfully found
		//fprintf(stderr, "# %dth steady state %d:%d\n", NSS, SSs, SSe);

		// find mean and var
		sum=y[SSs]; sum2=y[SSs]*y[SSs];
		for (j=SSs+1; j<=SSe; j++)  
			// cannot use i, since i is used to detect end of the while loop
		{
			sum += y[j];
			sum2 += y[j]*y[j];
		}
		mean[NSS] = sum/SSlength;
		sigma[NSS] = sqrt(sum2/SSlength - mean[NSS]*mean[NSS]);
		//fprintf(stderr, "# mean=%f simga=%f\n", mean, sigma);
		//fprintf(stderr, "# max possible trend = %f\n", 
		//		2.0*sqrt(3.0)*sigma/dbW);
		SS[2*NSS]=SSs;
		SS[2*(NSS++)+1]=SSe;
		imin=SSe+1;
		
		if (i==N) 
		{//fprintf(stderr, "# Only %d steady state(s)\n", NSS); 
			break;}
	}
	*pNSS=NSS;
}
