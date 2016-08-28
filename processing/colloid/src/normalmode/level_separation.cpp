#include "normalmode.h"

#include "colloid_base.h"
#include "lapack.h"
#include "miscellaneous.h"
#include "statistics.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;


void level_separation_w(int Gn, double *E, const char *file)
{
	const char subfix[]=".sp";
	char *filename=getfilename(file, subfix);

	FILE *spf=fopen(filename, "w");
	if (spf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename);
		exit (1);
	}
	
	// get omega
	double *omega=(double *)malloc((Gn-2)*sizeof(double));
	if (omega==NULL)
		pERROR("double *omega initialized error!");
	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);
	
	// speration of eigen frequencies, level repulsive or not.
	int nsp=Gn-3;
	double *sp=(double *)malloc(nsp*sizeof(double));
	if (sp==NULL)
		pERROR("double *sp initialized error!");

	for (i=0; i<nsp; i++)
		sp[i]=omega[i]-omega[i+1];
	
	/* mean separation */
	double msp=pairwise(nsp, sp)/nsp;
	double *sp2=(double *)malloc(nsp*sizeof(double));
	if (sp2==NULL)
		pERROR("double *sp2 initialized error!");
	for (i=0; i<nsp; i++)
	{
		sp2[i]=sp[i]-msp;
		sp2[i]*=sp2[i];
	}

	fprintf(spf, "# Level separation:\n");
	fprintf(spf, "# Mean speration   :   %6.6f\n", msp);
	fprintf(spf, "# Standard der.    :   %6.6f\n", 
			pairwise(nsp, sp2)/nsp);

	free(sp2);
	
	
	fprintf(spf, "\n# Level separation:\n");
	fprintf(spf, "# omega\tseparation\n");
	for (i=0; i<nsp; i++)
		fprintf(spf, "%6.6f\t%6.6f\n", (omega[i]+omega[i+1])/2, sp[i]);
	fclose(spf);
	free(omega);
	

	// distribution of level separation (statistics)
	const char subfix1[]=".dsp";
	char *filename1=getfilename(file, subfix1);

	FILE *dspf=fopen(filename1, "w");
	if (dspf==NULL)
	{
		fprintf(stderr, "# Error: file %s can not open!\n", filename1);
		exit (1);
	}

	const int nbin=100;

	double max=sp[0], min=sp[0];
	for (i=1; i<nsp; i++)
	{
		if( min > sp[i] )
			min = sp[i];
		else if ( max < sp[i] )
			max = sp[i];
	}
	double binsize=(1+1.e-6)*(max-min)/nbin;
	unsigned int *hist=(unsigned int *)calloc(nbin, sizeof(unsigned int));
	if (hist==NULL)
		pERROR("unsigned int *hist initialized error!");
	for (i=0; i<nsp; i++)
		++hist[(int)((sp[i]-min)/binsize)];

	fprintf(dspf, "\n# Distribution of level separation:\n");
	fprintf(dspf, "# separation\tPrabability\n");
	for (i=0; i<nbin; i++)
	{
		if (hist[i]!=0)
		{
			fprintf(dspf, "%6.6f\t%6.6f\n",
				    (i+0.5)*binsize+min,
				    (double)(hist[i])/binsize/nsp);
		}
	}
	fclose(dspf);

	free(sp);
	free(hist);
	free(filename);
	free(filename1);
}


void level_separation(colloid_base &ptid, const char *file, bool Remove_Drift)
{
	int Gn;
	double *E, *G;
	if ( readev(Gn, &E, &G, file) != 0 )
	{	
		printf("# Calling normal_mode!\n");
		normal_mode(ptid, file, Remove_Drift);
	}
	else
	{
		free(G);

		level_separation_w(Gn, E, file);
		free(E);
	}
}


// logarithmic binning
void level_separation_log(const char *file)
{
	int Gn;
	double *E;

	readE(Gn, &E, file);

	// get omega
	double *omega=(double *)malloc((Gn-2)*sizeof(double));
	POINTER_NULL(omega);

	int i, j;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);
	free(E);

	// speration of eigen frequencies, level repulsive or not.
	int nsp=Gn-3;
	double *sp=(double *)malloc(nsp*sizeof(double));
	POINTER_NULL(sp);

	for (i=0; i<nsp; i++)
		sp[i]=omega[i]-omega[i+1];
	free(omega);


	// distribution of level separation (statistics)
	char *filename=getfilename(file, ".dlsl");

	FILE *dspf=fopen(filename, "w");
	FILE_NULL(dspf, filename);	

	free(filename);
	
	
	/* mean separation */
	double msp=pairwise(nsp, sp)/nsp;
	double *sp2=(double *)malloc(nsp*sizeof(double));
	POINTER_NULL(sp2);

	for (i=0; i<nsp; i++)
	{
		sp2[i]=sp[i]-msp;
		sp2[i]*=sp2[i];
	}

	fprintf(dspf, "# Level separation:\n");
	fprintf(dspf, "# Mean speration   :   %6.6f\n", msp);
	fprintf(dspf, "# Standard der.    :   %6.6f\n", 
			pairwise(nsp, sp2)/nsp);

	free(sp2);

	fprintf(dspf, "\n# Distribution of level separation:\n");
	fprintf(dspf, "# Delta_omega probability histogram\n");


	sort(nsp, sp);
	double max_min=sp[nsp-1]/sp[0];

	for (i=1; i<nsp; i++)
		if (sp[i-1]>sp[i])
			printf("wrong\n");


	// get bins
	int bin [150];
	bin[0]=1;
	int nbin=1;
	for (i=0; i<150; i++)
	{
		j=(int)(pow(1.15,i)+0.5);
		if ( j>max_min )
			break;
		if (j==bin[nbin-1])
			continue;
		bin[nbin++]=j;
	}

	// since bin[i] may be very large for large i, so that
	// bin[i-1]*bin[i] maybe overflow
	// use sqrtbin to avoid that
	double *sqrtbin=(double *)malloc(nbin*sizeof(double));
	POINTER_NULL(sqrtbin);
	for (i=0; i<nbin; i++)
		sqrtbin[i]=sqrt((double)bin[i]);
	
	double s;
	unsigned int hist=1;
	j=1;
	for (i=1; i<nsp; i++)
	{
		s=sp[i]/sp[0];
		if ( s > bin[nbin-1] )
		{
			fprintf(dspf, "%6.6f %6.6f %d\n",
					sqrtbin[j-1]*sqrtbin[j]*sp[0],
					(double)(hist)/nsp/(bin[j]-bin[j-1])/sp[0], hist);
			break;
		}

		if ( s < bin[j] )
			++hist;
		else
		{
			fprintf(dspf, "%6.6f %6.6f %d\n",
					sqrtbin[j-1]*sqrtbin[j]*sp[0],
					(double)(hist)/nsp/(bin[j]-bin[j-1])/sp[0], hist);
			hist=1;
			while ( s >= bin[j] )
				++j;
		}
	}

	free(sqrtbin);
	fclose(dspf);
	free(sp);
}


void level_spacing(const char *file, double h, unfolding_t UNFOLDING)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);

	// get omega
	double *omega=(double *)malloc((Gn-2)*sizeof(double));
	POINTER_NULL(omega);

	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[Gn-1-i]);
	free(E);

	// unique omega
	int Nnew=1;
	for (i=1; i<Gn-2; i++)
		if ( omega[i] > omega[i-1] + 1.0e-7 )
			omega[Nnew++] = omega[i];

	double *X=unfolding(Nnew, omega, h, UNFOLDING);
	free(omega);

	const int Ns=Nnew-1;
	double *s=Malloc(double, Ns); POINTER_NULL(s);
	for (i=0; i<Ns; i++)
		s[i]=X[i+1]-X[i];
	free(X);

	double min=s[0], max=s[0];
	for (i=1; i<Ns; i++)
	{
		if ( min > s[i] )
			min = s[i];
		else if ( max < s[i] )
			max = s[i];
	}

	const int nbin=5*(int)sqrt(Ns);
	const double binsize=(max-min)*(1+1.0e-6)/nbin;
	unsigned int *hist=Calloc(unsigned int, nbin);
	POINTER_NULL(hist);
	for (i=0; i<Ns; i++)
		++hist[(int)((s[i]-min)/binsize)];

	//=============================================================
	// output results
	char *filename=getfilename(file, ".dls");
	FILE *dspf=fopen(filename, "w");
	FILE_NULL(dspf, filename);	
	free(filename);
	
	/* mean spacing */
	double ms=pairwise(Ns, s)/Ns;
	double *s2=(double *)malloc(Ns*sizeof(double));
	POINTER_NULL(s2);
	for (i=0; i<Ns; i++)
	{
		s2[i]=s[i]-ms;
		s2[i]*=s2[i];
	}

	fprintf(dspf, "# Level separation:\n");
	fprintf(dspf, "# Ratio of nondegenerate levels:   %.6f ( %d/%d )\n",
			      (double)Nnew/(Gn-2), Nnew, Gn-2);
	fprintf(dspf, "# Mean speration   :   %6.6f\n", ms);
	fprintf(dspf, "# Standard der.    :   %6.6f\n", 
			pairwise(Ns, s2)/Ns);
	free(s2);
	free(s);
	if ( UNFOLDING == KDE )
		fprintf(dspf, "# unfolding by kernel density estimation\n");
	else if ( UNFOLDING == LLR )
		fprintf(dspf, "# unfolding by local linear regression\n");

	fprintf(dspf, "\n# Distribution of level separation:\n");
	fprintf(dspf, "# Delta_omega_bar(unfolded) probability histogram\n");
	for (i=0; i<nbin; i++)
		fprintf(dspf, "%6.6f %6.6f %d\n",
					  min+(i+0.5)*binsize,
					  ((double)hist[i])/Ns/binsize, hist[i]);
	fclose(dspf);
	free(hist);
}


void level_rigidity(const char *file, double h, unfolding_t UNFOLDING)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);

	// get omega
	double *omega=(double *)malloc((Gn-2)*sizeof(double));
	POINTER_NULL(omega);

	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[Gn-1-i]);
	free(E);

	// unique omega
	int Nnew=1;
	for (i=1; i<Gn-2; i++)
		if ( omega[i] > omega[i-1] + 1.0e-7 )
			omega[Nnew++] = omega[i];

	double *X=unfolding(Nnew, omega, h, UNFOLDING);
	free(omega);
	
	char *Dfile = getfilename(file, ".D3");
	FILE *D = fopen(Dfile, "w");
	FILE_NULL(D, Dfile);
	free(Dfile);
	fprintf(D, "# L   Delta3\n");

	double *d3 = (double *)malloc(Nnew*sizeof(double));
	POINTER_NULL(d3);
	double *xx = (double *)malloc(Nnew*sizeof(double));
	POINTER_NULL(xx);
	double *xi = (double *)malloc(Nnew*sizeof(double));
	POINTER_NULL(xi);
	
	double sumx,sumx2, sum, D3;
	int L, n, x;
	//for (L=3; L<Nnew; L++)
	for (L=3; L<200; L++)
	{
		for (x=0; (X[x]+L) <= X[Nnew-1]; x++)
		{
			//========================================================
			// find all X in [X[x], X[x]+L]
			n=0;
			while ( X[x+n] <= (X[x]+L) )
			{
				xx[n]=X[x+n]-X[x]-L/2.0;
				++n;
			}

			sumx=pairwise(n, xx);
			for (i=0; i<n; i++)
			{
				xi[i]=(n-2*i-1)*xx[i];
				xx[i]=xx[i]*xx[i];
			}
			sumx2=pairwise(n, xx);
			sum=pairwise(n, xi);
			
			d3[x]=n*n/16.0-sumx*sumx/L/L+1.5*n*sumx2/L/L
				  -3*sumx2*sumx2/L/L/L/L
				  +sum/L;
		}
		D3=pairwise(x, d3)/x;
		fprintf(D, "%d %6.6f\n", L, D3);
	}
	fclose(D);
	free(xx);
	free(xi);
	free(d3);
	free(X);
}

	/*
	double max=sp[0], min=sp[0];
	for (i=1; i<nsp; i++)
	{
		if( min > sp[i] )
			min = sp[i];
		else if ( max < sp[i] )
			max = sp[i];
	}
	double max_min=max/min;

	unsigned int *hist=(unsigned int *)calloc((nbin-1), sizeof(unsigned int));
	POINTER_NULL(hist);
	for (i=0; i<nsp; i++)
	{
		++hist[(int)((sp[i]/min))];
	}
	*/

	/*
	fprintf(dspf, "\n# Distribution of level separation:\n");
	fprintf(dspf, "# separation\tPrabability\n");
	for (i=0; i<nbin; i++)
	{
		if (hist[i]!=0)
		{
			fprintf(dspf, "%6.6f\t%6.6f\n",
				    (i+0.5)*binsize+min,
				    (double)(hist[i])/binsize/nsp);
		}
	}
	fclose(dspf);

	free(sp);
	free(hist);
	free(filename);
	free(filename1);
	*/
