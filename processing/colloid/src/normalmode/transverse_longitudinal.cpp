#include "normalmode.h"

#include "colloid_base.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
	/*======================================================
	 * Eigenvalue manipulation
	 *====================================================*/
// assuming that the length scale is normalized by lattice constant, 
// so that the first Brillouin zone is simply as
// 		[-\pi, \pi]
void transverse_longitudinal(int& n, const char *file)
{
	// get mean position from .mp file
	int total_p;
	float *mp;
	const char mpsubfix[]=".mp";
	char *mpfilename=getfilename(file, mpsubfix);
	if ( readmp(total_p, &mp, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file \"%s\"\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);
	
	// get n'th normal mode from .ev file
	const char evsubfix[]=".ev";
	char *evfilename=getfilename(file, evsubfix);
	int Gn;
	double E;
	double *G;
	if ( readmode(n, Gn, E, &G, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", mpfilename);
		exit (1);
	}
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	E=1/sqrt(E);

	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");

	const double PI=acos(-1);

	const int nbin=61; // should set as odd
	double binsize=2*PI/nbin;
	double qlength, qx, qy, qr, Cos, Sin, qxe, qde, tR, tI, lR, lI;
	int i, j, k;
	float *mp_ptr=mp;
	double *G_ptr=G;
	double *tRe=(double *)malloc(total_p*sizeof(double));
	double *tIm=(double *)malloc(total_p*sizeof(double));
	double *lRe=(double *)malloc(total_p*sizeof(double));
	double *lIm=(double *)malloc(total_p*sizeof(double));
	if ( tRe==NULL )
		pERROR("double *tRe initialized wrong.");
	if ( tIm==NULL )
		pERROR("double *tIm initialized wrong.");
	if ( lRe==NULL )
		pERROR("double *lRe initialized wrong.");
	if ( lIm==NULL )
		pERROR("double *lIm initialized wrong.");

	double *T=(double *)malloc(nbin*nbin*sizeof(double));
	double *L=(double *)malloc(nbin*nbin*sizeof(double));
	// since the componets of -q is same as the components of q,
	// so only half is calculated.
	for (i=0; i<nbin; i++)
	{
		qx=(i+0.5)*binsize-PI;
		for (j=nbin/2; j<nbin; j++)
		{
			qy=(j+0.5)*binsize-PI;
			
			mp_ptr=mp;
			G_ptr=G;

			for (k=0; k<total_p; k++)
			{
				qr=qx*(*mp_ptr)+qy*(*(mp_ptr+1));
				mp_ptr+=2;
				
				Cos=cos(qr);
				Sin=sin(qr);

				qxe=qx*(*(G_ptr+1))-qy*(*G_ptr);
				qde=qx*(*G_ptr)+qy*(*(G_ptr+1));
				G_ptr+=2;

				tRe[k]=qxe*Cos;
				tIm[k]=qxe*Sin;

				lRe[k]=qde*Cos;
				lIm[k]=qde*Sin;
			}
			
			tR=pairwise(total_p, tRe);
			tI=pairwise(total_p, tIm);

			lR=pairwise(total_p, lRe);
			lI=pairwise(total_p, lIm);

			qlength=sqrt(qx*qx+qy*qy);

			T[i*nbin+j]=(tR*tR+tI*tI)/qlength;
			T[(nbin-1-i)*nbin+nbin-1-j]=T[i*nbin+j];
			
			L[i*nbin+j]=(lR*lR+lI*lI)/qlength;
			L[(nbin-1-i)*nbin+nbin-1-j]=L[i*nbin+j];
		}
	}
	free(mp);
	free(tRe);
	free(tIm);
	free(lRe);
	free(lIm);
	printf("# T/L component calculation done!\n");

	//===================================================
	// print out the data
	char TLsubfix [10];
	sprintf(TLsubfix, "_%d.T", n);
	char *Tfilename = getfilename(file, TLsubfix);
	FILE *Tfile=fopen(Tfilename, "w");
	if (Tfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", Tfilename);
		exit (1);
	}

	i=0;
	while (TLsubfix[i]!='\0')
		++i;
	TLsubfix[i-1]='L';
	char *Lfilename = getfilename(file, TLsubfix);
	FILE *Lfile=fopen(Lfilename, "w");
	if (Lfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", Lfilename);
		exit (1);
	}


	fprintf(Tfile, "# transverse component of %d'th normal mode\n"
			"# eigenvalue = %f\n", n, E);

	fprintf(Lfile, "# longitudinal component of %d'th normal mode\n"
			"# eigenvalue = %f\n", n, E);

	double *Tpointer=T, *Lpointer=L;
	for (i=0; i<nbin; i++)
	{
		qx=(i+0.5)*binsize-PI;
		for (j=0; j<nbin; j++)
		{
			qy=(j+0.5)*binsize-PI;
			fprintf(Tfile, "%f %f %f\n", qx, qy, *(Tpointer++));
			fprintf(Lfile, "%f %f %f\n", qx, qy, *(Lpointer++));
			//fprintf(Tfile, "%f ", *(Tpointer++));
			//fprintf(Lfile, "%f ", *(Lpointer++));
		}
		fprintf(Tfile, "\n");
		fprintf(Lfile, "\n");
	}
	fclose(Tfile);
	fclose(Lfile);

	free(T);
	free(L);

	//=========================================================
	// plot by GNUplot, using pipe
	
	/*
	FILE *Tpipe=popen("gnuplot -persist", "w");
	if (Tpipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"Tpipe\"!\n");
		exit (1);
	}

	fprintf(Tpipe, "reset\n"
				   "unset key\n"
				   "set pm3d map\n"
				   "set palette rgbformulae 33,13,10\n"
				   "set size ratio -1\n"
				   "set title \"transverse component of %d'th normal mode\\nomega=%f\"\n"
				   "set xlabel 'qx'\n"
				   "set ylabel 'qy'\n", n, E);

	fprintf(Tpipe, "set xrange [-%f:%f]\n"
				   "set yrange [-%f:%f]\n",
				   PI, PI, PI, PI);

	fprintf(Tpipe, "splot '%s' using 1:2:3 with image\n", Tfilename);

	fclose(Tpipe);


	FILE *Lpipe=popen("gnuplot -persist", "w");
	if (Lpipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"Lpipe\"!\n");
		exit (1);
	}

	fprintf(Lpipe, "reset\n"
				   "unset key\n"
				   "set pm3d map\n"
				   "set palette rgbformulae 33,13,10\n"
				   "set size ratio -1\n"
				   "set title \"longitudinal component of %d'th normal mode\\nomega=%f\"\n"
				   "set xlabel 'qx'\n"
				   "set ylabel 'qy'\n", n, E);

	fprintf(Lpipe, "set xrange [-%f:%f]\n"
				   "set yrange [-%f:%f]\n",
				   PI, PI, PI, PI);

	fprintf(Lpipe, "splot '%s' using 1:2:3 with image\n", Lfilename);

	fclose(Lpipe);
	*/


	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}
	char NULLsubfix [2];
	NULLsubfix[0]='\0';
	char *NULLfilename = getfilename(file, NULLsubfix);

	fprintf(pipe,  "reset\n"
				   "unset key\n"
				   "set pm3d map\n"
				   "set palette rgbformulae 33,13,10\n"
				   "set size ratio -1\n"
				   "set multiplot layout 1, 2 title "
				   "\"%s\\ntransverse and longitudinal components of"
				   " %d'th normal mode\\nomega=%f\"\n"
				   "set xlabel 'qx'\n"
				   "set ylabel 'qy'\n", NULLfilename, n, E);

	fprintf(pipe, "set xrange [-%f:%f]\n"
				  "set yrange [-%f:%f]\n",
				  PI, PI, PI, PI);

	fprintf(pipe, "set title \"transverse\"\n"
				  "splot '%s' using 1:2:3 with image\n", Tfilename);
	fprintf(pipe, "set title \"longitudinal\"\n"
				  "splot '%s' using 1:2:3 with image\n", Lfilename);

	fclose(pipe);

	free(Tfilename);
	free(Lfilename);
	free(NULLfilename);
}


void transverse_longitudinal_eps(int& n, const char *file)
{
	char TLsubfix [10];
	sprintf(TLsubfix, "_%d.T", n);
	char *Tfilename = getfilename(file, TLsubfix);
	FILE *pfile=fopen(Tfilename, "r");
	if (pfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", Tfilename);
		exit (1);
	}
	fclose(pfile);

	int i=0;
	while (TLsubfix[i]!='\0')
		++i;
	TLsubfix[i-1]='L';
	char *Lfilename = getfilename(file, TLsubfix);
	pfile=fopen(Lfilename, "r");
	if (pfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", Lfilename);
		exit (1);
	}
	fclose(pfile);

	// get Gn and E
	const char evsubfix [] = ".ev";
	char *evfilename = getfilename(file, evsubfix);
	pfile=fopen(evfilename, "r");
	if (pfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", evfilename);
		exit (1);
	}
	free(evfilename);
	
	int *header=(int *)malloc(2*sizeof(int));
	if ( fread(header, sizeof(header[0]), 2, pfile) != 2 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", evfilename);
		exit (1);
	}

	if ( header[0] != ('e'*256+'v') )
	{
		fprintf(stderr, "# Error: %s is not a .ev file!\n", evfilename);
		exit (1);
	}

	//int Gn=header[1];
	free(header);


	long int offset=(long int)(n*sizeof(double));
	if (fseek( pfile, offset, SEEK_CUR ) != 0 )
		pERROR("reading eigenvalue failed.");
	double * Ep=(double *)malloc(sizeof(double));

	if ( (int)fread(Ep, sizeof(double), 1, pfile) != 1 )
	{
		fprintf(stderr, "# Error: read from file %s failed!\n", evfilename);
		exit (1);
	}

	fclose(pfile);

	char epsTLsubfix [15];
	sprintf(epsTLsubfix, "_TL_%d.eps", n);
	char *eps = getfilename(file, epsTLsubfix);

	const double PI=acos(-1);

	FILE *pipe=popen("gnuplot", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}
	
	char NULLsubfix [2];
	NULLsubfix[0]='\0';
	char *NULLfilename = getfilename(file, NULLsubfix);

	fprintf(pipe,  "reset\n"
				   "set terminal postscript eps color enhanced\n"
				   "set output \"%s\"\n"
				   "unset key\n"
				   "set pm3d map\n"
				   "set palette rgbformulae 33,13,10\n"
				   "set size ratio -1\n"
				   "set multiplot layout 1, 2 title "
				   "\"%s\\ntransverse and longitudinal components of"
				   " %d'th normal mode\\nomega=%f\"\n"
				   "set xlabel 'qx'\n"
				   "set ylabel 'qy'\n", 
				   eps, NULLfilename, n, 1/sqrt(*Ep));

	free(Ep);

	fprintf(pipe, "set xrange [-%f:%f]\n"
				  "set yrange [-%f:%f]\n",
				  PI, PI, PI, PI);

	fprintf(pipe, "set title \"transverse\"\n"
				  "splot '%s' using 1:2:3 with image\n", Tfilename);
	fprintf(pipe, "set title \"longitudinal\"\n"
				  "splot '%s' using 1:2:3 with image\n", Lfilename);

	fclose(pipe);

	free(Tfilename);
	free(Lfilename);
	free(eps);
	free(NULLfilename);
}
