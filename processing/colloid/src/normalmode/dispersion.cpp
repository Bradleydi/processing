#include "normalmode.h"

#include "colloid_base.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

using namespace std;
	
// assuming that the length scale is normalized by lattice constant, 
// so that the first Brillouin zone is simply as
// 		[-\pi, \pi]
void dispersion(const char *file)
{
	/*=================================================================
	 * get mean position from .mp file
	 *=================================================================*/
	int total_p;
	float *mp;
	const char mpsubfix[]=".mp";
	char *mpfilename=getfilename(file, mpsubfix);
	if ( readmp(total_p, &mp, mpfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file \"%s\"\n",
				mpfilename);
		exit (1);
	}
	printf("# File \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);


	/*=================================================================
	 * get normal modes from .ev file
	 *=================================================================*/
	const char evsubfix[]=".ev";
	char *evfilename=getfilename(file, evsubfix);
	int Gn;
	double *E;
	double *G;
	if ( readev(Gn, &E, &G, evfilename) != 0 )
	{
		fprintf(stderr, "# Error: something wrong during reading file %s\n", 
				evfilename);
		exit (1);
	}
	printf("# File \"%s\" loaded!\n", evfilename); 
	free(evfilename);

	double *omega=(double *) malloc ((Gn-2)*sizeof(double));
	POINTER_NULL (omega);

	int i;
	for (i=0; i<Gn-2; i++)
		omega[i]=1/sqrt(E[i+2]);

	free (E);

	if ( Gn != 2*total_p )
		pERROR("reading files wrong.");

	const double PI=acos(-1);

	const int nbin=61; // should set as odd
	double binsize=2*PI/nbin;
	double qlength, qx, qy, qr, Cos, Sin, qxe, qde, tR, tI, lR, lI;
	int j, k;
	float *mp_ptr=mp;
	double *G_ptr=G;
	double *tRe=(double *)malloc(total_p*sizeof(double));
	double *tIm=(double *)malloc(total_p*sizeof(double));
	double *lRe=(double *)malloc(total_p*sizeof(double));
	double *lIm=(double *)malloc(total_p*sizeof(double));
	POINTER_NULL(tRe);
	POINTER_NULL(tIm);
	POINTER_NULL(lRe);
	POINTER_NULL(lIm);

	double *T=(double *)malloc(nbin*nbin*sizeof(double));
	POINTER_NULL(T);
	double *L=(double *)malloc(nbin*nbin*sizeof(double));
	POINTER_NULL(L);
	
	
	//===================================================
	// angle average
	const int rmax=nbin/2;
	const int nr=rmax+1;
	double *Tavg=(double *)calloc(nr, sizeof(double));
	POINTER_NULL(Tavg);
	double *Lavg=(double *)calloc(nr, sizeof(double));
	POINTER_NULL(Lavg);

	double *count=(double *)calloc(nr, sizeof(double));
	POINTER_NULL(count);

	double r, dR, data; 
	int R, x2;
	for (i=-rmax; i<=rmax; i++)
	{
		x2=i*i;
		for (j=-rmax; j<=rmax; j++)
		{
			r=sqrt(x2+j*j);
			R=(int)r;
			dR=r-R;

			if (R<rmax)
			{
				count[R]+=1-dR;
				count[R+1]+=dR;
			}
			else if (R==rmax)
				count[R]+=1-dR;
		}
	}
	// since the componets of -q is same as the components of q,
	// so only half is calculated.
	int omg;
	
	
	const char subfix []= ".dp";
	char *dpfilename = getfilename(file, subfix);
	FILE *dpfile=fopen(dpfilename, "w");
	//FILE *dpfile=fopen(dpfilename, "wb");
	FILE_NULL(dpfile, dpfilename);

	fprintf(dpfile, "# q omage T L\n");
	/*
	int *header = (int *)malloc (3*sizeof(int));
	header[0]='d'*256+'p';
	header[1]=Gn;
	header[2]=nr;

	if ( fwrite(header, sizeof(header[0]), 3, dpfile) != 3 )
	{
		fprintf(stderr, "# Error: write to file \"%s\" failed!\n", dpfilename);
		exit (1);
	}
	free(header);
	
	double *tmp = (double *) malloc (4*sizeof(double));
	POINTER_NULL(tmp);
	*/
	
	double *Gp=G+Gn*(Gn-1);
	for (omg=Gn-3; omg>=0; omg--)
	{
		for (i=0; i<nbin; i++)
		{
			qx=(i+0.5)*binsize-PI;
			for (j=rmax; j<nbin; j++) // rmax=nbin/2
			{
				qy=(j+0.5)*binsize-PI;
				
				mp_ptr=mp;
				G_ptr=Gp;

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

		Gp -= Gn;
		
		
		// angle average
		/*
		for (i=0; i<nr; i++)
		{
			Tavg[i]=0;
			Lavg[i]=0;
		}
		*/
		for (i=-rmax; i<=rmax; i++) // rmax=nbin/2
		{
			x2=i*i;
			for (j=-rmax; j<=rmax; j++)
			{
				r=sqrt(x2+j*j);
				R=(int)r;

				if ( R < rmax )
				{
					dR=r-R;

					k=(i+rmax)*nbin+j+rmax;
					
					data=T[k];
					Tavg[R]+=data*(1-dR);
					Tavg[R+1]+=data*dR;
					
					data=L[k];
					Lavg[R]+=data*(1-dR);
					Lavg[R+1]+=data*dR;
				}
				else if ( R == rmax )
				{
					dR=r-R;

					k=(i+rmax)*nbin+j+rmax;
					
					data=T[k];
					Tavg[R]+=data*(1-dR);
					
					data=L[k];
					Lavg[R]+=data*(1-dR);
				}
			}
		}

		for (i=0; i<nr; i++)
		{
			if (count[i] > 1.0e-6 )
			{
				fprintf(dpfile, "%f %f %f %f\n", 
						i*binsize, omega[omg],
						Tavg[i]/count[i],
						Lavg[i]/count[i]);
				Tavg[i]=0;
				Lavg[i]=0;
				

				/*
				tmp[0]=i*binsize;
				tmp[1]=omega[omg];
				tmp[2]=Tavg[i]/count[i];
				tmp[3]=Lavg[i]/count[i];


				if ( fwrite(tmp, sizeof(tmp[0]), 4, dpfile) != 4 )
				{
					fprintf(stderr, "# Error: write to file \"%s\" failed!\n",
							dpfilename);
					exit (1);
				}
				Tavg[i]=0;
				Lavg[i]=0;
				*/

			}
			else
			{
				printf("# wrong: count[%d]=%f\n", i, count[i]);
				Tavg[i]=0;
				Lavg[i]=0;
			}
		}
		if ( omg%10 == 0 )
			printf("# omg=%d\n", omg);
	}

	fclose(dpfile);
	free(dpfilename);

	free(tRe);
	free(tIm);
	free(lRe);
	free(lIm);
	free(T);
	free(L);
	free(count);
	free(Tavg);
	free(Lavg);
	//free(tmp);
	free(G);
	free(omega);
	free(mp);

	printf("# T/L component calculation done!\n");
}


void plot_dispersion(const char *file, int N)
{
	int Gn;
	double *E;
	readE(Gn, &E, file);
	
	double min=1/sqrt(E[Gn-1]);
	double max=1/sqrt(E[2]);

	if ( N <= 1 || N > Gn-2 )
		printf("# N is set to %d to include all freqencies.\n", Gn-2);
	else 
		max=1/sqrt(E[Gn-N]);
	free(E);

	double PI=acos(-1);

	char *dpfile = getfilename(file, ".dp");
	char *epsfile = getfilename(file, "_dp.eps");

	// read max of dp
	FILE *dpf = fopen (dpfile, "r");
	FILE_NULL(dpf, dpfile);

	const int buffersize=80;
	char *buffer = (char *)malloc(buffersize*sizeof(char));
	POINTER_NULL(buffer);

	fgets(buffer, buffersize, dpf);

	float *data = (float *)malloc(4*sizeof(float));
	POINTER_NULL (data);

	float TLmax=0;
	while ( fgets(buffer, buffersize, dpf) != NULL )
	{
		if ( sscanf(buffer, "%f %f %f %f", data, data+1, data+2, data+3) != 4)
		{
			fprintf(stderr, "# Error when reading \"%s\"\n", dpfile);
			exit (1);
		}

		if ( data[2] > TLmax )
			TLmax = data[2];
		if ( data[3] > TLmax )
			TLmax = data[3];
	}
	fclose(dpf);
	free(data);


	//=========================================================
	// plot by GNUplot, using pipe
	FILE *pipe=popen("gnuplot -persist", "w");
	if (pipe==NULL)
	{
		fprintf(stderr, "# Error: can not open pipe \"pipe\"!\n");
		exit (1);
	}

	fprintf(pipe, "reset\n"
				  "set terminal postscript eps color enhanced font 20\n"
				  "set output \"%s\"\n"
				  "unset key\n"
                  "set pm3d map\n"
				  "#set palette rgbformulae 33,13,10\n"
				  "#set palette grey\n"
				  "#set colorbox horizontal user origin 0.1,0.5 size .8,0.02\n"
				  "set xrange [0:%f]\n"
				  "set yrange [%f:%f]\n"
				  "set cbrange [0:%f] \n"
				  "set xtics(\"0\" 0, \"{/Symbol p}/4\" %f, "
				  "\"{/Symbol p}/2\" %f, \"3{/Symbol p}/4\" %f, "
				  "\"{/Symbol p}\" %f)\n"
				  "set xlabel \"a_0q\"\n"
				  "set ylabel \"{/Symbol W}\"\n"
				  "set multiplot layout 1,2\n"
				  "#set logscale y\n"
				  "#set logscale cb\n"
				  "unset colorbox\n"
				  "set origin 0.15, 0.0\n"
				  "splot \"%s\" using 1:2:3 with image palette\n"
				  "set colorbox\n"
				  "set origin 0.5, 0.0\n"
				  "#set size 0.5, 0.6\n"
				  "unset ylabel\n"
				  "#unset ytics\n"
				  "set format y \"\"\n"
				  "splot \"%s\" using 1:2:4 with image palette\n"
				  "#set origin 0.4, 0.0\n"
				  "set colorbox\n"
				  "unset multiplot\n"
				  "exit\n",
				  epsfile, PI, min, max, TLmax, 
				  PI/4, PI/2, 3*PI/4,PI,
				  dpfile, dpfile);

	fclose(pipe);
	free(dpfile);
	free(epsfile);
}
/*
void plot_dispersion(const char *file)
{
	const char *subfix = ".dp";
	char *dpfilename = getfilename(file, subfix);
	FILE *dpfile=fopen(dpfilename, "w");
	if (dpfile==NULL)
	{
		fprintf(stderr, "# Error: can not open file %s!\n", dpfilename);
		exit (1);
	}
	free(dpfilename);
	//===================================================
	// print out the data
	const char subfix [10];
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
*/	
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

/*
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
*/
