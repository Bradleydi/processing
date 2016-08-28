#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

/* Thu Nov  1 12:26:47 HKT 2012
 * soft spot
 * Actually, here I just calculate the sum of participation of a few
 * low frequency modes, and the regions with large low participation
 * are idenfitied by eyes as the soft spot.
 * Maybe a cluster program should be included, but this needs a program
 * to calculate the nearest neighbors of ellipses.
 * The definition NN of ellipses is complicated, and will not be implemented
 * in this version*/

/* Implementation details
 * Due to low freq. modes are important than high freq. modes, and they
 * also dominate the motions of the system. Hence I choose number of low
 * freq. modes that the contribution of them exeeds a given threshold, say, 
 * 90%, to use as the identification of softspots.
 * I don't want to implemented a program for read a given number of modes,
 * maybe this is much more effient.
 * Here in this version, I can
 * 1. try to make use of 'e_readmode', to read the mode one by one.
 * 2. open the .ev file here, and read the data
 * And I implemented method 2.
 * To make the reading easier, I just read all eigenvalues together first.
 *
 * The contribution can be calculated as
 * c_w = <(e_w \cdot u)^2/(u \cdot u)>,
 * i.e. the mean projection of displacement field on the mode
 * Even <e^Tu u^Te>= <e^T(uu^T)e>=e^T<uu^T>e=e^TCe=E,
 * we cannot directly use it, since 
 * <e^T (uu^T) e / u^Tu> !=<e^T (uu^T) e> / <u^Tu>
 * However, since we don't want to calculate the displacement u at all,
 * hence we want to avoid this method.
 * If we assume that <u^Tu>, which is N times of the vaiance of u, 
 * is a constant over time, then above inequility can be approximately viewed
 * as an equility, and then
 * c_w=E/<u^Tu>.
 * Since <u^Tu>=Tr<uu^T>=\sum E, hence c_w = E/(sum E).
 * So we can use the eigenvalues itself to get the contribution
 *
 * We can also make the output EPS figures contained in one file, but this
 * may be not good for showing, and also not good for including other files.
 */

void e_softspot(const char *file, double a, double b, double a0, 
		double threshold)
{	
	char *filename=getfilename(file, ".ev");
	FILE *evf=fopen(filename, "rb");
	FILE_NULL(evf, filename);

	int * header=(int *)malloc(2*sizeof(int)); POINTER_NULL(header);
	FreadN(header, 2, evf, filename);
	if ( header[0] != ('e'*256+'v'+'e') ) {
		fprintf(stderr, "# Error: '%s' is not a .ev file!\n", filename);
		exit (1);
	}
	int Gn=header[1];
	int total_p=Gn/3;
	free(header);
	
	double *E=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(E);
	double *G=(double *)malloc(Gn*sizeof(double)); POINTER_NULL(G);
	double *pT=(double *)calloc(total_p, sizeof(double)); POINTER_NULL(pT);
	double *pR=(double *)calloc(total_p, sizeof(double)); POINTER_NULL(pR);
	double *xG=G, *yG=xG+total_p, *aG=yG+total_p;
	// read eigenvalue
	FreadN(E, Gn, evf, filename);
	int i,j, N;
	double sum=E[0];
	for (i=1; i<Gn; i++) sum += E[i];
	// find the contributed modes
	i=Gn-1;
	double contribution=E[i--];
	threshold *= sum;
	while (contribution<=threshold)
		contribution += E[i--];
	N=i+1;
	long offset=(long int)(2*sizeof(int))+(long int)(Gn*(N+1)*sizeof(double));
	if ( fseek( evf, offset, SEEK_SET ) != 0 )
		pERROR("reading eigenvector failed.");

	for (i=N; i<Gn; i++)
	{
		// read eigen vector
		FreadN(G, Gn, evf, filename);
		for (j=0; j<total_p; j++)
		{
			pT[j] += xG[j]*xG[j] + yG[j]*yG[j];
			pR[j] += aG[j]*aG[j];
		}
	}
	Fclose(evf, filename);
	free(filename);
	free(E); free(G);

	//==================================================================
	// print the result to eps files
	//==================================================================

	//get mean position
	float *mp;
	char *mpfilename=getfilename(file, ".mp");
	e_readmp(&total_p, &mp, mpfilename);
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);

	if ( Gn != 3*total_p ) pERROR("reading files wrong.");

	a/=a0;
	b/=a0;
	//===================================================
	// get boundingbox, i.e. max & min of x, y
	float *ymp=mp+total_p, *amp=ymp+total_p;
	float xmin=*mp, xmax=xmin, ymin=*ymp, ymax=ymin;
	for (i=1; i<total_p; i++) {
		if ( xmin > mp[i] ) xmin = mp[i];
		else if ( xmax < mp[i] ) xmax = mp[i];
		if ( ymin > ymp[i] ) ymin = ymp[i];
		else if ( ymax < ymp[i] ) ymax = ymp[i];
	}
	//printf("%%%%BoundingBox: %f %f %f %f\n", xmin, ymin, xmax, ymax);
	//===================================================
	// get max & min of pT and pR
	/*
	double pTmin=*pT, pTmax=pTmin, pRmin=*pR, pRmax=pRmin;
	for (i=1; i<total_p; i++) {
		if ( pTmin > pT[i] ) pTmin = pT[i];
		else if ( pTmax < pT[i] ) pTmax = pT[i];
		if ( pRmin > pR[i] ) pRmin = pR[i];
		else if ( pRmax < pR[i] ) pRmax = pR[i];
	}*/
	// windsorize pT and pR
	double *copy = Malloc(double, total_p); POINTER_NULL(copy);
	double pTmin, pTmax, pRmin, pRmax;
	windsorize_wcp(total_p, pT, copy, 0.001, &pTmin, &pTmax);
	windsorize_wcp(total_p, pR, copy, 0.001, &pRmin, &pRmax);
	printf("OK\n");
	free(copy);
	printf("OK\n");
	
	//=========================================================
	// plot  PT
	char *epsfilename = getfilename(file, "_SPT.eps");
	FILE *fp = fopen(epsfilename, "w");
	FILE_NULL(fp, epsfilename);

	char string0 []="%!PS-Adobe-2.0 EPSF-2.0\n%%Title:";
	fprintf(fp, "%s %s\n", string0, epsfilename);
	time_t rawtime;
	time (&rawtime);
	struct tm *timeinfo=localtime( &rawtime );
	fprintf(fp, "%%%%Creator: e_softspot\n"
			    "%%%%CreationDate: %s", asctime (timeinfo));
	// it seems that asctime has a implicate \n
	fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", 
			    xmin-a, ymin-a, 
				xmax+a, ymax+a);
	fprintf(fp, "%% mode number : %d\n%% threshold : %f\n", 
			Gn-N, threshold/sum);
	char string1 []="\
%%EndComments\n\n\
%%BeginProlog\n\
/n {newpath} bind def\n\
/m {moveto} bind def\n\
/l {lineto} bind def\n\
/s {stroke} bind def\n\
/R {sethsbcolor} bind def\n\
/EllipseDict 8 dict def\n\
EllipseDict /mtrx matrix put\n\
/Ellipse {\n\
  /orientation exch def\n\
  /yrad exch def\n\
  /xrad exch def\n\
  /y exch def\n\
  /x exch def\n\
  /savematrix mtrx currentmatrix def\n\
  x y translate orientation rotate xrad yrad scale 0 0 1 0 360 arc\n\
  closepath fill\n\
  savematrix setmatrix\n\
} def\n\
%%EndProlog\n\
\n\
/EllipseBegin {EllipseDict begin /EllipseEnteredState save def} def\n\
/EllipseEnd {EllipseEnteredState restore end} def\n\
\n\
save\n\
EllipseBegin\n\
gsave\n";
	fprintf(fp, "%s", string1);

	float r2d = 180.0/acos(-1.0);
	double color, blue=4.0/6.0;
	//=============================================
	// plot ellipses & pT
	for (i=0; i<total_p; i++) {
		// make it the larger, the redder
		color=(1.0-(pT[i]-pTmin)/(pTmax-pTmin))*blue;
		fprintf(fp, "%2.2f 1 1 R n %f %f %f %f %f Ellipse\n",
				color, mp[i], ymp[i], a, b, r2d*amp[i]);
	}
	fprintf(fp, "grestore\nEllipseEnd\n"
				"showpage\n%%%%EOF");
	Fclose(fp, epsfilename);
	free(epsfilename);


	//=========================================================
	// plot  PT
	epsfilename = getfilename(file, "_SPR.eps");
	fp = fopen(epsfilename, "w");
	FILE_NULL(fp, epsfilename);

	//char string0 []="%!PS-Adobe-2.0 EPSF-2.0\n%%Title:";
	fprintf(fp, "%s %s\n", string0, epsfilename);
	time (&rawtime);
	timeinfo=localtime( &rawtime );
	fprintf(fp, "%%%%Creator: e_softspot\n"
			    "%%%%CreationDate: %s", asctime (timeinfo));
	// it seems that asctime has a implicate \n
	fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", 
			    xmin-a, ymin-a, 
				xmax+a, ymax+a);
	fprintf(fp, "%% mode number : %d\n%% threshold : %f\n", 
			Gn-N, threshold/sum);
	/*
	char string1 []="\
%%EndComments\n\n\
%%BeginProlog\n\
/n {newpath} bind def\n\
/m {moveto} bind def\n\
/l {lineto} bind def\n\
/s {stroke} bind def\n\
/R {sethsbcolor} bind def\n\
/EllipseDict 8 dict def\n\
EllipseDict /mtrx matrix put\n\
/Ellipse {\n\
  /orientation exch def\n\
  /yrad exch def\n\
  /xrad exch def\n\
  /y exch def\n\
  /x exch def\n\
  /savematrix mtrx currentmatrix def\n\
  x y translate orientation rotate xrad yrad scale 0 0 1 0 360 arc\n\
  closepath fill\n\
  savematrix setmatrix\n\
} def\n\
%%EndProlog\n\
\n\
/EllipseBegin {EllipseDict begin /EllipseEnteredState save def} def\n\
/EllipseEnd {EllipseEnteredState restore end} def\n\
\n\
save\n\
EllipseBegin\n\
gsave\n";
*/
	fprintf(fp, "%s", string1);

	//=============================================
	// plot ellipses & pR
	for (i=0; i<total_p; i++) {
		color=(1.0-(pR[i]-pRmin)/(pRmax-pRmin))*blue;
		fprintf(fp, "%2.2f 1 1 R n %f %f %f %f %f Ellipse\n",
				color, mp[i], ymp[i], a, b, r2d*amp[i]);
	}
	fprintf(fp, "grestore\nEllipseEnd\n"
				"showpage\n%%%%EOF");
	Fclose(fp, epsfilename);
	free(epsfilename);
	free(mp); free(pT); free(pR);
}
