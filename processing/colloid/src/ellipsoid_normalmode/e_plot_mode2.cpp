#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

void e_plot_mode2(const char *file, double a, double b, double a0, int n)
{
	if (n<0) pERROR("wrong mode id");
	// get mean position from .mp file
	int total_p;
	float *mp;
	char *mpfilename=getfilename(file, ".mp");
	e_readmp(&total_p, &mp, mpfilename);
	printf("# file \"%s\" loaded!\n", mpfilename); 
	free(mpfilename);

	// check mode n whether valid
	if (n>=3*total_p) {
		fprintf(stderr, "# Error: max mode id = %d\n", 3*total_p-1);
		exit (1);
	}
	
	// get n'th normal mode from .ev file
	char *evfilename=getfilename(file, ".ev");
	int Gn;
	double E;
	double *G;
	//n=total_p*3-1-n;
	e_readmode(total_p*3-1-n, &Gn, &E, &G, evfilename);
	printf("# file \"%s\" loaded!\n", evfilename); 
	free(evfilename);
	E=1.0/sqrt(E);

	if ( Gn != 3*total_p ) pERROR("reading files wrong.");

	a/=a0;
	b/=a0;
	//===================================================
	// get boundingbox, i.e. max & min of x, y
	int i, total_p2=total_p*2;
	float *ymp=mp+total_p, *thetamp=mp+total_p2;
	float xmin=*mp, xmax=xmin, ymin=*ymp, ymax=ymin;
	for (i=1; i<total_p; i++) {
		if ( xmin > mp[i] ) xmin = mp[i];
		else if ( xmax < mp[i] ) xmax = mp[i];
		if ( ymin > ymp[i] ) ymin = ymp[i];
		else if ( ymax < ymp[i] ) ymax = ymp[i];
	}
	//printf("%%%%BoundingBox: %f %f %f %f\n", xmin, ymin, xmax, ymax);
	//===================================================
	// get max abs(e_angle)
	// set max abs(e_angle) as redest color
	// and set average(e_tran)= (1-total(e_angle^2))/N as length a/a0
	// to draw arrows, and as the average arrow length.
	// Maybe set average arrow length is not a good idea, since one may be
	// very long, and so leave the bounding box, and ugly.
	// So set max arrow length as (a/a0)*1.25, and this may be better.
	//
	// it seems that information is too much, so it's better to
	// remove very short arrows to make figure clear.
	double *yG=G+total_p, *thetaG=G+total_p2;
	double max_etrans=sqrt(G[0]*G[0]+yG[0]*yG[0]), ei;
	double max_eag=*thetaG, min_eag=*thetaG;
	double eL=sqrt(G[0]*G[0]+yG[0]*yG[0]);
	for (i=1; i<total_p; i++) {
		ei=sqrt(G[i]*G[i]+yG[i]*yG[i]);
		if ( max_etrans < ei )
			max_etrans = ei;
		eL += ei;
		ei=thetaG[i];
		if (ei > max_eag ) max_eag=ei;
		else if (ei < min_eag ) min_eag=ei;
	}
	if (-min_eag > max_eag) max_eag=-min_eag;
	const double max_array_length=3.0*a;
	double xyratio1=max_array_length/max_etrans;
	
	// By practice I found that the idea of setting max length is not good,
	// since in some case, most are almost same, then all will be too large
	// and in some case, one is very large, and cannot be shown
	/*
	double *yG=G+total_p, *thetaG=G+total_p2;
	double ei, eL=sqrt(G[0]*G[0]+yG[0]*yG[0]);
	double max_eag=*thetaG, min_eag=*thetaG;
	for (i=1; i<total_p; i++) {
		eL += sqrt(G[i]*G[i]+yG[i]*yG[i]);
		ei=thetaG[i];
		if (ei > max_eag ) max_eag=ei;
		else if (ei < min_eag ) min_eag=ei;
	}
	if (-min_eag > max_eag) max_eag=-min_eag;
	*/
	const double ave_array_length=0.5*a;
	double xyratio2=ave_array_length*total_p/eL;
	
	//It seems that use both the max length and average is better
	double xyratio=( xyratio1 < xyratio2 ) ? xyratio1 : xyratio2;
	
	char *longarrow=Calloc(char, total_p); POINTER_NULL(longarrow);
	double arrow_threshold = 0.3*a/xyratio;
	for (i=1; i<total_p; i++) {
		ei=sqrt(G[i]*G[i]+yG[i]*yG[i]);
		if ( ei > arrow_threshold ) longarrow[i]=1;
	}


	//=========================================================
	// plot as eps
	char modesubfix [128];
	sprintf(modesubfix, "_mode_%05d.eps", n);
	char *epsfilename = getfilename(file, modesubfix);
	FILE *fp = fopen(epsfilename, "w");
	FILE_NULL(fp, epsfilename);

	char string0 []="%!PS-Adobe-2.0 EPSF-2.0\n%%Title:";
	fprintf(fp, "%s %s\n", string0, epsfilename);
	time_t rawtime;
	time (&rawtime);
	struct tm *timeinfo=localtime( &rawtime );
	fprintf(fp, "%%%%Creator: e_plot_mode\n"
			    "%%%%CreationDate: %s", asctime (timeinfo));
	// it seems that asctime has a implicate \n
	fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", 
			    xmin-a, ymin-a, 
				xmax+a, ymax+a);
	fprintf(fp, "%% mode ID : %d\n%% freq. :  %f\n", n, E);
	// this for mode ID cannot put before Boundingbox, or the box size
	// will be wrong.
	char string1 []="\
%%EndComments\n\n\
%%BeginProlog\n\
/n {newpath} bind def\n\
/m {moveto} bind def\n\
/l {lineto} bind def\n\
/s {stroke} bind def\n\
/R {setrgbcolor} bind def\n\
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
	double color;
	//=============================================
	// plot ellipses & orientational mode
	for (i=0; i<total_p; i++) {
		color=thetaG[i]/2.0/max_eag+0.5;
		fprintf(fp, "%2.2f 0 %2.2f R n %f %f %f %f %f Ellipse\n",
				color, 1.0-color,
				mp[i], ymp[i], a, b, r2d*thetamp[i]);
					//x, y, a, b, ang);
	}
	fprintf(fp, "grestore\nEllipseEnd\n\n");
	//=============================================
	// plot translational mode
	// arrow head length  &  arrow width
	double L=a/6.0, W=L/sqrt(3.0); // 60 degree arrow
	double Cos, Sin, x2, y2;
	// linewidth
	double lw = b/2.5;
	fprintf(fp, "gsave\n%2.2f setlinewidth\n1 1 0 setrgbcolor\n", lw);
	for (i=0; i<total_p; i++) {
		if ( longarrow[i] )
		{
		ei=sqrt(G[i]*G[i]+yG[i]*yG[i]);
		Cos=G[i]/ei;
		Sin=yG[i]/ei;
		x2=mp[i]+G[i]*xyratio, // xend
		y2=ymp[i]+yG[i]*xyratio, // yend
		/* This method doesn't work, since different angle of lineto
		 * the linejoint is different. So we should use method two
		 * to make it symmetric.
		fprintf(fp, "n %f %f m %f %f l %f %f l %f %f l %f %f l s\n",
				mp[i], ymp[i], // x0 y0
				x2, y2,  // xy end
				x2-L*Cos-W*Sin,   // x_up_arrow
				y2-L*Sin+W*Cos,  
				x2, y2,  // xy end
				x2-L*Cos+W*Sin,   // x_down_arrow
				y2-L*Sin-W*Cos);
				*/
		/* this method also have problems
		fprintf(fp, "n %f %f m %f %f l %f %f l\n%f %f m %f %f l %f %f l s\n",
				mp[i], ymp[i], // x0 y0
				x2, y2,  // xy end
				x2-L*Cos-W*Sin,   // x_up_arrow
				y2-L*Sin+W*Cos,
				mp[i], ymp[i],
				x2, y2,  // xy end
				x2-L*Cos+W*Sin,   // x_down_arrow
				y2-L*Sin-W*Cos);
	*/
		fprintf(fp, "n %f %f m %f %f l %f %f m %f %f l %f %f l s\n",
				mp[i], ymp[i], // x0 y0
				x2, y2,  // xy end
				x2-L*Cos-W*Sin,   // x_up_arrow
				y2-L*Sin+W*Cos,  
				x2, y2,  // xy end
				x2-L*Cos+W*Sin,   // x_down_arrow
				y2-L*Sin-W*Cos);
		}
	}
	fprintf(fp, "grestore\nrestore\nshowpage\n%%%%EOF");
	Fclose(fp, epsfilename);
	free(epsfilename);
	free(mp); free(G);
}
