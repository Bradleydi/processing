#include "normalmode.h"

#include "miscellaneous.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <ctime>

/* freedom_t fdm: freedom to plot, translation, or rotation */
void e_plot_mode4(const char *file, double a, double b, double a0, int n, 
		freedom_t fdm, double factor)
{
	if (fdm) // rotation
		e_plot_mode_rotat(file, a, b, a0, n, factor);
	else // translation
		e_plot_mode_trans(file, a, b, a0, n, factor);
}

void e_plot_mode_trans(const char *file, double a, double b, double a0, int n, 
		double factor)
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
	float *ymp=mp+total_p, *Omp=mp+total_p2;
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
	double *yG=G+total_p, *OG=G+total_p2;
	double max_eT=sqrt(G[0]*G[0]+yG[0]*yG[0]), ei;
	double max_eO=*OG, min_eO=*OG;
	double eL=sqrt(G[0]*G[0]+yG[0]*yG[0]);
	for (i=1; i<total_p; i++) {
		ei=sqrt(G[i]*G[i]+yG[i]*yG[i]);
		if ( max_eT < ei )
			max_eT = ei;
		eL += ei;
		ei=OG[i];
		if (ei > max_eO ) max_eO=ei;
		else if (ei < min_eO ) min_eO=ei;
	}
	if (-min_eO > max_eO) max_eO=-min_eO;
	const double max_array_length=3.0*a;
	double xyratio1=max_array_length/max_eT;
	
	// By practice I found that the idea of setting max length is not good,
	// since in some case, most are almost same, then all will be too large
	// and in some case, one is very large, and cannot be shown
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
	xyratio *= factor;


	//=========================================================
	// plot as eps
	char modesubfix [128];
	sprintf(modesubfix, "_modeT_%05d.eps", n);
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

	double PI = acos(-1.0);
	float r2d = 180.0/PI;
	double color, lighten=0.4;  // the smaller lighten, the lighter
	//=============================================
	// plot ellipses, no orientational mode.
	// the color represent the orientation of ellipses
	for (i=0; i<total_p; i++) {
		//color=OG[i]/2.0/max_eO+0.5;
		color = acos(cos(2.0*Omp[i]))/PI;
		//fprintf(stderr, "%d %f\n", i, color);
		if (color<0.0 || color > 1.0) pERROR("Wrong mean orientation\n");
		fprintf(fp, "%2.2f %2.2f %2.2f R n %f %f %f %f %f Ellipse\n",
				1.0-lighten*color, 1.0-lighten, 1.0-lighten*(1.0-color),
				mp[i], ymp[i], a*0.85, b*0.85, 
				r2d*(Omp[i]+0.5*factor*(OG[i]/max_eO))); // rotation
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
	fprintf(fp, "gsave\n%2.2f setlinewidth\n0 0 0 setrgbcolor\n", lw);
	for (i=0; i<total_p; i++) {
		if ( longarrow[i] )
		{
		ei=sqrt(G[i]*G[i]+yG[i]*yG[i]);
		Cos=G[i]/ei;
		Sin=yG[i]/ei;
		x2=mp[i]+G[i]*xyratio; // xend
		y2=ymp[i]+yG[i]*xyratio; // yend
		if (xyratio<0) { Cos=-Cos; Sin=-Sin; }
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


void e_plot_mode_rotat(const char *file, double a, double b, double a0, int n, 
		double factor)
{
}
