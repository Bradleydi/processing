#include "normalmode.h"
#include "io.h"
#include "miscellaneous.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>

void e_plot_config(const char *file, int n, double a, double b,
		bool Tracked)
{
	char *gdffile = getfilename (file, ".gdf");
	colloid_base p;
	readgdf(p, gdffile);

	float *p_ptr=p.get_array_pointer();
	int *size=p.get_size();
	int ti=size[0]-1;
	if (Tracked) --ti;
	if ( ti < 3 )
		pERROR("data should contain x y angle.\n");
	int i, t;
	int si=size[0]-ti;

	float *pointer=p_ptr;
	float xmin=*pointer, xmax=xmin, ymin=*(pointer+1), ymax=ymin;
	float t1=*(pointer+ti), t2=t1;
	for (i=1; i<size[1]; i++) {
		if ( xmin > *pointer ) xmin = *pointer;
		else if ( xmax < *pointer ) xmax = *pointer;
		
		if ( ymin > *(pointer+1) ) ymin = *(pointer+1);
		else if ( ymax < *(pointer+1) ) ymax = *(pointer+1);

		pointer+=ti;
		if ( t1 > *pointer ) t1 = *pointer;
		else if ( t2 < *pointer ) t2 = *pointer;
		pointer+=si;
	}
	const int total_t = (int)(t2-t1) + 1;
	const int tmin = (int)t1;
	printf("# total frame               = %d\n", total_t);
	if (tmin!=0) printf("# \tframe index [ %d : %d ]\n", tmin, (int)t2);

	if ( n < tmin || n > (int)t2 ) pERROR("input wrong frame index");

	//============================================
	// plot the data
	char subfix [20];
	sprintf(subfix, "_%05d.eps", n);
	char *epsfile=getfilename(file, subfix);
	FILE *fp = fopen(epsfile, "w");
	FILE_NULL(fp, epsfile);
	char string0 []="%!PS-Adobe-2.0 EPSF-2.0\n%%Title:";
	fprintf(fp, "%s %s\n", string0, epsfile);
	time_t rawtime;
	time (&rawtime);
	struct tm *timeinfo=localtime( &rawtime );
	fprintf(fp, "%%%%Creator: e_plot_config\n"
			    "%%%%CreationDate: %s", asctime (timeinfo));
	// it seems that asctime has a implicate \n
	fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", 
			    xmin-a, ymin-a, xmax+a, ymax+a);
	char string1 []="\
%%EndComments\n\
%%BeginProlog\n\
/n {newpath} bind def\n\
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
EllipseBegin\n\
0.6 setgray\n";
	fprintf(fp, "%s\n", string1);

	pointer=p_ptr;
	float r2d = 180.0/acos(-1.0);
	for (i=0; i<size[1]; i++) {
		t=(int)(*(pointer+ti)); // t
		if (t==n)
			fprintf(fp, "n %f %f %f %f %f Ellipse\n",
					*(pointer), *(pointer+1), a, b, r2d*(*(pointer+2)));
					//x, y, a, b, ang);
		pointer += size[0];
	}

	fprintf(fp, "EllipseEnd\nshowpage\n%%%%EOF");
	Fclose(fp, epsfile);
	free(epsfile);
	delete [] size;
}
