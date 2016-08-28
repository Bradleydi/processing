#include "libtiff.h"

#include <cstdio>

using namespace std;

void readtiff(colloid_tiff& tiff, char *filename)
{
	TIFF *tif = TIFFOpen (filename, "r");
	if (tif==0)
	{
		fprintf(stderr, "# Error: cannot open tiff file %s\n", filename);
		exit (1);
	}

	uint32 w, h;

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);

	size_t npixels=w*h;
	uint32* raster=(uint32 *) _TIFFmalloc(npixels *sizeof(uint32));
	if (raster==NULL)
	{
		fprintf(stderr, "# Error: raster==NULL\n");
		exit (1);
	}

	if (TIFFReadRGBAImage(tif, width, height, raster, 0))
	{
	}
	else
	{
		fprintf(stderr, "# Error: can not access pixels\n");
		exit (1);
	}
	_TIFFfree(raster);
	TIFFClose(tif);
}
