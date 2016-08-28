#include "Magick.h"

#include "colloid_base.h"

#include <Magick++.h>
#include <sstream>

using namespace Magick;

void readimage(colloid_base& cb, const char * filename, const bool * not_in_use)
{
	Image image;
	image.read(filename);

	size_t col=image.columns(), row=image.rows();
	
	if (not_in_use)
		cb.reserve_memory(col, row);
	float * cb_ptr=cb.get_array_pointer();

	int n=col*row;
	unsigned char * pixels=new unsigned char [n];
	image.write(0,0, col, row, "I", CharPixel, pixels);
	
	for (int i=0; i<n; i++)
		*(cb_ptr++)=float(*(pixels+i));

	delete [] pixels;
}


void showposition(colloid_base& p, const int factor, const char * filename)
{
	//const int factor=1;
	int * size=p.get_size();
	float *pointer=p.get_array_pointer();

	float *pp=pointer;
	int i;
	float xmin=*pp, xmax=*pp, ymin=*(pp+1), ymax=*(pp+1);
	pp+=size[0];
	for (i=1; i<size[1]; i++)
	{
		xmin=( xmin <= *pp ) ? xmin : *pp;
		xmax=( xmax >= *pp ) ? xmax : *pp;
		ymin=( ymin <= *(pp+1) ) ? ymin : *(pp+1);
		ymax=( ymax >= *(pp+1) ) ? ymax : *(pp+1);
		pp+=size[0];
	}

	int col=int(xmax-xmin+1.5)*factor;
	int row=int(ymax-ymin+1.5)*factor;

	//#ifdef DEBUG
	//cout << col << 'x' << row << endl;
	//#endif

	//
	Image image(Geometry(col+100, row+100), "white");
	image.type(TrueColorType);
	image.modifyImage();

	//pp=pointer;
	//for (i=0; i<size[1]; i++)
	//{
	//	image.pixelColor(int(*pp-xmin+0.5+50), int(*(pp+1)-ymin+0.5+50),"red");
	//	pp+=size[0];
	//}
	Pixels view(image);

	PixelPacket *pixels=view.get(50, 50, col, row);
	pp=pointer;
	Color red("red");
	//int x, y;
	for (i=0; i<size[1]; i++)
	{
		//x=int(*pp-xmin+0.5);
		//y=int(*(pp+1)-ymin+0.5);
		//cout << x << ' ' << y << endl;
		*(pixels+(int(*pp-xmin+0.5)+int(*(pp+1)-ymin+0.5)*col)*factor)=red;
		pp+=size[0];
	}
	//cout << i << endl;
	view.sync();
	//
	/* A new method not debuged
	char * pixels=new char [col*row];
	for (i=0; i<col*row; i++)
		pixels[i]=255;

	pp=pointer;
	for (i=0; i<size[1]; i++)
	{
		pixels[int(*pp-xmin+0.5)*col+int(*(pp+1)-ymin+0.5)]=0;
		pp+=size[0];
	}
	Image image;
	image.read(col, row, "I", CharPixel, pixels);
	*/

	// draw axis
	
	Color black("black");
	for(i=0; i<col; i++)
	{
		*(pixels+i*factor)=black;
		*(pixels+(i+col)*factor)=black;
		*(pixels+(i+(row-2)*col)*factor)=black;
		*(pixels+(i+(row-1)*col)*factor)=black;
	}
	for(i=0; i<row; i++)
	{
		*(pixels+i*col*factor)=black;
		*(pixels+(i*col+1)*factor)=black;
		*(pixels+(i*col+(col-2))*factor)=black;
		*(pixels+(i*col+(col-1))*factor)=black;
	}

	view.sync();

	// add text
	image.strokeColor("black");
	image.fillColor("black");
	image.strokeWidth(1);
	image.font("");
	image.fontPointsize(20);

	std::stringstream ss;
	ss << int(xmin);
	image.draw( DrawableText( 50, 40, ss.str()));
	ss.str("");
	ss << int(xmax);
	image.draw( DrawableText( col+20, 40, ss.str()));

	image.fontPointsize(20);
	ss.str("");
	ss << int(ymin);
	image.draw( DrawableText( 5, 70, ss.str()));
	ss.str("");
	ss << int(ymax);
	image.draw( DrawableText( 10, row+50, ss.str()));

	//cout << "sync\n";
	image.write(filename);
	//cout << "write\n";

	delete [] size;
	//delete [] pixels;
}
