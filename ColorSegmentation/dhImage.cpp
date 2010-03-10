
#include "dhImage.h"
#include <iostream>

namespace dh
{

//------- Image functions:
void Image::report_coords(int r, int c)
{ 
	std::cout << "		At coordinates [ r = " << r << ", c = " << c << "]";
}

//------- Color Image Functions
void Color_Image::define (int in_rows, int in_cols, _RGB color )
{ 
	rows = in_rows;
	cols = in_cols;

	int r, c;

	allocate_image();
	
	FOR_EACH_ROW
	{ 
		FOR_EACH_COL
		{
			 set(r, c, color);
		}
	 }
}

//------- PPM Image Functions
void PPM_Image::allocate_image()
{
	int r, ax;
 
	mem = new IntensityType** [3];
 
	FOR_EACH_AXIS
	{ 
		mem[ax] = new IntensityType* [rows];
		FOR_EACH_ROW
	    { 
			mem[ax][r] = new IntensityType[cols];
	    }
	}
}

PPM_Image::~PPM_Image() 
{ 
	if ( !empty() )
	{ 
		int r, ax;
		FOR_EACH_AXIS
	    { 
			FOR_EACH_ROW
		    { 
				delete mem[ax][r]; 
			}
			delete mem[ax];
		}
		delete mem;
	}
}




} //end namespace dh

