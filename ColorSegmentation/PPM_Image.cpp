#include "PPM_Image.h"

namespace dh
{

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
