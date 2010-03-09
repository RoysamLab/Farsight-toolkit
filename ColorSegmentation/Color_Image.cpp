#include "Color_Image.h"

namespace dh
{

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

} //end namespace dh
