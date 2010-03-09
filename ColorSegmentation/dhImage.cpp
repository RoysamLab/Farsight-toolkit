
#include "dhImage.h"
#include <iostream>

namespace dh
{

void Image::report_coords(int r, int c)
{ 
	std::cout << "		At coordinates [ r = " << r << ", c = " << c << "]";
}

}

