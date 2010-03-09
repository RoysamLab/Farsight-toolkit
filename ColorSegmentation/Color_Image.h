#ifndef _COLOR_IMAGE_H_
#define _COLOR_IMAGE_H_

#include "dhImage.h"
#include "dhColors.h"

namespace dh
{

class Color_Image : public dh::Image
{ 
public:
	virtual void define (int in_rows, int in_cols, _RGB color = RGB_BLACK );
	
	virtual _RGB P (int r, int c) = 0;
	virtual _RGB P (RC_coords coords){ return ( P( coords.r, coords.c ) ); }
		
	virtual IntensityType P (int r, int c, ColorAxis axis) = 0;
		
	inline IntensityType P (RC_coords coords, ColorAxis axis)
		 { return ( P( coords.r, coords.c, axis ) ); }

	virtual void set(int r, int c, _RGB color) = 0;
		
	inline void set(RC_coords coords, _RGB color)
		 { set( coords.r, coords.c,  color ); }
		
    virtual void set(int r, int c, ColorAxis axis, IntensityType value) = 0;
		
	inline void set(RC_coords coords, ColorAxis axis, IntensityType value)
		 { set ( coords.r, coords.c, axis, value); }
 };

} // end namespace dh
		
#endif
