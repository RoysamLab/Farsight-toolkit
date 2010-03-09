#ifndef _PPM_IMAGE_H_
#define _PPM_IMAGE_H_

#include "Color_Image.h"

namespace dh
{

#define FOR_EACH_AXIS for ( ax=0; ax<=2; ax++)

class PPM_Image : public Color_Image
{
public:    
      
	PPM_Image() : mem(NULL) { rows = 0;  cols = 0; }
				
	inline PPM_Image(int in_rows, int in_cols, _RGB color = _RGB(0, 0, 0) )
		 { define( in_rows, in_cols, color ); }

	~PPM_Image();
		
	inline _RGB P (int r, int c);
      
	//This should be defined in Color_Image
	inline _RGB P (RC_coords coords);
	inline IntensityType P (int r, int c, ColorAxis axis);
	inline IntensityType P (RC_coords coords, ColorAxis axis);
	inline void set(int r, int c, _RGB color);
  	//This should be defined in Color_Image
	inline void set(RC_coords coords, _RGB color);
	inline void set(int r, int c, ColorAxis axis, IntensityType value);
	inline void set(RC_coords coords, ColorAxis axis, IntensityType value);
		
	inline int empty() { return (mem == NULL); }
		
protected:
	void allocate_image();
	IntensityType *** mem;
 
};

//------------------------------------
//  INLINE FUNCTIONS
//------------------------------------

inline _RGB PPM_Image::P (int r, int c)
{ 
	if ( r < 0 || r >= rows || c < 0 || c >= cols )
	{ 
		std::cerr<<"PPM_Image.P: Invalid Coordinates";
		 report_coords(r, c);
	}
	return ( _RGB(mem[aR][r][c], mem[aG][r][c], mem[aB][r][c]) );
}

//This should be defined in Color_Image
inline _RGB PPM_Image::P (RC_coords coords)
 { return ( P( coords.r, coords.c ) ); }


inline IntensityType PPM_Image::P (int r, int c, ColorAxis axis)
{ 
	if ( r < 0 || r >= rows || c < 0 || c >= cols )
	{ 
		std::cerr<<"PPM_Image.P: Invalid Coordinates";
		report_coords(r, c);
	}
	if ( axis < 0 || axis > 2 )
		std::cerr<<"PPM_Image.P: Invalid Color Axis";
	return ( mem[axis][r][c] );
}

//This should be defined in Color_Image
inline IntensityType PPM_Image::P (RC_coords coords, ColorAxis axis)
 { return ( P( coords.r, coords.c, axis ) ); }
 
inline void PPM_Image::set(int r, int c, _RGB color)
{ 
	if ( r < 0 || r >= rows || c < 0 || c >= cols )
	{ 
		std::cerr<<"PPM_Image.set(r, c, color): Invalid Coordinates";
		report_coords(r, c);
	}
	mem[aR][r][c] = color.R;
	mem[aG][r][c] = color.G;
	mem[aB][r][c] = color.B;
}
    
//This should be defined in Color_Image
inline void PPM_Image::set(RC_coords coords, _RGB color)
 { set( coords.r, coords.c,  color ); }

inline void PPM_Image::set(int r, int c, ColorAxis axis, IntensityType value)
{ 
	if ( r < 0 || r >= rows || c < 0 || c >= cols )
	{ 
		std::cerr<<"PPM_Image.set(r, c, axis, intensity): Invalid Coordinates";
		report_coords(r, c);
	}
	if ( axis < 0 || axis > 2 )
		std::cerr<<"PPM_Image.set(r, c, axis, intensity): Invalid Color Axis";
	mem[axis][r][c] = value;
}

//This should be defined in Color_Image
inline void PPM_Image::set(RC_coords coords, ColorAxis axis, IntensityType value)
 { set ( coords.r, coords.c, axis, value); } 


} //end namespace dh
#endif
