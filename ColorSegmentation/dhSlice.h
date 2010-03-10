#ifndef _dhSLICE_H_
#define _dhSLICE_H_

#include "dhImage.h"
// RLI Slices are currently not scalable.

namespace dh
{

class Col_Slice;
class RGB_Slice;
class RL_Slice;
class RI_Slice;
class LI_Slice;
class Slice3D_Space;

enum AxSign { POS = 0, NEG = 1 };
const int axis_dir[2] = { 1, -1 };
const IntensityType axis_min[2] = { 0, 255 };
const IntensityType axis_max[2] = { 255, 0 };
const int axis_thickness = 2;

class Col_Slice : public PPM_Image
{ 
public:
	const int scale;
	const int size;
		
	int overflow_count;

	Col_Slice( int scale_in, _RGB bkgrnd )
    : PPM_Image( 256*scale_in + 2, 256*scale_in + 2, bkgrnd ),
			  scale(scale_in), size(256*scale + 2),
			  overflow_count(0)
	{
	}
		
	virtual RC_coords color_coords( _RGB clr ) = 0;
		
	void plot_point ( _RGB clr, _RGB false_clr = RGB_BLACK );
	void plot_cross ( _RGB clr, _RGB false_clr = RGB_BLACK, int size = 1, 
		                  bool border = false, _RGB border_clr = RGB_GRAY(1) );

	void inc_point( int r, int c, ColorAxis ax, IntensityType cc, int incr );
	void dec_point( int r, int c, ColorAxis ax, IntensityType cc, int incr );

	void inc_freq_at ( int r, int c, int incr );
	void inc_freq( _RGB clr, int incr );

	_RGB point_color( _RGB clr );	
	void cut_intensity( float factor, _RGB clr );

	void print_overflow();
};

class RGB_Slice : public Col_Slice
{ 
public:
	const ColorAxis proj_axis, x_axis, y_axis;
	const IntensityType slice_at;
	const AxSign x_sign, y_sign;

	RGB_Slice( int scale_in, _RGB bkgrnd,
		                     ColorAxis proj_axis_in, int slice_at_in,
		                     ColorAxis x_axis_in, AxSign x_sign_in,
		                     ColorAxis y_axis_in, AxSign y_sign_in,
									bool draw_border = true )
    : Col_Slice( scale_in, bkgrnd ),
			  proj_axis(proj_axis_in), slice_at(slice_at_in),
			  x_axis(x_axis_in), x_sign(x_sign_in),
			  y_axis(y_axis_in), y_sign(y_sign_in)
	{
		if ( draw_border )
		{ 
			draw_frame();
		}
	}

	void draw_frame();
    RC_coords color_coords( _RGB clr );			 
};

class Slice3D_Space
{ 
public:
	Col_Slice* s1;
	Col_Slice* s2;
	Col_Slice* s3;
		
	Slice3D_Space( Col_Slice* s1_in, Col_Slice* s2_in, Col_Slice* s3_in )
	: s1(s1_in), s2(s2_in), s3(s3_in)
	{
	}
		
	inline void plot_point( _RGB clr, _RGB false_clr = RGB_BLACK )
	{ 
		s1->plot_point (clr, false_clr);
		s2->plot_point (clr, false_clr);
		s3->plot_point (clr, false_clr);
	}
		 
	inline void plot_cross ( _RGB clr, _RGB false_clr = RGB_BLACK, int size = 1, 
							bool border = false, _RGB border_clr = RGB_GRAY(1) )
	{ 
		s1->plot_cross (clr, false_clr, size, border, border_clr);
		s2->plot_cross (clr, false_clr, size, border, border_clr);
		s3->plot_cross (clr, false_clr, size, border, border_clr);
	}
		
	inline void inc_freq ( _RGB clr, int incr )
	{ 
		s1->inc_freq(clr, incr);
		s2->inc_freq(clr, incr);
		s3->inc_freq(clr, incr);
	}
		
	inline void cut_intensity( float factor, _RGB clr )
	{ 
		s1->cut_intensity(factor, clr);
		s2->cut_intensity(factor, clr);
		s3->cut_intensity(factor, clr);
	}
		
	inline void print_overflow()
	{ 
		s1->print_overflow();
		s2->print_overflow();
		s3->print_overflow();
	}	
};

class RL_Slice : public RGB_Slice
{ 
public:
 
    RC_coords HSI_color_coords ( _RGB clr )
	{ 
		RLI RLI_clr = (RLI)clr;
		return( RC_coords( scale * (255 - RLI_clr.L), scale * RLI_clr.R + 1 ) );
	}

    void plot_circle_point ( _RGB clr )
	{ 
		set( HSI_color_coords( clr ), clr );
	}
		
    void draw_circle()
	{
		int j;

		_RGB color = _RGB(255, 0, 0);
		for ( j=0; j<255; j++ )
		{ 
			plot_circle_point( color );
			color.G++;
		}
		for ( j=0; j<255; j++ )
		{ 
			plot_circle_point( color );
			color.R--;
		}
		for ( j=0; j<255; j++ )
		{ 
			plot_circle_point( color );
			color.B++;
		}
		for ( j=0; j<255; j++ )
		{ 
			plot_circle_point( color );
			color.G--;
		}
		for ( j=0; j<255; j++ )
		{ 
			plot_circle_point( color );
			color.R++;
		}
		for ( j=0; j<256; j++ )
		{ 
			plot_circle_point( color );
			color.B--;
		}
		// Fill in left sixth: G-B diagonal
		color = _RGB(0, 255, 0);
		for ( j=0; j<256; j++ )
		{ 
			plot_circle_point( color );
			color.B++; color.G--;
		}
	}

	RL_Slice ( )
		: RGB_Slice( 1, _RGB(0, 0, 0), aB, 128, aR, POS, aG, NEG, false )
	{
	}
};


class RI_Slice : public RGB_Slice
{ 
public:
	void draw_axes ()
	{ 
		// This will only work for standard size!!  
		int i, r, c;
		// Draw grey axis
		for ( i = 0, c = 126 + axis_thickness / 2; i < axis_thickness; i++, c++ )
        { 
			for ( r = 1; r < 257; r++ )
            { 
				if ( P(r, c) == RGB_BLACK)
	            { 
					set( r, c, RGB_GRAY( 256 - r ) ); 
				}
            }
		}
			
		// Draw red axis
		for ( i = 0, r = 126 + axis_thickness / 2; i < axis_thickness; i++, r++ )
	    { 
			for ( c = 1; c < 256; c++ )
            { 
				 if ( P(r, c) == RGB_BLACK)
		         { 
					 set( r, c, _RGB ( RLI ( c - 1, 127, 255/3 + 255/3 * (255 - ( c - 1 )) / 255 ) ) );
				 }
			}
	     }
	}

    RI_Slice ( )
	: RGB_Slice( 1, _RGB(0, 0, 0), aG, 128, aR, POS, aB, NEG, false )
	{
	}
};

class LI_Slice : public RGB_Slice
{ 
public:
	void draw_axes ()
	{ 
		// This will only work for standard size!!  
		int i, r, c;
		// Draw grey axis
		for ( i = 0, c = 126 + axis_thickness / 2; i < axis_thickness; i++, c++ )
	    { 
			for ( r = 1; r < 257; r++ )
            { 
				if ( P(r, c) == RGB_BLACK)
	            { 
					set( r, c, RGB_GRAY( 256 - r ) ); 
				}
	        }
        }

        // Draw lime axis
		for ( i = 0, r = 126 + axis_thickness / 2; i < axis_thickness; i++, r++ )
      	{ 
			for ( c = 1; c < 256; c++ )
            { 
				if ( P(r, c) == RGB_BLACK)
		        { 
					set( r, c, _RGB ( RLI ( 127, c - 1, 127) ) ); 
				}
			}
		}
	}

	LI_Slice ( )
	: RGB_Slice( 1, _RGB(0, 0, 0), aR, 128, aG, POS, aB, NEG, false )
	{
	}
};

} // end namespace dh
#endif
