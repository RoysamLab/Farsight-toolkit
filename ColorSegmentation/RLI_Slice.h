#ifndef _RLI_SLICE_H_
#define _RLI_SLICE_H_

#include "RGB_Slice.h"
// RLI Slices are currently not scalable.

namespace dh
{

const int axis_thickness = 2;

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
