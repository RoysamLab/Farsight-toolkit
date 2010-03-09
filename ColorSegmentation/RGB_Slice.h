#ifndef _RGB_SLICE_H_
#define _RGB_SLICE_H

#include "Col_Slice.h"

namespace dh
{

class RGB_Slice : public Col_Slice
{ 
public:
	const ColorAxis proj_axis, x_axis, y_axis;
	const IntensityType slice_at;
	const AxSign x_sign, y_sign;

	void draw_frame()
	{   
		int r, c, i, x, y;
		IntensityType dv;
		int rcmax = size - 1;
			
		set(0, 0, proj_axis, slice_at);
		set(0, 0, x_axis, axis_min[x_sign]);
		set(0, 0, y_axis, axis_min[y_sign]);

		set(0, rcmax, proj_axis, slice_at);
		set(0, rcmax, x_axis, axis_max[x_sign]);
		set(0, rcmax, y_axis, axis_min[y_sign]);

		set(rcmax, rcmax, proj_axis, slice_at);
		set(rcmax, rcmax, x_axis, axis_max[x_sign]);
		set(rcmax, rcmax, y_axis, axis_max[y_sign]);

		set(rcmax, 0, proj_axis, slice_at);
		set(rcmax, 0, x_axis, axis_min[x_sign]);
		set(rcmax, 0, y_axis, axis_max[y_sign]);
		  
		dv = axis_min[x_sign];
		for( c = 0; c < 256; c++ )
		{ 
			for ( i = 1; i < scale+1; i++ )
			{ 
				y = c*scale+i;
				set ( 0, y, proj_axis, slice_at );
				set ( 0, y, x_axis, dv );
				set ( 0, y, y_axis, axis_min[y_sign] );
					
				set ( rcmax, y, proj_axis, slice_at );
				set ( rcmax, y, x_axis, dv );
				set ( rcmax, y, y_axis, axis_max[y_sign] );

			}
			dv = dv + axis_dir[x_sign]; 
		}
			
		dv = axis_min[y_sign];			
		for( r = 0; r < 256; r++ )
		{ 
			for ( i = 1; i < scale+1; i++ )
			{ 
				x = r*scale+i;
				set ( x, 0, proj_axis, slice_at );
				set ( x, 0, x_axis, axis_min[x_sign] );
				set ( x, 0, y_axis, dv );

				set ( x, rcmax, proj_axis, slice_at );
				set ( x, rcmax, x_axis, axis_max[x_sign] );
				set ( x, rcmax, y_axis, dv );
			}
			dv = dv + axis_dir[y_sign]; 
		}
	}

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

    RC_coords color_coords( _RGB clr )
	{ 
		return ( RC_coords (
							(y_sign == NEG ? 255 - clr.AxisColor(y_axis) :
			                              clr.AxisColor(y_axis) ) + 1,
							(x_sign == NEG ? 255 - clr.AxisColor(x_axis) :
			                              clr.AxisColor(x_axis) ) + 1
					        ) );
	}								 
};

} //end namespace dh
#endif
