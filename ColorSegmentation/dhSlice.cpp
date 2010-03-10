#include "dhSlice.h"
// RLI Slices are currently not scalable.

namespace dh
{

	void Col_Slice::plot_point ( _RGB clr, _RGB false_clr )
	{ 
		if( false_clr == RGB_BLACK )
		{ 
			set( color_coords( clr ), clr ); 
		}
		else
		{ 
			set( color_coords( clr ), false_clr); 
		}
	}

	void Col_Slice::plot_cross ( _RGB clr, _RGB false_clr , int size , 
		                  bool border, _RGB border_clr )
	{ 
		RC_coords coords = color_coords( clr );
		_RGB plot_clr = ( false_clr == RGB_BLACK ) ? clr : false_clr ;
        set( coords, plot_clr );
		for(int i = 1; i <= size + 1; i++ )
		{ 
			set( coords + RC_coords( i, 0), plot_clr );
			if (border)
			{ 
				set( coords + RC_coords( i, 1), border_clr );
				set( coords + RC_coords( i, -1), border_clr );
			}
				
			set( coords + RC_coords(-i, 0), plot_clr );
			if (border)
			{ 
				set( coords + RC_coords(-i, 1), border_clr );
				set( coords + RC_coords(-i, -1), border_clr );
            }
            
			set( coords + RC_coords( 0, i), plot_clr );
			if (border)
			{
				set( coords + RC_coords( 1, i), border_clr );
				set( coords + RC_coords(-1, i), border_clr );
			}
			  
			set( coords + RC_coords( 0,-i), plot_clr );
			if (border)
			{ 
				set( coords + RC_coords( 1,-i), border_clr );
				set( coords + RC_coords(-1,-i), border_clr );
			}
		}
			
		if (border)
		{ 
			set( coords + RC_coords( size+1, 0), border_clr );
			set( coords + RC_coords( -(size+1), 0), border_clr );
			set( coords + RC_coords( 0, size+1), border_clr );
			set( coords + RC_coords( 0, -(size+1)), border_clr );
        }
	}

	void Col_Slice::inc_point( int r, int c, ColorAxis ax, IntensityType cc, int incr )
	{ 
		int newc = (int)cc + incr;
		if (newc > 255)
		{ 
			set ( r, c, ax, 255 );
			inc_freq_at ( r, c, newc - 255 );
		}
		else
		{ 
			set ( r, c, ax, newc );
		}
	}

	void Col_Slice::dec_point( int r, int c, ColorAxis ax, IntensityType cc, int incr )
	{ 
		int newc = (int)cc - incr;
		if (newc < 0)
		{	
			set ( r, c, ax, 0 );
			inc_freq_at ( r, c, -newc );
		}
		else
		{ 
			set ( r, c, ax, newc );
		}
	}

	void Col_Slice::inc_freq_at ( int r, int c, int incr )
	{ 
		_RGB cc = P(r, c);
			
		if (cc.B == 255)
		{ 
			if (cc.R == 255)
			{ 
				if (cc.G == 255)
				{ 
					overflow_count++; 
				}
				else
				{ 
					inc_point ( r, c, aG, cc.G, incr ); 
				}
			}
			else
			{ 
				if (cc.G == 0)
				{ 
					inc_point ( r, c, aR, cc.R, incr ); 
				}
				else
				{ 
					dec_point ( r, c, aG, cc.G, incr ); 
				}
			}
		}
		else
		{ 
			if (cc.R == 0)
			{ 
				if (cc.G == 0)
				{ 
					set ( r, c, aR, 255 ); 
				}
				else
				{ 
					inc_point ( r, c, aB, cc.B, incr ); 
				}
			}
			else
			{ 
				if (cc.G == 255)
				{ 
					dec_point ( r, c, aR, cc.R, incr ); 
				}
				else
				{ 
					inc_point ( r, c, aG, cc.G, incr ); 
				}
			}
		}
	}

	void Col_Slice::inc_freq( _RGB clr, int incr )
	{ 
		RC_coords coords = color_coords( clr );
		inc_freq_at ( coords.r, coords.c, incr );
	}

	_RGB Col_Slice::point_color( _RGB clr )
	{ 
		return ( P( color_coords( clr ) ) );
	}
		
	void Col_Slice::cut_intensity( float factor, _RGB clr )
	{ 
		plot_point ( clr, point_color( clr ) / factor );
	}

	void Col_Slice::print_overflow()
	{ 
		if (overflow_count > 0)
		{ std::cout << "    White points overflowed "
			        << overflow_count << " times." << std::endl;
		}
	}

	void RGB_Slice::draw_frame()
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

	RC_coords RGB_Slice::color_coords( _RGB clr )
	{ 
		return ( RC_coords (
							(y_sign == NEG ? 255 - clr.AxisColor(y_axis) :
			                              clr.AxisColor(y_axis) ) + 1,
							(x_sign == NEG ? 255 - clr.AxisColor(x_axis) :
			                              clr.AxisColor(x_axis) ) + 1
					        ) );
	}	

} // end namespace dh