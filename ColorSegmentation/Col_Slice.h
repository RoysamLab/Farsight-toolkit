#ifndef _COL_SLICE_H_
#define _COL_SLICE_H_

#include "PPM_Image.h"

namespace dh
{

enum AxSign { POS = 0, NEG = 1 };

const int axis_dir[2] = { 1, -1 };
const IntensityType axis_min[2] = { 0, 255 };
const IntensityType axis_max[2] = { 255, 0 };

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
		
	void plot_point ( _RGB clr, _RGB false_clr = RGB_BLACK )
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

	void plot_cross ( _RGB clr, _RGB false_clr = RGB_BLACK, int size = 1, 
		                  bool border = false, _RGB border_clr = RGB_GRAY(1) )
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

	void inc_point( int r, int c, ColorAxis ax, IntensityType cc, int incr )
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

	void dec_point( int r, int c, ColorAxis ax, IntensityType cc, int incr )
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

	void inc_freq_at ( int r, int c, int incr )
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

	void inc_freq( _RGB clr, int incr )
	{ 
		RC_coords coords = color_coords( clr );
		inc_freq_at ( coords.r, coords.c, incr );
	}

	_RGB point_color( _RGB clr )
	{ 
		return ( P( color_coords( clr ) ) );
	}
		
	void cut_intensity( float factor, _RGB clr )
	{ 
		plot_point ( clr, point_color( clr ) / factor );
	}

	void print_overflow()
	{ 
		if (overflow_count > 0)
		{ std::cout << "    White points overflowed "
			        << overflow_count << " times." << std::endl;
		}
	}
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



} //end namespace dh
#endif
