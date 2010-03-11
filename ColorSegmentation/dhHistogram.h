#ifndef _dhHistogram_H
#define _dhHistogram_H

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dhColors.h"
#include "dhArray3D.h"
#include "extrema.cpp"

namespace dh
{

#define FOR_AXIS(X) for ( X = 0; X < 128; X++ )

// Note: The histogram is 128 x 128 x 128, but still uses "RGB" classes
// to store the coordinates.  This is conceptually iffy, since RGB's
// actually go to 256 everywhere alse, and the values are often represent
// RLI values instead of RGB values.

class Histogram
{	
protected:
	Array3D histogram_array;
	long int*** a;
      
	Array3D processing_array;
    long int*** sa;

	void mark_point_and_nbrs(IntensityType x, IntensityType y, IntensityType z);

public:
	int max_freq;
	//_RGB mode;
	IntensityType mode[3];
	const int inc_scale;
		
	int rmax, rmin,
		gmax, gmin,
		bmax, bmin;

	_RGB modeAsRGB(){return _RGB(mode[0],mode[1],mode[2]);};
	RLI modeAsRLI(){return RLI(mode[0],mode[1],mode[2]);};

	Histogram(int inc_scale_in = 1) 
		:   max_freq(0), 
			inc_scale(inc_scale_in), 
			rmax(126), rmin(1), 
			gmax(126), gmin(1), 
			bmax(126), bmin(1), 
			histogram_array(128, 128, 128),
			processing_array(128, 128, 128, 0)
	{
		mode[0] = 0;
		mode[1] = 0;
		mode[2] = 0;
		a = histogram_array.a;
		sa = processing_array.a;
	}

	void dump();
	void find_bounding_box();
    void smooth();
	void delete_secondary_blobs();

	char v(int r, int g, int b, bool suppress_warning = false ) const
	{ 
		if ( r < 0 || r > 127
		  || g < 0 || g > 127
		  || b < 0 || b > 127 )
		{ 
			if(!suppress_warning)
			{ 
				 std::cout << "The histogram coordinates (" << r << ", " << g << ", " << b << ") are invalid!!" << std::endl; 
			}
			return 0;
		}	
		else
		{ 
			return (char)a[r][g][b]; 
		}
    }
		
	void set(int r, int g, int b, char val)
	{
		if ( r < 0 || r > 127
		  || g < 0 || g > 127
		  || b < 0 || b > 127 )
		{ 
			std::cout << "The histogram coordinates (" << r << ", " << g << ", " << b << ") are invalid!!" << std::endl; 
		}
		else
		{ 
			a[r][g][b] = val; 
		}	
	}

	void inc_element(int r, int g, int b)
	{ 
		if ( r < 0 || r > 127
		  || g < 0 || g > 127
		  || b < 0 || b > 127 )
		{ 
			std::cout << "The histogram coordinates (" << r << ", " << g << ", " << b << ") are invalid!!" << std::endl; 
		}
		else    
		{ 
			if ( a[r][g][b] == 65535 )
			{ 
				std::cerr<<"Histogram increment overflow!!!"; 
			}
			else
			{ 
				a[r][g][b] += inc_scale; 
			}
		}
    }

	void inc( _RGB c )
	{ 
		inc_element( c.R / 2, c.G / 2, c.B / 2 ); 
	}

	void inc( RLI c )
	{
		inc_element(c.R /2, c.L /2, c.I /2 );
	}
		
	long int RB_proj_at ( char r, char b )
	{
		int total = 0;
		char g;
		FOR_AXIS(g)
		{ 
			total += a[r][g][b];
		}
		return(total);
	}

	long int RG_proj_at ( char r, char g )
	{
		int total = 0;
		char b;
		FOR_AXIS(b)
		{ 
			total += a[r][g][b];
		}
		return(total);
	}

	long int GB_proj_at ( char g, char b )
	{ 
		int total = 0;
		char r;
		FOR_AXIS(r)
		{ 
			total += a[r][g][b];
		}
		return(total);
	}
 };

}//end namespace dh
#endif
