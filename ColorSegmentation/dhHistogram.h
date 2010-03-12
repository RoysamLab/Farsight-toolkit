#ifndef _dhHistogram_H
#define _dhHistogram_H

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dhColors.h"

namespace dh
{



static const int histSize = 128;
#define FOR_AXIS(X) for ( X = 0; X < histSize; X++ )

// Note: The histogram is 128 x 128 x 128, but still uses "RGB" classes
// to store the coordinates.  This is conceptually iffy, since RGB's
// actually go to 256 everywhere alse, and the values are often represent
// RLI values instead of RGB values.

class Array3D
{
public:
	int d1, d2, d3;
	long int*** a;

	Array3D(int d1_in, int d2_in, int d3_in, int init = 1, long int init_value = 0);
	~Array3D();

	void clear_to(long int val);
};

class Histogram
{	
public:
	int max_freq;
	long int mode[3];
	const int inc_scale;
		
	int rmax, rmin,
		gmax, gmin,
		bmax, bmin;

	_RGB modeAsRGB(){return _RGB((RGBType)mode[0],(RGBType)mode[1],(RGBType)mode[2]);};
	RLI modeAsRLI(){return RLI((RGBType)mode[0],(RGBType)mode[1],(RGBType)mode[2]);};

	Histogram(int inc_scale_in = 1);

	void find_bounding_box();
    void smooth();
	void delete_secondary_blobs();
	void dump();

	long int v(int d1, int d2, int d3, bool suppress_warning = false ) const;	
	void set(int d1, int d2, int d3, long int val);
	void inc_element(int d1, int d2, int d3);

	long int D1_proj_at(int d2, int d3);
	long int D2_proj_at(int d1, int d3);
	long int D3_proj_at(int d1, int d3);
	
	void inc( _RGB c )
	{ 
		inc_element( c.R / 2, c.G / 2, c.B / 2 ); 
	};

	void inc( RLI c )
	{
		inc_element( c.R / 2, c.L / 2, c.I / 2 );
	};

protected:
	Array3D histogram_array;
	long int ***a;	//The main array pointer
	Array3D processing_array;
	long int ***sa; //The processing array

	bool check_bounds(int d1, int d2, int d3, bool suppress_warning = false) const;
	void mark_point_and_nbrs(int d1, int d2, int d3);

 };

inline int min( int a, int b )
 { return( ( a <= b ) ? a : b ); }
 
inline int max( int a, int b )
 { return( ( a >= b ) ? a : b ); }

inline int min( int a, int b, int c )
 { return( ( a <= b ) ? min( a, c ) : min( b, c ) ); }

inline int max( int a, int b, int c )
 { return( ( a >= b ) ? max( a, c ) : max( b, c ) ); }

}//end namespace dh
#endif
