#ifndef _dhHistogram_H
#define _dhHistogram_H

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConnectedComponentImageFilter.h"

#include "dhColors.h"

namespace dh
{

static const int histSize = 256;
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
	typedef long int PixelType;
	typedef itk::Image< PixelType, 3> HistImageType;

	int max_freq;
	long int mode[3];
	const int inc_scale;
		
	int d1max, d1min, d2max, d2min, d3max, d3min;	//These are the bounding boxes of the histogram (used for convolution)

	_RGB modeAsRGB(){return _RGB((RGBType)mode[0],(RGBType)mode[1],(RGBType)mode[2]);};
	RLI modeAsRLI(){return RLI((RGBType)mode[0],(RGBType)mode[1],(RGBType)mode[2]);};

	Histogram(int inc_scale_in = 1);

	void find_bounding_box();
    void smooth();
	void delete_secondary_blobs();
	void save_as(const char * fname);

	PixelType v(int d1, int d2, int d3) const;
	void set(int d1, int d2, int d3, long int val);
	void inc_element(int d1, int d2, int d3);

	void projection_extrema(int dir, long int &max, long int &min);
	long int proj_at(int dir, int da, int db);
	
	void inc( _RGB c )
	{ 
		inc_element( c.R, c.G, c.B); 
	};

	void inc( RLI c )
	{
		inc_element( c.R, c.L, c.I );
	};

protected:
	HistImageType::Pointer histImage;

	bool check_bounds(int d1, int d2, int d3) const;

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
