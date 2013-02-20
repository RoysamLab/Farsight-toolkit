#ifndef _dhCLASSIFIERS_H_
#define _dhCLASSIFIERS_H_

#include <math.h>
#include "dhColors.h"

namespace dh
{

class Classifier
{ 
public:
	static float euclidean_dist
		(_RGB c1, _RGB c2, float rw = 1.0, float gw = 1.0, float bw = 1.0)
	{
		float rd = ( c2.R - c1.R ) * rw;
		float gd = ( c2.G - c1.G ) * gw;
		float bd = ( c2.B - c1.B ) * bw;
		return( sqrt( rd*rd + gd*gd + bd*bd ) );
	}

	static float euclidean_dist
	     (RLI c1, RLI c2, float rw = 1.0, float lw = 1.0, float iw = 1.0)
	{
		  float rd = ( c2.R - c1.R ) * rw;
		  float ld = ( c2.L - c1.L ) * lw;
		  float id = ( c2.I - c1.I ) * iw;
		  return( sqrt( rd*rd + ld*ld + id*id ) );
	}

	static float ortho_euclidean_dist
	     (const HSI c1, const HSI c2, float hw = 1.0, float sw = 1.0, float iw = 1.0)
	{
		  float hd = ( c2.H - c1.H ) * hw;
		  float sd = ( c2.S - c1.S ) * sw;
		  float id = ( c2.I - c1.I ) * iw;
		  return( sqrt( hd*hd + sd*sd + id*id ) );
	}
	  
	static float true_euclidean_dist
	     (const HSI c1, const HSI c2, float hsw = 1.0, float iw = 1.0)
	{
		  float xd = c2.S * cos(c2.H) - c1.S * cos(c1.H);
		  float yd = c2.S * sin(c2.H) - c1.S * sin(c1.H);
		  float xyd = hsw * sqrt( xd*xd + yd*yd );
		  float id = iw * ( c2.I - c1.I );
		  return( sqrt( xyd*xyd + id*id ) );
	}
};

} // end namespace dh
#endif
