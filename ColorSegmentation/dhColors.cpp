#include "dhColors.h"

namespace dh
{

//=======================================================================
//===  RGB (Red - Green - Blue) COLOR CLASS
//=======================================================================
RGBType _RGB::AxisColor(RGBType ax)
{
	switch (ax)
	{
		case 0:
			return(R);
		case 1:
			return(G);
		case 2 :
			return(B);
	}
	std::cerr<<"RGB.axis_color: Invalid Color Axis." << std::endl;
	return((RGBType)255);
}

int operator==(const _RGB& c1, const _RGB& c2)
 { return ((c1.R == c2.R && c1.G == c2.G && c1.B == c2.B) ? 1 : 0);}

int operator!=(const _RGB& c1, const _RGB& c2)
 { return ((c1.R != c2.R || c1.G != c2.G || c1.B != c2.B) ? 1 : 0);}

_RGB operator*(const _RGB& c, double scale_factor)
{ 
	int new_R = (int)(scale_factor * c.R);
	int new_G = (int)(scale_factor * c.G);
	int new_B = (int)(scale_factor * c.B);
	if (    new_R < 0 || new_R > 255
	     || new_G < 0 || new_G > 255
		 || new_B < 0 || new_B > 255 )
	{ 
		std::cerr<<"RGB::operator*(RGB, double): result out of range." << std::endl;
		return ( _RGB( (RGBType)255, (RGBType)255, (RGBType)255 ) );
	}
	else
	{ 
		return ( _RGB((RGBType)new_R, (RGBType)new_G, (RGBType)new_B) ); 
	}
}

_RGB operator/(const _RGB& c, double scale_factor)
 { return ( c * ( 1 / scale_factor ) ); }


std::ostream &operator<<(std::ostream &output, const _RGB &c)
 { output << "[RGB: " << (int)c.R << " " << (int)c.G << " " << (int)c.B << "]";
   return output;
 }

/*
istream &operator>>(istream &input, _RGB &c)
 { Str inbuf;
   input >> inbuf;
	
	if( !input || inbuf != "[RGB:" )
	 { std::cerr<<"RGB Input does not have a valid header!"; }
	else
	 { int r_in = 0; int g_in = 0; int b_in = 0;
	   input >> r_in >> g_in >> b_in;
		c.R = r_in;
		c.G = g_in;
		c.B = b_in;
		 cout << r_in << ";"<< g_in << ";" << b_in << endl;
		 if ( !input 
		     || r_in < 0 || r_in > 255
			  || g_in < 0 || g_in > 255
			  || b_in < 0 || b_in > 255 )
		 { std::cerr<<"RGB Input has invalid data!"; }
		else
		 {
		   input.ignore( 1, ']' );
		   if ( !input )
		    { std::cerr<<"RGB Input has no closing ']'!"; }
		 }
	 }
	return input;
 }
*/

//=======================================================================
//===  HSI (Hue - Saturation - Intensity) COLOR CLASS
//=======================================================================

int operator==(const HSI& c1, const HSI& c2)
 { return (c1.H == c2.H && c1.S == c2.S && c1.I == c2.I);}

int operator!=(const HSI& c1, const HSI& c2)
 { return (c1.H != c2.H || c1.S != c2.S || c1.I != c2.I);}
		
HSI operator*(const HSI& c, double scale_factor)
 { float new_I = (float)(scale_factor * c.I);
	if ( new_I < 0 || new_I > 1.0 )
	{ std::cerr<<"HSI::operator*(HSI, double): result out of range." << std::endl;
	   return ( HSI( 0, 0, 1 ) );
	 }
   else
	 { return ( HSI(c.H, c.S, new_I) ); }
 }

HSI operator/(const HSI& c, double scale_factor)
 { return ( c * ( 1 / scale_factor ) ); }


std::ostream &operator<<(std::ostream &output, const HSI &c)
 { output << "[HSI: " << (int)((c.H * 180 / M_PI)+.5) << " deg  "
          // using precision(3) here changes the value from 0 to 60 for white!?!
          << c.S << "  " << c.I << "]";
	return output;
 }



//=======================================================================
//===  RLI (Red - Lime - Intensity) COLOR CLASS
//=======================================================================
int operator==(const RLI& c1, const RLI& c2)
 { return (c1.R == c2.R && c1.L == c2.L && c1.I == c2.I);}

int operator!=(const RLI& c1, const RLI& c2)
 { return (c1.R != c2.R || c1.L != c2.L || c1.I != c2.I);}

/* RLI double should be like RGB double ?
RLI operator*(const RLI& c, double scale_factor)
 { 
	 int new_I = (int)(scale_factor * c.I);
	if ( 
		new_I < 0 || new_I > 255 )
	 { std::cerr<<"RLI::operator*(RLI, double): result out of range.";
	   return ( RLI( (IntensityType)0, (IntensityType)0, (IntensityType)255 ) );
	 }
   else
	 { return ( RLI(c.R, c.L, (float)new_I) ); }
 }
*/
RLI operator*(const RLI& c, double scale_factor)
{ 
	int new_R = (int)(scale_factor * c.R);
	int new_L = (int)(scale_factor * c.L);
	int new_I = (int)(scale_factor * c.I);
	if (    new_R < 0 || new_R > 255
	     || new_L < 0 || new_L > 255
		 || new_I < 0 || new_I > 255 )
	{ 
		std::cerr<<"RLI::operator*(RLI, double): result out of range." << std::endl;
		return ( RLI( (RLIType)255, (RLIType)255, (RLIType)255 ) );
	}
	else
	{ 
		return ( RLI((RLIType)new_R, (RLIType)new_L, (RLIType)new_I) ); 
	}
}

RLI operator/(const RLI& c, double scale_factor)
 { return ( c * ( 1 / scale_factor ) ); }



std::ostream &operator<<(std::ostream &output, const RLI &c)
 { output << "[RLI: " << (int)c.R << " " << (int)c.L << " " << (int)c.I << "]";
	return output;
 }



//=======================================================================
//===  XYZ GENERIC SIGNED 3-D COORDINATE CLASS
//=======================================================================
int operator==(const XYZ& c1, const XYZ& c2)
 { return ((c1.X == c2.X && c1.Y == c2.Y && c1.Z == c2.Z) ? 1 : 0);}

int operator!=(const XYZ& c1, const XYZ& c2)
 { return ((c1.X != c2.X || c1.Y != c2.Y || c1.Z != c2.Z) ? 1 : 0);}

XYZ operator+(const XYZ& c1, const XYZ& c2)
 { return ( XYZ( c1.X + c2.X, c1.Y + c2.Y, c1.Z + c2.Z ) );
 }

XYZ operator-(const XYZ& c1, const XYZ& c2)
 { return ( XYZ( c1.X - c2.X, c1.Y - c2.Y, c1.Z - c2.Z ) );
 }

XYZ operator*(const XYZ& c, const double scale_factor)
 { return ( XYZ( c.X * scale_factor,
                 c.Y * scale_factor,
					  c.Z * scale_factor ) );
 }

XYZ operator/(const XYZ& c, const double scale_factor)
 { return ( c * ( 1 / scale_factor ) ); }

//-------- XYZ as a vector: vector operations --------
double dot (const XYZ& c1, const XYZ& c2)
 { return ( (c1.X * c2.X) + (c1.Y * c2.Y) + (c1.Z * c2.Z) ); }
 
XYZ cross (const XYZ& c1, const XYZ& c2)
 { return ( XYZ( (c1.Y * c2.Z) - (c1.Z * c2.Y),
                 (c1.Z * c2.X) - (c1.X * c2.Z),
					  (c1.X * c2.Y) - (c1.Y * c2.X) ) );
 }

double magnitude (const XYZ& c)
 { return ( sqrt ( (c.X * c.X) + (c.Y * c.Y) + (c.Z * c.Z) ) ); }

std::ostream &operator<<(std::ostream &output, const XYZ &c)
 { output << "[XYZ: " << c.X << " " << c.Y << " " << c.Z << "]";
   return output;
 }

XYZ& XYZ::operator=(XYZ c)
 { X = c.X;
   Y = c.Y;
	Z = c.Z;
   return *this;
 }

XYZ& XYZ::operator=(_RGB c)
 { X = c.R;
   Y = c.G;
	Z = c.B;
   return *this;
 }

XYZ::operator _RGB() const
{ 
	if (    X < 0 || X > 255
	     || Y < 0 || Y > 255
		  || Z < 0 || Z > 255 )
	{ 
		std::cerr<<"XYZ::operator RGB(): XYZ value [XYZ " << X << " " << Y << " " << Z << "] is not in range." << std::endl;
		return ( _RGB( (RGBType)255, (RGBType)255, (RGBType)255 ) );
	}
	else
	{ 
		return ( _RGB( (RGBType)X, (RGBType)Y, (RGBType)Z ) ); 
	}
}

//=======================================================================
//== CONVERSION OPERATORS
//=======================================================================

const float athird = (float)(1.0/3.0);
const float sqrtthree = (float)sqrt( 3.0 );
const float sqrt2o3 = (float)sqrt(2.0/3.0);

/* OLD HSI CONVERTION
_RGB::operator HSI() const
 {
  float total = (float)R + (float)G + (float)B;
  float i = athird * total / 255;
	 
  if ( is_gray(*this) )
   { return( HSI_GRAY(i) );
   }
  else
   {
 	  float r = (float)R / total;
	  float g = (float)G / total;
	  float b = (float)B / total;

	  float s = 1 - 3 * ((r<g) ? ((r<b) ? r : b ) : ((g<b) ? g : b ));

	  float h;
	  float rmt = r - athird;
	  float gmt = g - athird;
	  float bmt = b - athird;

	  float at = sqrt( rmt*rmt + gmt*gmt + bmt*bmt );
	  float bt = athird * (rmt + rmt - gmt - bmt);
     
	  float swt = (bt / at / sqrt2o3);
	  if ( swt > 1 ) swt = 1;
	  if ( swt < -1 ) swt = -1;
	  h = acos( swt );

	  if(b > g)
	   { h = (float)(2*M_PI - h); }

	  return ( HSI ( h, s, i ) );
   }
 }
*/
// BETTER IMPLEMENTATION:??
_RGB::operator HSI() const
 {
  float total = (float)R + (float)G + (float)B;
  float i = athird * total / 255;
	 
  if ( is_gray(*this) )
   { 
	   return( HSI_GRAY(i) );
   }
  else
   {
 	  float r = (float)R / total;
	  float g = (float)G / total;
	  float b = (float)B / total;

	  float s = 1 - 3 * ((r<g) ? ((r<b) ? r : b ) : ((g<b) ? g : b ));

	  double angle = 0.5 * ( (r-g) + (r-b) ) / (sqrt((r-g)*(r-g)+(r-b)*(g-b)) );
	
	  double h;
	  if(b<=g)
	  {
		h = acos(angle);
	  }
	  else
	  {
		h = 2*M_PI - acos(angle);
	  }

	  return ( HSI ( (HSIType)h, s, i ) );
   }
 }
 

_RGB::operator XYZ() const
 { return ( XYZ ( R, G, B ) );
 }
 
HSI::operator _RGB() const
 { double r, g, b;
	if ( H > 0 && H <= 2 * M_PI / 3.0 )
	 {	b = athird * ( 1 - S );
	   r = athird * ( 1 + ( S * cos(H) ) / cos( M_PI / 3.0 - H ) );
		g = 1 - (r + b);
	 }
	else if ( H > 2 * M_PI / 3.0 && H <= 4 * M_PI / 3.0 )
	 { double h = H - 2 * M_PI / 3.0;
	   r = athird * ( 1 - S );
	   g = athird * ( 1 + ( S * cos(h) ) / cos( M_PI / 3.0 - h ) );
		b = 1 - (r + g);
	 }
	else
	 { double h = H - 4 * M_PI / 3.0;
	   g = athird * ( 1 - S );
	   b = athird * ( 1 + ( S * cos(h) ) / cos( M_PI / 3.0 - h ) );
		r = 1 - (g + b);
	 }
	float sf = 3 * I * 255;
	int ri = (int)(sf*r +.5);
	int gi = (int)(sf*g +.5);
	int bi = (int)(sf*b +.5);
	if(   ri < 0   || bi < 0   || gi < 0
	   || ri > 255 || bi > 255 || gi > 255 )
	 { std::cout<<"WARNING: HSI::operator RGB(): returned values out of range!!! ["\
			<< ri << " " << gi << " " << bi << "]" << std::endl;
	   return ( _RGB(  ( ri < 0 ? 0 : ( ri > 255 ? 255 : ri )),
		               ( gi < 0 ? 0 : ( gi > 255 ? 255 : gi )),
					      ( bi < 0 ? 0 : ( bi > 255 ? 255 : bi )) ) );
	 }
	return ( _RGB( ri, gi, bi ));
 }

HSI::operator RLI() const
 {
   float cr = cos ( H );
   float cl = sin ( H );

	float l, r;

   #if HSI_CYLINDER
	   l = S * cl;
		r = S * cr;
	#else
	   float icpt = S * (1 - 2 * fabs( .5 - I )) / sqrtthree;
	
	 	// Something is a little off here... the corners don't match up.
	   if ( H >= 0 && H < 2 * M_PI / 3.0 )
	 	 {	l = icpt / ( 1 + ( cr / ( sqrtthree * cl ) ) );
	 		r = ( cr / cl ) * l;
	 	 }
	 	else if ( H >= 2 * M_PI / 3.0 && H < 4 * M_PI / 3.0 )
	 	 { r = - icpt;
	 	   l = - cl * icpt / cr;
	 	 }
	 	else
	 	 { l = icpt / ( -1 + ( cr / ( sqrtthree * cl ) ) );
	 		r = ( cr / cl ) * l;
	 	 }
	#endif
	return ( RLI( (int)(r * 127.0 + 127),
	              (int)(l * 127.0 + 127), 
	              (int)(I * 255 + .5) ) );
 }

RLI::operator HSI() const
{ 
	#if HSI_CYLINDER
		double nr = ( R - 127 ) / 127.0;
		double nl = ( L - 127 ) / 127.0;
		
		double h;
		if ( nr == 0 )  // On l-axis
		{ 
			 h = (nl >= 0) ? M_PI / 2 : 3 * M_PI / 2;
		}
		else
		{ 
			 h = (nr >= 0) ?  ((nl >= 0) ? atan( nl / nr ) : atan( nl / nr ) + 2 * M_PI)
		                  : atan( nl / nr ) + M_PI;
		}

		HSIType s = (HSIType)sqrt( nr*nr + nl*nl );
		HSIType i = (HSIType)(I / 255.0);
      
		return ( HSI( (HSIType)h, s, i ) );
	#else
	   std::cerr<<"The RLI -> HSI transform is not implemented for conical space.";
	#endif
 }

RLI::operator XYZ() const
 { return ( XYZ ( R, L, I ) );
 }

} //end namespace dh
