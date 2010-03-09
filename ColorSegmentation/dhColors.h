#ifndef _dhCOLORS_H_
#define _dhCOLORS_H_

#include <math.h>
#include <iostream>

namespace dh
{
	
#ifndef M_PI
#define M_PI 3.14159265358979
#endif

class _RGB;
class HSI;
class RLI;
class XYZ;

#define RGB_GRAY(level) _RGB( level, level, level )
#define HSI_GRAY(level) HSI( 0, 0, level )

enum ColorAxis { aR = 0, aG = 1, aB = 2 };
typedef char IntensityType; //For RGB and RLI space


//=======================================================================
//===  RGB (Red - Green - Blue) COLOR CLASS
//=======================================================================
class _RGB
{ 
public:
	IntensityType R;
	IntensityType G;
    IntensityType B;

	_RGB(IntensityType r, IntensityType g, IntensityType b) 
		: R(r), G(g), B(b) {};	

	_RGB() : R(0), G(0), B(0) {};
		
	IntensityType AxisColor(ColorAxis ax);

	friend int operator==(const _RGB& c1, const _RGB& c2);
	friend int operator!=(const _RGB& c1, const _RGB& c2);

	friend _RGB operator*(const _RGB& c, double scale_factor);
	friend _RGB operator/(const _RGB& c, double scale_factor);
      
	//friend ostream &operator<<(ostream &, const RGB &);

	operator HSI() const;
	operator RLI() const;
	operator XYZ() const;
		
	inline RLI mapRGBtoRLI();
 };

//=======================================================================
//===  HSI (Hue - Saturation - Intensity) COLOR CLASS
//=======================================================================
class HSI
{ 
public:
    float H;   // Hue: Radians from red, +angle = +frequency
    float S;   // Scale 0 to 1.0
    float I;   // Scale 0 to 1.0

    HSI (float h, float s, float i) : H(h), S(s), I(i) {};
	HSI () : H(0), S(0), I(0) {};
	
	friend int operator==(const HSI& c1, const HSI& c2);
	friend int operator!=(const HSI& c1, const HSI& c2);
		
	friend HSI operator*(const HSI& c, double scale_factor);
	friend HSI operator/(const HSI& c, double scale_factor);
		
	//friend ostream &operator<<(ostream &, const HSI &);
   
	operator _RGB() const;
	operator RLI() const;
 };

//=======================================================================
//===  RLI (Red - Lime - Intensity) COLOR CLASS
//=======================================================================
class RLI
{ 
public:
    IntensityType R;    // All coordinates are on scale 0 - 255.
    IntensityType L;
    IntensityType I;

    RLI (IntensityType r, IntensityType l, IntensityType i) : R(r), L(l), I(i) {};
    RLI () : R(0), L(0), I(0) {};

	friend int operator==(const RLI& c1, const RLI& c2);
	friend int operator!=(const RLI& c1, const RLI& c2);

	friend HSI operator*(const RLI& c, double scale_factor);
	friend RLI operator/(const RLI& c, double scale_factor);

    //friend ostream &operator<<(ostream &, const RLI &);
   
    operator HSI() const;
	operator _RGB() const;

    inline _RGB mapRLItoRGB();
 };

inline _RGB::operator RLI() const { return( (RLI)((HSI)(*this)) ); }
inline RLI::operator _RGB() const { return( (_RGB)((HSI)(*this)) ); }

//=======================================================================
//===  XYZ GENERIC SIGNED 3-D COORDINATE CLASS
//=======================================================================
class XYZ
{ 
public:
    double X;
    double Y;
    double Z;
        
    XYZ (double x, double y, double z) : X(x), Y(y), Z(z) {};	
	XYZ () : X(0), Y(0), Z(0) {};
		
	friend int operator==(const XYZ& c1, const XYZ& c2);
	friend int operator!=(const XYZ& c1, const XYZ& c2);

    friend XYZ operator+(const XYZ& c1, const XYZ& c2);
	friend XYZ operator-(const XYZ& c1, const XYZ& c2);
	friend XYZ operator*(const XYZ& c, const double scale_factor);
	friend XYZ operator/(const XYZ& c, const double scale_factor);
      
	//friend ostream &operator<<(ostream &, const RGB &);
      
	XYZ& operator=(XYZ);
	XYZ& operator=(_RGB);
	operator _RGB() const;
 };


//=======================================================================
// Color Constants
//=======================================================================
const IntensityType BLACK = 0;
const IntensityType WHITE = (IntensityType)255;

const _RGB RGB_BLACK = RGB_GRAY (0);
const _RGB RGB_WHITE = RGB_GRAY ((IntensityType)255);
const HSI HSI_BLACK = HSI_GRAY (0);
const HSI HSI_WHITE = HSI_GRAY (1);

//=======================================================================

inline bool is_gray( _RGB c ) { return (bool)( c.R == c.G && c.R == c.B ); }
inline bool is_gray( HSI c ) { return (bool)( c.S == 0 ); }
inline bool is_gray( RLI c ) { return (bool)( c.R == 0  && c.L == 0 ); }

double dot (const XYZ& c1, const XYZ& c2);
XYZ cross (const XYZ& c1, const XYZ& c2);
double magnitude (const XYZ& c);
} // end namespace dh

#endif
