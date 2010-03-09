#ifndef _RGB_ATYPE_H_
#define _RGB_ATYPE_H_

#include "dhColors.h"
#include <iostream>
#include <string>

namespace dh
{

class Named_RGB
{ 
public:
	std::string name;
	_RGB clr;
	   
	Named_RGB( std::string name_in = "BLACK", _RGB clr_in = RGB_BLACK )
	: name(name_in), clr(clr_in) {}
};

const int num_named_colors = 11;
const int first_named_hue = 3;

const Named_RGB named_colors[num_named_colors] =
{ 
	Named_RGB( "BLACK", RGB_BLACK ),
	Named_RGB( "WHITE", RGB_WHITE ),
	Named_RGB( "GRAY", RGB_GRAY(128) ),
	Named_RGB( "RED", _RGB(255, 0, 0) ),
	Named_RGB( "ORANGE", _RGB(255, 128, 0) ),
	Named_RGB( "YELLOW", _RGB(255, 255, 0) ),
	Named_RGB( "LIME", _RGB(128, 255, 0) ),
	Named_RGB( "GREEN", _RGB(0, 255, 0) ),
	Named_RGB( "AQUA", _RGB(0, 255, 255) ),
	Named_RGB( "BLUE", _RGB(0, 0, 255) ),
	Named_RGB( "PURPLE", _RGB(255, 0, 255) )
};

class RGB_Atype
{ 
public:
    _RGB red;
	_RGB blue;
	_RGB bkgrnd;
		
	std::string red_name;
	std::string blue_name;
	std::string bkgrnd_name;
		
	_RGB red_approx;
	_RGB blue_approx;
	_RGB bkgrnd_approx;
		
	RGB_Atype ()
	: red(RGB_BLACK), blue (RGB_BLACK), bkgrnd(RGB_BLACK)
	{
	}

	RGB_Atype ( _RGB red_in, _RGB blue_in, _RGB bkgrnd_in, 
	            std::string red_name_in = "RED", std::string blue_name_in = "BLUE", 
				std::string bkgrnd_name_in = "WHITE")
	: red(red_in), blue (blue_in), bkgrnd(bkgrnd_in), 
      red_name(red_name_in), blue_name(blue_name_in),
	  bkgrnd_name(bkgrnd_name_in)
	{ 
		find_named_RGBs();
	}
	
	_RGB find_named_RGB(std::string name);
	void find_named_RGBs();
	void set_named_RGBs_to_indeces( int red_ind, int blue_ind, int bkgrnd_ind );
 };

} //end namespace dh
#endif
	
