#include "RGB_Atype.h"

namespace dh
{

_RGB RGB_Atype::find_named_RGB(std::string name)
{ 
	for ( int i = 0; i < num_named_colors; i++ )
    {
		if ( name == named_colors[i].name )
	    { 
			return ( named_colors[i].clr ); 
		}
	}
	std::cerr<<"Bad Configuration Parameter:\n " << name << " is not a recognized color.";
	return(RGB_BLACK);
}

void RGB_Atype::find_named_RGBs()
{ 
	red_approx = find_named_RGB(red_name);
	blue_approx = find_named_RGB(blue_name);
	bkgrnd_approx = find_named_RGB(bkgrnd_name);
}

void RGB_Atype::set_named_RGBs_to_indeces( int red_ind, int blue_ind, int bkgrnd_ind )
{ 
	red_name = named_colors[red_ind].name;
	blue_name = named_colors[blue_ind].name;
	bkgrnd_name = named_colors[bkgrnd_ind].name;
		
	red_approx = named_colors[red_ind].clr;
	blue_approx = named_colors[blue_ind].clr;
	bkgrnd_approx = named_colors[bkgrnd_ind].clr;
}

} //end namespace dh
	
