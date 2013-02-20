#include "dhEvalState.h"

namespace dh
{

double eval_state(const _RGB& p1, const _RGB& p2, const _RGB& p3)
{
	switch(EVAL_FUNCTION)
	{
	case 1:
		return min_inv_square_dist_eval(p1,p2,p3);
		break;
	case 2:
		return max_perimeter_eval(p1,p2,p3);
		break;
	case 3:
		return area_eval(p1,p2,p3);
		break;
	default:
		return 0.0;
	}
	
}


//////////////////// EVALUATION FUNCTIONS ///////////////////////////
double min_inv_square_dist_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 )
{ 
	double d12 = Classifier::euclidean_dist( p1, p2 );
	double d23 = Classifier::euclidean_dist( p2, p3 );
	double d13 = Classifier::euclidean_dist( p1, p3 );
	if( d12 == 0 || d23 == 0 || d13 == 0 )
	{ 
		return NEG_HUGE; 
	}
	else
	{ 
		return -(   ( 1 / pow ( d12, 2 ) )  
                + ( 1 / pow ( d23, 2 ) )
	             + ( 1 / pow ( d13, 2 ) )
	           );
    }
}

double max_perimeter_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 )
{ 
	 return (   Classifier::euclidean_dist( p1, p2 )
            + Classifier::euclidean_dist( p2, p3 )
				+ Classifier::euclidean_dist( p1, p3 ) );
}

double area_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 )
{ 
	XYZ a = p1;
	XYZ b = p2;
	XYZ c = p3;
	
	if ( a == c )
	 { return (0); }

	XYZ ac = c - a;
	XYZ ab = b - a;
	
	XYZ d = a + ( ac * ( dot ( ac, ab ) / dot ( ac, ac ) ) );
   
	//double t = dot ( ac, ab ) / dot ( ac, ac );

	double area = .5 * magnitude ( ac ) * magnitude ( d - b );
	
	return ( area );
}

} //end namespace dh
