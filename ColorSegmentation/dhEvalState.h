#ifndef _EVALSTATE_H_
#define _EVALSTATE_H_

#include "dhColors.h"
#include "dhClassifiers.h"

namespace dh
{
	
#define EVAL_FUNCTION 3

const double POS_HUGE = 1E100;
const double NEG_HUGE = -1E100;

//////////////////// EVALUATION FUNCTIONS ///////////////////////////
double min_inv_square_dist_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 );
double max_perimeter_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 );
double area_eval( const _RGB& p1, const _RGB& p2, const _RGB& p3 );

////////////////////////////////////////////////////////////////////
////// Choose default evaluation function here !!! ////////
double eval_state(const _RGB& p1, const _RGB& p2, const _RGB& p3);
////////////////////////////////////////////////////////////////////

} //end namespace dh
#endif
