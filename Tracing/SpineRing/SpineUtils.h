#ifndef SPINEUTILS_H
#define SPINEUTILS_H

#include <vnl/vnl_math.h>
// the following function defs are defined in vnl_math.h
// if vnl_math is not desired, simply implement them 
// as inline functions here
#define MAX   vnl_math_max
#define MIN   vnl_math_min
#define EDIST(x1,x2) sqrt(vnl_vector_ssd((x1), (x2)))

#define NUMMAX(x1, x2)  ((x1)>=(x2)? (x1) : (x2))
#define NUMMIN(x1, x2)  ((x1)>(x2)? (x2) : (x1))
////#define DEBUGPRINT PrintSelf()

#define VERBOSE(x) std::cout<<x<<std::endl;
#define VERBOSE2(x,y) std::cout<<x<<y<<std::endl;
#define VERBOSE2NL(x,y) std::cout<<x<<y;
#endif
