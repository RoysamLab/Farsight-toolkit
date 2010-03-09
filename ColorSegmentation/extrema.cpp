#ifndef _EXTREMA_CXX_
#define _EXTREMA_CXX_

inline int min( int a, int b )
 { return( ( a <= b ) ? a : b ); }
 
inline int max( int a, int b )
 { return( ( a >= b ) ? a : b ); }

inline int min( int a, int b, int c )
 { return( ( a <= b ) ? min( a, c ) : min( b, c ) ); }

inline int max( int a, int b, int c )
 { return( ( a >= b ) ? max( a, c ) : max( b, c ) ); }

#endif
