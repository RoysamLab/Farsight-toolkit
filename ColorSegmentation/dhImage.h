#ifndef _dhIMAGE_H_
#define _dhIMAGE_H_

namespace dh
{

#define FOR_EACH_ROW for ( r=0; r<rows; r++)
#define FOR_EACH_COL for ( c=0; c<cols; c++)

class Image
{ 
public:
      
	virtual void allocate_image() = 0;

	virtual void report_coords(int r, int c);

	inline int num_rows() { return (rows); }
	inline int num_cols() { return (cols); }
	
	virtual int empty() = 0;
	
protected:
	int rows;
	int cols;
 };

class RC_coords
{ //friend ostream &operator<<(ostream &, const RC_coords &);
   
public:
	int r;
    int c;
      
    RC_coords () : r(0), c(0) {};
    RC_coords ( int in_r, int in_c ) : r(in_r), c(in_c) {};
      
    friend RC_coords operator+(RC_coords, RC_coords);
    friend int operator==(RC_coords, RC_coords);
    friend int operator!=(RC_coords, RC_coords);  
};

//------- Define operators for RC coords ----------
inline RC_coords operator+(RC_coords base, RC_coords delta)
 { return RC_coords( base.r + delta.r, base.c + delta.c ); }
   
inline int operator==(RC_coords a, RC_coords b)
 { return ( a.r == b.r && a.c == b.c ); }

inline int operator!=(RC_coords a, RC_coords b)
 { return ( a.r != b.r || a.c != b.c ); }

} // end namespace dh

#endif
