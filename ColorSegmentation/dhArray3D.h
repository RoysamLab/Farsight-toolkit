#ifndef _SHORT_3D_ARRAY_CXX_
#define _SHORT_3D_ARRAY_CXX

namespace dh
{

class Array3D
{ 
public:
	int d1, d2, d3;
	long int*** a;
	
	void clear_to (long int val)
	{ 
		for ( int i=0; i<d1; i++ )
		{ 
			for ( int j=0; j<d2; j++ )
		    { 
				for ( int k=0; k<d3; k++ )
				{ 
					 a[i][j][k] = val;
				}
			}
		 }
	 }

	Array3D(int d1_in, int d2_in, int d3_in,
	               int init = 1, long int init_value = 0)
		: d1 (d1_in), d2(d2_in), d3(d3_in)
	{ 
		int i, j;
	   
		a = new long int**[d1];
		for ( i=0; i<d1; i++ )
		{ 
			a[i] = new long int*[d2];
		}
		
		for ( i=0; i<d1; i++ )
	    { 
			for ( j=0; j<d2; j++ )
		    { 
				a[i][j] = new long int[d3];
			}
		}
      
		if ( init )
		{ 
			clear_to(init_value);
		}
	}
	
	~Array3D()
	{ 
		int i, j;
	   	for ( i=0; i<d1; i++ )
	    { 
			for ( j=0; j<d2; j++ )
		    { 
				delete[] a[i][j];
			}
		 }
		
		for ( i=0; i<d1; i++ )
		{ 
			delete[] a[i];
		}

		delete[] a;
    }
	
};

} //end namespace dh
#endif
