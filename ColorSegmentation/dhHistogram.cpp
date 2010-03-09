#include "dhHistogram.h"

namespace dh
{

void RGBHistogram::dump()
{ 
	int i, j, k;
	char inp[20];
	int planeinc = 1;
	for ( i = 0; i < 128; i += planeinc )
	{ 
		for ( j = 127; j >= 0; j-- )
		{ 
			std::cout << std::endl;
			FOR_AXIS(k)
			{ 
				if (a[i][j][k] == 0)
			    { 
					std::cout << "."; 
				}
				else
				{ 
					int n = (int)floor( (double)log( (double)a[i][j][k] ) / (double)log((double)2) );
					if ( n>=0 && n<=9 )
					{ std::cout << n; }
					else
					{ std::cout << '#'; }
				}
			}
		}
		std::cout << std::endl << "PLANE RED = " << i << "  ";
		std::cin >> inp;
		if ( strlen(inp) != 0 )
		{ planeinc = atoi( inp ); }
	}
}

void RGBHistogram::smooth()
{ 
	//std::cout << "Smoothing..." << flush;
	int sk[3][3][3] =	{ { { 3, 4, 3 }, 
	 						{ 4, 6, 4 }, 
							{ 3, 4, 3 } }, 
							 
						  { { 4, 6, 4 }, 
							{ 6,12, 6 }, 
							{ 4, 6, 4 } }, 
							
						  { { 3, 4, 3 }, 
							{ 4, 6, 4 }, 
							{ 3, 4, 3 } } };

	const int sk_sum = 120;
	int i, j, k, x, y, ixp, jyp;
	long sum;
   
	processing_array.clear_to(0);

	//-------------------------------------------------------------------
	for ( i = rmin; i <= rmax; i++ )
	 { //std::cout << "." << flush;
	   for ( j = gmin; j <= gmax; j++ )
		 { for ( k = bmin; k < bmax; k++ )
			 { sum = 0;
			   for ( x = 0; x < 3; x++ )
				 { ixp = i + x - 1;
				   for ( y = 0; y < 3; y++ )
					 { jyp = j + y - 1;
					   //for ( z = 0; z < 3; z++ )
						// { sum += sk[x][y][z] *
						//			 a[ixp][jyp][k+z-1];
						// }
						sum +=
						     sk[x][y][0] * a[ixp][jyp][k-1]
							+ sk[x][y][1] * a[ixp][jyp][k]
							+ sk[x][y][2] * a[ixp][jyp][k+1];
					 }
				 }
			   sa[i][j][k] = sum / sk_sum;
			 }
		 }
	 }
   std::cout << std::endl;
  
   FOR_AXIS(i)
	 { FOR_AXIS(j)
		 { FOR_AXIS(k)
			 { a[i][j][k] = sa[i][j][k];
			   if ( a[i][j][k] > max_freq )
				 { max_freq = a[i][j][k];
				   mode = _RGB(i, j, k);
				 }
			 }
		 }
	 }
}
 
void RGBHistogram::find_bounding_box()
 { 
   std::cout << "Finding bounding box..." << std::endl;
	
	int i, j, k;
	const int minhits1 = 24;    // MAGIC NUMBER
	const int minhits2 = 24;    // MAGIC NUMBER
	const int minhits3 = 16;    // MAGIC NUMBER
	//-------------------------------------------------------------------
	// Efficiency trick: find 3-D bounding box for i, j, k such that
	// there are at least minhits hits in each dimention
	// Do red, then blue, then green.
   //Find rmin
	int th = 0;
	for( i = 1; i < 127; i++ )
	 { for ( j = 1; j < 127; j++ )
	    { for ( k = 1; k < 127; k++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits1 ) goto rmindone;
				 }
			 }
		 }
	 }
	rmindone:
	rmin = max( i-1, 1 );
	// Find rmax
	th = 0;
	for( i = 126; i > 0; i-- )
	 { for ( j = 1; j < 127; j++ )
	    { for ( k = 1; k < 127; k++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits1 ) goto rmaxdone;
				 }
			 }
		 }
	 }
	rmaxdone:
	rmax = min( i+1, 126 );

	// Find bmin
	th = 0;
	for( k = 1; k < 127; k++ )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( j = 1; j < 127; j++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits2 ) goto bmindone;
				 }
			 }
		 }
	 }
	bmindone:
	bmin = max( k-1, 1 );

	// Find bmax
	th = 0;
	for( k = 126; k > 0; k-- )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( j = 1; j < 127; j++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits2 ) goto bmaxdone;
				 }
			 }
		 }
	 }
	bmaxdone:
	bmax = min( k+1, 126 );

	// Find gmin
	th = 0;
	for( j = 1; j < 127; j++ )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( k = bmin; k <= bmax; k++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits3 ) goto gmindone;
				 }
			 }
		 }
	 }
	gmindone:
	gmin = max( j-1, 1 );

	// Find jmax
	th = 0;
	for( j = 126; j > 0; j-- )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( k = bmin; k <= bmax; k++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits3 ) goto gmaxdone;
				 }
			 }
		 }
	 }
	gmaxdone:
	gmax = min( j+1, 126 );
 }

void RGBHistogram::mark_point_and_nbrs(IntensityType x, IntensityType y, IntensityType z)
 { if ( sa[x][y][z] == 0 && a[x][y][z] > 0 )
    { sa[x][y][z] = 1;
	   
	   mark_point_and_nbrs(x+1, y, z);
		mark_point_and_nbrs(x-1, y, z);
		mark_point_and_nbrs(x, y+1, z);
		mark_point_and_nbrs(x, y-1, z);
		mark_point_and_nbrs(x, y, z+1);
		mark_point_and_nbrs(x, y, z-1);
	 }
 }

void RGBHistogram::delete_secondary_blobs()
 { 
   // This will remove all parts of the histogram that are not
	// 6-connected (in 3-d) with the blob containing the mode.

	std::cout << "Removing secondary blobs..." << std::endl;
	
   processing_array.clear_to(0);
	mark_point_and_nbrs(mode.R, mode.G, mode.B);

	int i, j, k;
	FOR_AXIS(i)
	 { FOR_AXIS(j)
		 { FOR_AXIS(k)
			 { if (sa[i][j][k] == 0)
			    { a[i][j][k] = 0; }
			 }
		 }
	 }
 }

} //end namespace dh
