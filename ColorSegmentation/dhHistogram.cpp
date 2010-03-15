 #include "dhHistogram.h"

namespace dh
{

Histogram::Histogram(int inc_scale_in) 
	:max_freq(0), 
	inc_scale(inc_scale_in), 
	rmax(histSize-2), rmin(1), 
	gmax(histSize-2), gmin(1), 
	bmax(histSize-2), bmin(1), 
	histogram_array(histSize, histSize, histSize),
	processing_array(histSize, histSize, histSize, 0)
{
		mode[0] = 0;
		mode[1] = 0;
		mode[2] = 0;
		a = histogram_array.a;
		sa = processing_array.a;
}

bool Histogram::check_bounds(int d1, int d2, int d3, bool suppress_warning) const
{
	if ( d1 < 0 || d1 >= histSize
		  || d2 < 0 || d2 >= histSize
		  || d3 < 0 || d3 >= histSize )
	{ 
		if(!suppress_warning)
		{ 
			std::cout << "The histogram coordinates (" << d1 << ", " << d2 << ", " << d3 << ") are invalid!!" << std::endl; 
		}
		return false;
	}	
	else
	{ 
		return true; 
	}
}

//Get the count at d1,d2,d3
long int Histogram::v(int d1, int d2, int d3, bool suppress_warning) const
{ 
	if(check_bounds(d1,d2,d3,suppress_warning))
	{ 
		return a[d1][d2][d3]; 
	}
	return 0;
}

void Histogram::set(int d1, int d2, int d3, long int val)
{
	if(check_bounds(d1,d2,d3))
	{
		a[d1][d2][d3] = val; 
	}	
}

//Increment the count at (d1,d2,d3)
void Histogram::inc_element(int d1, int d2, int d3)
{ 
	if(check_bounds(d1,d2,d3))
	{
		if ( a[d1][d2][d3] == 65535 )
		{ 
			std::cerr<<"Histogram increment overflow!!!"; 
		}
		else
		{ 
			a[d1][d2][d3] += inc_scale; 
		}
	}
}

void Histogram::projection_extrema(int dir, long int &max, long int &min)
{
	max = 0;
	min = LONG_MAX;
	int da,db;
	FOR_AXIS(da)
	{ 
		FOR_AXIS(db)
		{
			long int cnt = proj_at(dir,da,db);
			cnt > max ? max = cnt : max = max;
			cnt < min ? min = cnt : min = min;
		}
	}
}

long int Histogram::proj_at(int dir, int da, int db)
{
	long int total = 0;
	int d;
	switch(dir)
	{
	case 1:
		FOR_AXIS(d)
		{ 
			total += a[d][da][db];
		}
		break;
	case 2:
		FOR_AXIS(d)
		{ 
			total += a[da][d][db];
		}
		break;
	case 3:
		FOR_AXIS(d)
		{ 
			total += a[da][db][d];
		}
		break;
	}
	return(total);
}

void Histogram::smooth()
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
				 { 
					 max_freq = a[i][j][k];
				   //mode = _RGB(i, j, k);
					 mode[0] = i;
					 mode[1] = j;
					 mode[2] = k;
				 }
			 }
		 }
	 }
}
 
void Histogram::find_bounding_box()
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
	for( i = 1; i < histSize-1; i++ )
	 { for ( j = 1; j < histSize-1; j++ )
	    { for ( k = 1; k < histSize-1; k++ )
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
	for( i = histSize-2; i > 0; i-- )
	 { for ( j = 1; j < histSize-1; j++ )
	    { for ( k = 1; k < histSize-1; k++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits1 ) goto rmaxdone;
				 }
			 }
		 }
	 }
	rmaxdone:
	rmax = min( i+1, histSize-2 );

	// Find bmin
	th = 0;
	for( k = 1; k < histSize-1; k++ )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( j = 1; j < histSize-1; j++ )
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
	for( k = histSize-2; k > 0; k-- )
	 { for ( i = rmin; i <= rmax; i++ )
	    { for ( j = 1; j < histSize-1; j++ )
          { if(a[i][j][k])
			    { th += a[i][j][k];
				   if( th > minhits2 ) goto bmaxdone;
				 }
			 }
		 }
	 }
	bmaxdone:
	bmax = min( k+1, histSize-2 );

	// Find gmin
	th = 0;
	for( j = 1; j < histSize-1; j++ )
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
	for( j = histSize-2; j > 0; j-- )
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
	gmax = min( j+1, histSize-2 );
 }

void Histogram::dump()
{ 
	int i, j, k;
	char inp[20];
	int planeinc = 1;
	for ( i = 0; i < histSize; i += planeinc )
	{ 
		for ( j = histSize-1; j >= 0; j-- )
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

void Histogram::delete_secondary_blobs()
{ 
	// This will remove all parts of the histogram that are not
	// 6-connected (in 3-d) with the blob containing the mode.

	std::cout << "Removing secondary blobs..." << std::endl;
	
	processing_array.clear_to(0);	
	mark_point_and_nbrs(mode[0], mode[1], mode[2]);

	int i, j, k;
	FOR_AXIS(i)
	{ 
		FOR_AXIS(j)
		{ 
			FOR_AXIS(k)
			{ 
				if (sa[i][j][k] == 0)
				{ 
					a[i][j][k] = 0; 
				}
			 }
		 }
	 }
 }

void Histogram::mark_point_and_nbrs(int d1, int d2, int d3)
{ 
	if ( sa[d1][d2][d3] == 0 && a[d1][d2][d3] > 0 )
    { 
		sa[d1][d2][d3] = 1;
	   
		mark_point_and_nbrs(d1+1, d2, d3);
		mark_point_and_nbrs(d1-1, d2, d3);
		mark_point_and_nbrs(d1, d2+1, d3);
		mark_point_and_nbrs(d1, d2-1, d3);
		mark_point_and_nbrs(d1, d2, d3+1);
		mark_point_and_nbrs(d1, d2, d3-1);
	}
}
	
//*****************************************************************
// ARRAY 3D

Array3D::Array3D(int d1_in, int d2_in, int d3_in, int init, long int init_value)
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
	
Array3D::~Array3D()
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


void Array3D::clear_to (long int val)
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

} //end namespace dh
