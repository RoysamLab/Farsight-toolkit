 #include "dhHistogram.h"

namespace dh
{

Histogram::Histogram(int inc_scale_in) 
: inc_scale(inc_scale_in) 
{
	mode[0] = 0;
	mode[1] = 0;
	mode[2] = 0;
	max_freq = 0;

	//Used in find bounding box and smoothing:
	d1min = 1;
	d2min = 1;
	d3min = 1;
	d1max = histSize-2;
	d2max = histSize-2;
	d3max = histSize-2;

	histImage = HistImageType::New();
	HistImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	histImage->SetOrigin( origin ); 

	HistImageType::IndexType start = { { 0,0,0 } };
	HistImageType::SizeType  size = { { histSize, histSize, histSize } };
	HistImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	histImage->SetRegions( region );           
	histImage->Allocate();        
	histImage->FillBuffer(0);             
	histImage->Update();
}

void Histogram::save_as(const char * fname)
{
	if(!histImage)
		return;

	typedef itk::Image< unsigned char, 3 > UcharImageType;
	typedef itk::RescaleIntensityImageFilter< HistImageType, UcharImageType > RescaleFlUcType;
	RescaleFlUcType::Pointer rescale = RescaleFlUcType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( histImage );
	rescale->Update();

	typedef itk::ImageFileWriter< UcharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( fname );
	writer->SetInput( rescale->GetOutput() );
	writer->Update();
}

bool Histogram::check_bounds(int d1, int d2, int d3) const
{
	if ( d1 < 0 || d1 >= histSize || d2 < 0 || d2 >= histSize || d3 < 0 || d3 >= histSize )
	{ 
		return false;
	}	
	else
	{ 
		return true; 
	}
}

//Get the count at d1,d2,d3
Histogram::PixelType Histogram::v(int d1, int d2, int d3) const
{ 
	if(check_bounds(d1,d2,d3))
	{ 
		HistImageType::IndexType index = { { d1, d2, d3 } };
		return histImage->GetPixel(index);
		//return a[d1][d2][d3]; 
	}
	return 0;
}

void Histogram::set(int d1, int d2, int d3, PixelType val)
{
	if(check_bounds(d1,d2,d3))
	{
		HistImageType::IndexType index = { { d1, d2, d3 } };
		histImage->SetPixel(index, val);
		//a[d1][d2][d3] = val; 
	}	
}

//Increment the count at (d1,d2,d3)
void Histogram::inc_element(int d1, int d2, int d3)
{ 
	if(check_bounds(d1,d2,d3))
	{
		HistImageType::IndexType index = { { d1, d2, d3 } };
		PixelType val = histImage->GetPixel(index);
		//if ( a[d1][d2][d3] >= LONG_MAX )
		if( val >= LONG_MAX )
		{ 
			std::cerr<<"Histogram increment overflow!!!"; 
		}
		else
		{ 
			histImage->SetPixel(index, val+inc_scale);
			//a[d1][d2][d3] += inc_scale; 
		}
	}
}

void Histogram::projection_extrema(int dir, long int &max, long int &min)
{
	max = 0;
	min = LONG_MAX;
	int a, b;
	FOR_AXIS(a)
	{ 
		FOR_AXIS(b)
		{
			long int cnt = proj_at(dir,a,b);
			cnt > max ? max = cnt : max = max;
			cnt < min ? min = cnt : min = min;
		}
	}
}

long int Histogram::proj_at(int dir, int da, int db)
{
	long int total = 0;

	HistImageType::IndexType start = { { 0, 0, 0 } };
	HistImageType::SizeType  size = { { 1, 1, 1 } };
	if(dir == 1)
	{
		size[0] = histSize;
		start[1] = da;
		start[2] = db;
	}
	else if(dir == 2)
	{
		start[0] = da;
		size[1] = histSize;
		start[2] = db;
	}
	else if(dir == 3)
	{
		start[0] = da;
		start[1] = db;
		size[2] = histSize;
	}
	HistImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	typedef itk::ImageRegionIterator< HistImageType > IteratorType;
	IteratorType iterator( histImage, region );
	for( iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator )
	{
		total += iterator.Get();
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
   
	std::cerr << "  Smoothing...";

	Array3D processing_array(histSize, histSize, histSize, 1, 0);
	long int ***sa = processing_array.a;

	//-------------------------------------------------------------------
	for ( int i = d1min; i <= d1max; i++ )
	{ 
		//std::cout << "." << flush;
		for ( int j = d2min; j <= d2max; j++ )
		{ 
			for ( int k = d3min; k < d3max; k++ )
			{ 
				//*************
				long int sum = 0;
				for ( int x = 0; x < 3; x++ )
				{ 
					int ixp = i + x - 1;
					for ( int y = 0; y < 3; y++ )
					{ 
						int jyp = j + y - 1;
						for ( int z = 0; z < 3; z++ )
						{
							int kzp = k + z - 1;
							HistImageType::IndexType index = { { ixp, jyp, kzp } };
							long int val = histImage->GetPixel(index);
							sum += sk[x][y][z] * val;
						}
						//sum +=
						//     sk[x][y][0] * a[ixp][jyp][k-1]
						//	+ sk[x][y][1] * a[ixp][jyp][k]
						//	+ sk[x][y][2] * a[ixp][jyp][k+1];
					}
				}
				sa[i][j][k] = sum / sk_sum;
				//***************
			}
		}
	}
	
	int i, j, k;
	FOR_AXIS(i)
	{ 
		FOR_AXIS(j)
		{ 
			FOR_AXIS(k)
			{ 
				HistImageType::IndexType index = { { i, j, k } };
				histImage->SetPixel(index, sa[i][j][k]);
				//a[i][j][k] = sa[i][j][k];
				if ( sa[i][j][k] > max_freq )
				{ 
					max_freq = sa[i][j][k];
				   //mode = _RGB(i, j, k);
					mode[0] = i;
					mode[1] = j;
					mode[2] = k;
				}
			}
		}
	}

	std::cerr << "...Done" << std::endl;

}
 
void Histogram::find_bounding_box()
{ 
	std::cout << "  Finding bounding box...";
	
	const int minhits1 = 24;    // MAGIC NUMBER (d1)
	const int minhits2 = 24;    // MAGIC NUMBER (d3)
	const int minhits3 = 16;    // MAGIC NUMBER (d2)

	//-------------------------------------------------------------------
	// Efficiency trick: find 3-D bounding box for i, j, k such that
	// there are at least minhits hits in each dimention
	// Do red, then blue, then green.
	
	//Find d1min
	int th = 0;
	int i, j, k;
	for( i = 1; i < histSize-1; i++ )
	{ for ( j = 1; j < histSize-1; j++ )
	  { for ( k = 1; k < histSize-1; k++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits1 ) 
					goto d1mindone;
		    }
	    }
	  }
	}
	d1mindone:
	d1min = max( i-1, 1 );

	// Find d1max
	th = 0;
	for( i = histSize-2; i > 0; i-- )
	{ for ( j = 1; j < histSize-1; j++ )
	  { for ( k = 1; k < histSize-1; k++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits1 ) 
					goto d1maxdone;
			}
	    }
	  }
	}
	d1maxdone:
	d1max = min( i+1, histSize-2 );

	// Find d3min
	th = 0;
	for( k = 1; k < histSize-1; k++ )
	{ for ( i = d1min; i <= d1max; i++ )
	  { for ( j = 1; j < histSize-1; j++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits2 ) 
					goto d3mindone;
			}
	    }
	  }
	}
	d3mindone:
	d3min = max( k-1, 1 );

	// Find d3max
	th = 0;
	for( k = histSize-2; k > 0; k-- )
	{ for ( i = d1min; i <= d1max; i++ )
	  { for ( j = 1; j < histSize-1; j++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits2 ) 
					goto d3maxdone;
			}
	    }
	  }
	}
	d3maxdone:
	d3max = min( k+1, histSize-2 );

	// Find d2min
	th = 0;
	for( j = 1; j < histSize-1; j++ )
	{ for ( i = d1min; i <= d1max; i++ )
	  { for ( k = d3min; k <= d3max; k++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits3 ) 
					goto d2mindone;
			}
		}
      }
	}
	d2mindone:
	d2min = max( j-1, 1 );

	// Find d2max
	th = 0;
	for( j = histSize-2; j > 0; j-- )
	{ for ( i = d1min; i <= d1max; i++ )
	  { for ( k = d3min; k <= d3max; k++ )
        { 
			HistImageType::IndexType index = { { i, j, k } };
			long int val = histImage->GetPixel(index);
			if(val)
			{ 
				th += val;
				if( th > minhits3 ) 
					goto d2maxdone;
				 }
			 }
		 }
	 }
	d2maxdone:
	d2max = min( j+1, histSize-2 );

	std::cout << "...Done" << std::endl;
}

void Histogram::delete_secondary_blobs()
{ 
	// This will remove all parts of the histogram that are not
	// 6-connected (in 3-d) with the blob containing the mode.

	std::cout << "  Removing secondary blobs...";
	
	//Find connected component image:
	typedef itk::ConnectedComponentImageFilter< HistImageType, HistImageType > CCFilterType;
	CCFilterType::Pointer ccfilter = CCFilterType::New();
	ccfilter->SetInput(histImage);
	try
	{
		ccfilter->Update();
	}
	catch( itk::ExceptionObject & excep )
    {
		std::cerr << "    Connected Components: exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return;
    }

	HistImageType::Pointer ccImage = ccfilter->GetOutput();

	//Get ID of the object containing the mode:
	HistImageType::IndexType index = { { mode[0], mode[1], mode[2] } };
	unsigned char id = ccImage->GetPixel(index);

	//Set the histogram count to 0 for all colors not in the main object:
	typedef itk::ImageRegionIterator< HistImageType > IteratorType;
	IteratorType histIterator( histImage, histImage->GetRequestedRegion() );
	IteratorType ccIterator( ccImage, ccImage->GetRequestedRegion() );
	for( histIterator.GoToBegin(), ccIterator.GoToBegin(); !histIterator.IsAtEnd(); ++histIterator, ++ccIterator )
	{
		if(ccIterator.Get() != id)
			histIterator.Set(0);
	}

	std::cout << "...Done" << std::endl;
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
