#include "ColorSegmentation.h"

//Constructor
ColorSegmentation::ColorSegmentation(RGBImageType::Pointer input)
{
	rgb_input = input;

	rli_image = NULL;

	red_weights = NULL;
	blue_weights = NULL;

	IGNORE_BACKGROUND = false;
	LIGHT_BACKGROUND = false; // normally background is black
	TESTING = false;
	GEN_PROJ = false;
}

void ColorSegmentation::TransformToRLI()
{
	if(!rgb_input)
		return;

	int size1 = rgb_input->GetLargestPossibleRegion().GetSize()[0];
	int size2 = rgb_input->GetLargestPossibleRegion().GetSize()[1];
	int size3 = rgb_input->GetLargestPossibleRegion().GetSize()[2];

	rli_image = RLIImageType::New();

	RLIImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	rli_image->SetOrigin( origin );

	RLIImageType::IndexType start = { 0,0,0 };
	RLIImageType::SizeType  size = { size1, size2, size3 };
	RLIImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	rli_image->SetRegions( region ); 
	rli_image->Allocate();                       

	typedef itk::ImageRegionIterator< RLIImageType > IteratorType;
	IteratorType iteratorRLI ( rli_image, rli_image->GetRequestedRegion() );

	typedef itk::ImageRegionConstIterator< RGBImageType > RGBIteratorType;
	RGBIteratorType iteratorRGB( rgb_input, rgb_input->GetRequestedRegion() );

	for( iteratorRLI.GoToBegin(), iteratorRGB.GoToBegin();
	      !iteratorRLI.IsAtEnd() && !iteratorRGB.IsAtEnd();
	      ++iteratorRLI, ++iteratorRGB
		  )
	{ 
		RGBPixelType p_rgb = iteratorRGB.Get();

		dh::_RGB rgb(p_rgb[0], p_rgb[1], p_rgb[2]);
		dh::RLI rli = (dh::RLI)rgb;

		RLIPixelType p_rli;
		p_rli[0] = rli.R;
		p_rli[1] = rli.L;
		p_rli[2] = rli.I;
		iteratorRLI.Set( p_rli );
	}

	if(TESTING)
	{
		typedef itk::ImageFileWriter< RLIImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();

		writer->SetFileName( "rli.tif" );
		writer->SetInput( rli_image );
		writer->Update();
	}
}

void ColorSegmentation::SetArchetypalColors(dh::RLI r, dh::RLI b, dh::RLI w)
{
	archTypRED = r;
	archTypBLUE = b;
	archTypBACK = w;
}

void ColorSegmentation::SetArchetypalColors(dh::_RGB r, dh::_RGB b, dh::_RGB w)
{
	archTypRED = (dh::RLI)r;
	archTypBLUE = (dh::RLI)b;
	archTypBACK = (dh::RLI)w;
}

//Run Pixel classification
void ColorSegmentation::FindArchetypalColors()
{
	if(!rli_image)
		return;

	// ========== Create 3D Histogram =========
	dh::Histogram * hist = new dh::Histogram(2);

	typedef itk::ImageRegionConstIterator< RLIImageType > IteratorType;
	IteratorType it( rli_image, rli_image->GetRequestedRegion() );
	for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		RLIPixelType p_rli = it.Get();
		dh::RLI P = dh::RLI( p_rli[0], p_rli[1], p_rli[2] );	//Creating a 3D "pixel"
		hist->inc( P );
	}

	hist->find_bounding_box();
	hist->smooth();
	hist->delete_secondary_blobs();

	std::cout << "Mode " << hist->modeAsRLI() * 2 << " has max frequency " << hist->max_freq << std::endl;

	//===================================================================
	// =========== Find Color Archetypes ================================

	dh::SeedGrid seed_grid( hist, LIGHT_BACKGROUND );

	dh::_RGB a1, a2, bkgrnd; //REMEMBER That a1 and a2 are actually RLI values

	seed_grid.find_seeds( a1, a2, bkgrnd );

	std::cout << " Seed1 = " << a1.mapRGBtoRLI() * 2 << std::endl;
	std::cout << " Seed2 = " << a2.mapRGBtoRLI() * 2 << std::endl;
	std::cout << " Bgrnd = " << bkgrnd.mapRGBtoRLI() * 2 << std::endl;

	seed_grid.find_most_distinct_colors( a1, a2, bkgrnd );

	if( a1.B > a2.B )
	{
		archTypRED = a1.mapRGBtoRLI()*2; 
		archTypBLUE = a2.mapRGBtoRLI()*2;
	}
	else
	{
		archTypRED = a2.mapRGBtoRLI()*2; 
		archTypBLUE = a1.mapRGBtoRLI()*2;
	}
 
	archTypBACK = bkgrnd.mapRGBtoRLI()*2;

	std::cout << "FOUND ARCHETYPES: " << std::endl;
	std::cout << " RED =        " << archTypRED << std::endl;
	std::cout << " BLUE =       " << archTypBLUE << std::endl;
	std::cout << " BACKGROUND = " << archTypBACK << std::endl;

	delete hist;
}

void ColorSegmentation::ComputeClassWeights()
{ 
	std::cout << "USING ARCHETYPES: " << std::endl;
	std::cout << " RED =        " << archTypRED << std::endl;
	std::cout << " BLUE =       " << archTypBLUE << std::endl;
	std::cout << " BACKGROUND = " << archTypBACK << std::endl;

	//Preliminary computations:
	float R_B_axis_len = dh::Classifier::euclidean_dist (archTypRED, archTypBLUE );
	float Y_B_axis_len = dh::Classifier::euclidean_dist (archTypBACK, archTypBLUE );
	float Y_R_axis_len = dh::Classifier::euclidean_dist (archTypBACK, archTypRED );
	
	float slide_wt = 0.7;//0.7; //Set to magic number
	//Set slide_wt to deal with differing spread
	//slide_wt = Y_R_axis_len / ( Y_B_axis_len + Y_R_axis_len );

	//----------- Find decision planes ---------------
	// 1. Find D (split point)
	dh::XYZ r = archTypRED;
	dh::XYZ b = archTypBLUE;
	dh::XYZ w = archTypBACK;
	const double spw = 0.50;
	dh::XYZ d =  (b *(1-spw)+ r*spw);
	dh::XYZ split_point = d;
	
	// 2. Find p (decision plan for 2 colors - represented by normal vector)
	dh::XYZ wb = b - w;
	dh::XYZ wr = r - w;
	dh::XYZ wd = d - w;
	dh::XYZ wbwr = dh::cross( wb, wr );
	dh::XYZ p = dh::cross ( wbwr, wd );
	dh::XYZ decision_plane = p / magnitude ( p );	//decision plane normal vector

	p = dh::cross ( wbwr, wr );
	dh::XYZ red_plane = p / magnitude ( p );	//red plane normal vector
	double red_sp_dist = dh::dot ( (d - r), red_plane ); //scalar projection of RD onto red plane unit normal vector
	//This is the distance between d and the red plane (negative = inside triangle)

	p = dh::cross ( wbwr, wb );
	dh::XYZ blue_plane = p / magnitude ( p );	//blue plane normal vector
	double blue_sp_dist = dh::dot ( (d - b), blue_plane );//scalar projection of BD onto blue plane unit normal vector
	//Distance between d and the blue plane (positive = inside triangle)

	if(TESTING)
	{ 
		std::cerr << "Distance between red & blue: "
				  << R_B_axis_len << std::endl
			      << "Distance between bkgrnd & blue: "
				  << Y_B_axis_len << std::endl
			      << "Distance between bkgrnd & red: "
				  << Y_R_axis_len << std::endl;

		std::cerr //<< "split_point: " << split_point << std::endl
				//<< "decision_plane: " << decision_plane << std::endl
				<< "red_sp_dist: " << red_sp_dist << std::endl
				<< "blue_sp_dist: " << blue_sp_dist << std::endl;
  
		std::cerr << "Using slide_wt: " << slide_wt << std::endl;
	}

	//------------------------------------------------
	// FIND WEIGHTS FOR EACH PIXEL:
	//**************************************
	int size1 = rgb_input->GetLargestPossibleRegion().GetSize()[0];
	int size2 = rgb_input->GetLargestPossibleRegion().GetSize()[1];
	int size3 = rgb_input->GetLargestPossibleRegion().GetSize()[2];

	red_weights = UcharImageType::New();
	blue_weights = UcharImageType::New();

	UcharImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	red_weights->SetOrigin( origin ); 
	blue_weights->SetOrigin( origin );

	UcharImageType::IndexType start;
	start[0] = 0; // first index on X
	start[1] = 0; // first index on Y
	start[2] = 0; // first index on Z
	UcharImageType::SizeType  size;
	size[0] = size1; // size along X
	size[1] = size2; // size along Y
	size[2] = size3; // size along Z
	UcharImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	red_weights->SetRegions( region ); 
	blue_weights->SetRegions( region );
	red_weights->Allocate();            
	blue_weights->Allocate();
	red_weights->FillBuffer(0);         
	blue_weights->FillBuffer(0);
	red_weights->Update();              
	blue_weights->Update();

	typedef itk::ImageRegionIterator< UcharImageType > UcharIteratorType;
	UcharIteratorType iterator1 ( red_weights, red_weights->GetRequestedRegion() );
	UcharIteratorType iterator2 ( blue_weights, blue_weights->GetRequestedRegion() );
	
	typedef itk::ImageRegionConstIterator< RLIImageType > rliIteratorTypeConst;
	rliIteratorTypeConst rli_it( rli_image, rli_image->GetRequestedRegion() );
	
	for ( rli_it.GoToBegin(), iterator1.GoToBegin(), iterator2.GoToBegin();
	      !rli_it.IsAtEnd() && !iterator1.IsAtEnd() && !iterator2.IsAtEnd();
	      ++rli_it, ++iterator1, ++iterator2 )
	{
		RLIPixelType p_rli = rli_it.Get();
		dh::RLI pixel = dh::RLI( p_rli[0], p_rli[1], p_rli[2] );	//Creating a 3D "pixel"

		float red_dist = 2.0 * ( 1.0 - slide_wt ) *
			dh::Classifier::euclidean_dist( pixel, archTypRED );	//Distance to red archtype
		
		float blue_dist = 2.0 * slide_wt *
			dh::Classifier::euclidean_dist( pixel, archTypBLUE );//Distance to blue archtype
		
		float bkgrnd_dist = dh::Classifier::euclidean_dist( pixel, archTypBACK );//Distance to background archtype

		Pixel_Class pixel_class;

		//certainty is a measure of how close to the decision plane the point is
		//does not take into account how close it may be to the background color
		float certainty;  	// 1 = certain; 0 = unknown

		//Distance from pixel to decision plane (pos=red, neg=blue):
		double s_plane_dist = dh::dot ( ((dh::XYZ)pixel - split_point), decision_plane );

		if ( s_plane_dist >= 0 ) 
		{ 
			// Pixel is "red"
			pixel_class = RED_CELL;
			//Distance from pixel to red plane (pos=outside, neg=inside)
			double r_plane_dist = dh::dot ( ((dh::XYZ)pixel - (dh::XYZ)archTypRED), red_plane );
			// This could be done with signum functions to make it faster...
			if ( r_plane_dist * red_sp_dist < 0 ) //If signs are different the point is outside the triangle
			{
				certainty = 1.0;
			}
			else
			{
				certainty = s_plane_dist / ( s_plane_dist + fabs(r_plane_dist) );
				//certainty = 1.0;
			}
		}
		else
		{ 
			// Pixel is "blue"
			pixel_class = BLUE_CELL;
			double b_plane_dist = dh::dot ( ((dh::XYZ)pixel - (dh::XYZ)archTypBLUE), blue_plane );
			// This could be done with signum functions to make it faster...
			if ( b_plane_dist * blue_sp_dist < 0 )
			{
				certainty = 1.0;
			}
			else
			{
				s_plane_dist = - s_plane_dist;
				certainty = s_plane_dist / ( s_plane_dist + fabs(b_plane_dist) );
				//certainty = 1.0;
			}
		}

		float bkgrnd_wt = 1 - certainty;

		
		if(IGNORE_BACKGROUND)
		{
			// If pixel is near background color, ignore it
			if(pixel_class == RED_CELL)
			{ 
				pixel_class = (red_dist * (1 - bkgrnd_wt))
			                < (bkgrnd_dist * bkgrnd_wt)
							? RED_CELL : BKGD_FIELD; 
			}
			else
			{ 
				pixel_class = (blue_dist * (1 - bkgrnd_wt))
			                < (bkgrnd_dist * bkgrnd_wt)
			   				? BLUE_CELL : BKGD_FIELD; 
			}

		}
		else
		{
			// If pixel is near background color, discount it by cubing certainty factor (which is < 1.0)
			if( ( pixel_class == RED_CELL && (red_dist * (1 - bkgrnd_wt) < bkgrnd_dist * bkgrnd_wt))
				||(pixel_class == BLUE_CELL && (blue_dist * (1 - bkgrnd_wt) < bkgrnd_dist * bkgrnd_wt))
				) 
			{
				certainty = pow(certainty, 3);
			}
		}
		

		/*
		if( pixel_class != RED_CELL && pixel_class != BLUE_CELL )
		{
			iterator1.Set( blue_dist/(blue_dist+red_dist)*certainty );  //Function1
			iterator2.Set( red_dist/(blue_dist+red_dist) *certainty );  //Function2
		}
		else
		{
			iterator1.Set( (blue_dist+bkgrnd_dist)/(blue_dist+red_dist+bkgrnd_dist)*certainty );  //Function3
			iterator2.Set( (red_dist+bkgrnd_dist) /(blue_dist+red_dist+bkgrnd_dist)*certainty );  //Function4
		}
		*/
		
		switch(pixel_class)
		{
		case BKGD_FIELD:
			break;
		case RED_CELL:
			iterator1.Set((unsigned char)((1-certainty)*255));
			//iterator1.Set(255);
			break;
		case BLUE_CELL:
			iterator2.Set((unsigned char)((1-certainty)*255));
			//iterator2.Set(255);
			break;
		default:
			std::cerr << "INVALID PIXEL CLASS" << std::endl;
		}
	}

	if(TESTING)
	{
		typedef itk::ImageFileWriter< UcharImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( "red_wts.tif" );
		writer->SetInput( red_weights );
		writer->Update();
		WriterType::Pointer writer1 = WriterType::New();
		writer1->SetFileName( "blue_wts.tif" );
		writer1->SetInput( blue_weights );
		writer1->Update();
	}
	
}

/*
dh::Slice3D_Space * ColorSegmentation::GenerateProjection(dh::Histogram * hist, dh::SeedGrid *seed_grid)
{
	//================== Create Histogram Projections ==================
	std::cout << "Creating histogram images..." << std::endl;
	const int incr = 6;  // Was 3; changed 11/4/97

	dh::Col_Slice* RB_Hist = new dh::RI_Slice();
	dh::Col_Slice* RG_Hist = new dh::RL_Slice();
	dh::Col_Slice* GB_Hist = new dh::LI_Slice();
	
	dh::Slice3D_Space * hist_imgs = new dh::Slice3D_Space( RB_Hist, RG_Hist, GB_Hist );

	//--- Show bounding box ---
	dh::_RGB bbox_min(2 * hist->rmin, 2 * hist->gmin, 2 * hist->bmin);
	dh::_RGB bbox_max(2 * hist->rmax, 2 * hist->gmax, 2 * hist->bmax);

	hist_imgs->plot_cross( bbox_min, dh::RGB_GRAY(30), 2 );
	hist_imgs->plot_cross( bbox_max, dh::RGB_GRAY(50), 2 );
	
	//--- Increment histogram projections ---
	typedef itk::ImageRegionConstIterator< UcharImageType > UcharIteratorTypeConst;
	UcharIteratorTypeConst pix_buf1( red_image, red_image->GetRequestedRegion() );
	UcharIteratorTypeConst pix_buf2( lime_image, lime_image->GetRequestedRegion() );
	UcharIteratorTypeConst pix_buf3( intensity_image, intensity_image->GetRequestedRegion() );
	
	for ( pix_buf1.GoToBegin(),	pix_buf2.GoToBegin(), pix_buf3.GoToBegin();
	      !pix_buf1.IsAtEnd() && !pix_buf2.IsAtEnd() && !pix_buf3.IsAtEnd();
	      ++pix_buf1, ++pix_buf2, ++pix_buf3 )
	{
		dh::_RGB pixel = dh::_RGB( pix_buf1.Get(), pix_buf2.Get(), pix_buf3.Get() );
		hist_imgs->inc_freq( pixel, incr );
	}

	hist_imgs->print_overflow();

    //--- Superimpose smoothed histogram ---
	std::cout << "Superimposing Smoothed Histogram..." << std::endl;
	
	int i,j;
	FOR_AXIS(i)
	{ 
		FOR_AXIS(j)
		{ 
			if ( hist->RB_proj_at(i, j) > 0 )
			{ 
				super_point( RB_Hist, dh::_RGB(2*i, 128, 2*j) );
	 			super_point( RB_Hist, dh::_RGB(2*i+1, 128, 2*j) );
	 	 		super_point( RB_Hist, dh::_RGB(2*i, 128, 2*j+1) );
	 			super_point( RB_Hist, dh::_RGB(2*i+1, 128, 2*j+1) );
	 		}
	 	}
	}
	
	FOR_AXIS(i)
	{ 
		FOR_AXIS(j)
		{ 
			if ( hist->RG_proj_at(i, j) > 0 )
	 		{ 
				super_point( RG_Hist, dh::_RGB(2*i, 2*j, 128) );
	 			super_point( RG_Hist, dh::_RGB(2*i+1, 2*j, 128) );
	 			super_point( RG_Hist, dh::_RGB(2*i, 2*j+1, 128) );
	 			super_point( RG_Hist, dh::_RGB(2*i+1, 2*j+1, 128) );
	 		}
	 	}
	}

	FOR_AXIS(i)
	{ 
		FOR_AXIS(j)
		{ 
			if ( hist->GB_proj_at(i, j) > 0 )
	 		{ 
				super_point( GB_Hist, dh::_RGB(128, 2*i, 2*j) );
	 			super_point( GB_Hist, dh::_RGB(128, 2*i+1, 2*j) );
	 			super_point( GB_Hist, dh::_RGB(128, 2*i, 2*j+1) );
	 			super_point( GB_Hist, dh::_RGB(128, 2*i+1, 2*j+1) );
	 		}
	 	}
	}

	//Plot Seeds
	std::list<dh::_RGB> sds = seed_grid->s;
	std::list<dh::_RGB>::iterator it;
	for( it=sds.begin(); it != sds.end(); it++ )
	{
		dh::_RGB pcl = *it;
		dh::_RGB point = ((dh::XYZ)hist->modeAsRGB() 
			+ ((dh::XYZ)pcl-(dh::XYZ)seed_grid->center)
			* seed_grid->sample_dist)*2;
		hist_imgs->plot_point( point, dh::_RGB(252, 200, 252) );
	}
	
	return hist_imgs;
	

}
*/
/*
void ColorSegmentation::super_point( dh::Col_Slice* slice, dh::_RGB clr )
{ 
	if ( slice->point_color( clr ) == dh::RGB_BLACK )
	{ 
		slice->plot_point( clr, dh::RGB_GRAY(32) ); 
	}
	else
	{ 
		slice->cut_intensity( 2, clr ); 
	}
}
*/
/*
UcharImageType::Pointer ColorSegmentation::ComputeBinary(int num_bins, int num_in_fg, bool fgrnd_dark)
{
	
	if(!intensity_image)
		return NULL;

	//Create histogram:
	typedef itk::Statistics::ScalarImageToHistogramGenerator< UcharImageType > HistogramGeneratorType;
	typedef HistogramGeneratorType::HistogramType HistogramType;
	HistogramGeneratorType::Pointer histoGenerator = HistogramGeneratorType::New();
	histoGenerator->SetNumberOfBins( 128 );
	histoGenerator->SetInput( intensity_image );
	histoGenerator->Compute();

	//Estimate thresholds:
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;
	CalculatorType::Pointer calculator = CalculatorType::New();
	calculator->SetNumberOfThresholds( num_bins );
	calculator->SetInputHistogram( histoGenerator->GetOutput() );
	calculator->Update();

	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();

	//Estimate Binarization thresholds:
	int lowerThreshold,upperThreshold;
	if( fgrnd_dark )	//Do I want to make the foregound of the binary dark:
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.end();
 		for( int i=0; i<(num_bins-num_in_fg+1); ++i ) 
			--itNum;
		upperThreshold = static_cast<UcharPixelType>(*itNum);
		lowerThreshold = itk::NumericTraits<UcharPixelType>::min();
	} 
	else 
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
		for( int i=0; i<(num_bins-num_in_fg); ++i ) 
			++itNum;
		lowerThreshold = static_cast<UcharPixelType>(*itNum);
		upperThreshold = itk::NumericTraits<UcharPixelType>::max();
	}

	typedef itk::BinaryThresholdImageFilter< UcharImageType, UcharImageType >  ThreshFilterType;
	ThreshFilterType::Pointer threshfilter = ThreshFilterType::New();
	threshfilter->SetOutsideValue( 0 );
	threshfilter->SetInsideValue( (int)itk::NumericTraits<UcharPixelType>::max() );
	threshfilter->SetInput( intensity_image );
	threshfilter->SetLowerThreshold( lowerThreshold );
	threshfilter->SetUpperThreshold( upperThreshold );
	threshfilter->Update();

	if(TESTING)
	{
		typedef itk::ImageFileWriter< UcharImageType > WriterType;
		WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( "bin.tif" );
		writer->SetInput( threshfilter->GetOutput() );
		writer->Update();
	}

	return threshfilter->GetOutput();
}
*/