#include "ColorSegmentation.h"

//Constructor
ColorSegmentation::ColorSegmentation(RGBImageType::Pointer input, int fore_ground_dark, int number_bins, int number_in_fg)
{
	rgb_input = input;
	foreground_dark = fore_ground_dark;
	number_of_bins = number_bins;
	number_in_foreground = number_in_fg;
	bin_done = 0;

	RLI_MODE = true;
	GEN_PROJ = false;

	hist = NULL;
}

//Destructor
ColorSegmentation::~ColorSegmentation()
{
}

UcharImageType::Pointer ColorSegmentation::get_binary()
{
	if( bin_done ) 
		return intial_binarization;
	else
	{ 
		std::cerr<<"Run Binarization First\n";
		return NULL;
	}
}

//Run Initial Binarization
void ColorSegmentation::RunInitialBinarization()
{
	typedef itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorType;
	typedef itk::ImageRegionConstIterator< RGBImageType > RGBIteratorType;

	typedef itk::Statistics::ScalarImageToHistogramGenerator< FloatImageType > ScalarImageToHistogramGeneratorType;
	typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;
	typedef itk::BinaryThresholdImageFilter< FloatImageType, UcharImageType >  ThreshFilterType;

	int size1 = rgb_input->GetLargestPossibleRegion().GetSize()[0];
	int size2 = rgb_input->GetLargestPossibleRegion().GetSize()[1];
	int size3 = rgb_input->GetLargestPossibleRegion().GetSize()[2];

	FloatImageType::Pointer intensity_image1 = FloatImageType::New();
	FloatImageType::Pointer red_image1 = FloatImageType::New();
	FloatImageType::Pointer lime_image1 = FloatImageType::New();

	FloatImageType::PointType origin1, origin2, origin3;
	origin1[0] = 0; origin2[0] = 0; origin3[0] = 0;
	origin1[1] = 0; origin2[1] = 0; origin3[1] = 0;
	origin1[2] = 0; origin2[2] = 0; origin3[2] = 0;
	red_image1->SetOrigin( origin1 ); 
	lime_image1->SetOrigin( origin2 ); 
	intensity_image1->SetOrigin( origin3 );

	FloatImageType::IndexType start1, start2, start3;
	start1[0] = 0; start2[0] = 0; start3[0] = 0;  // first index on X
	start1[1] = 0; start2[1] = 0; start3[1] = 0;  // first index on Y
	start1[2] = 0; start2[2] = 0; start3[2] = 0;  // first index on Z
	FloatImageType::SizeType  size11, size22, size33;
	size11[0] = size1; size22[0] = size1; size33[0] = size1;  // size along X
	size11[1] = size2; size22[1] = size2; size33[1] = size2;  // size along Y
	size11[2] = size3; size22[2] = size3; size33[2] = size3;  // size along Z
	FloatImageType::RegionType region1, region2, region3;
	region1.SetSize( size11 );  region2.SetSize( size22 );  region3.SetSize( size33 );
	region1.SetIndex( start1 ); region2.SetIndex( start2 ); region3.SetIndex( start3 );
	red_image1->SetRegions( region1 ); 
	lime_image1->SetRegions( region2 ); 
	intensity_image1->SetRegions( region3 );
	red_image1->Allocate();            
	lime_image1->Allocate();            
	intensity_image1->Allocate();
	red_image1->FillBuffer(0);         
	lime_image1->FillBuffer(0);         
	intensity_image1->FillBuffer(0);
	red_image1->Update();              
	lime_image1->Update();              
	intensity_image1->Update();

	IteratorType iterator1 ( red_image1,       red_image1      ->GetRequestedRegion() );
	IteratorType iterator2 ( lime_image1,      lime_image1     ->GetRequestedRegion() );
	IteratorType iterator3 ( intensity_image1, intensity_image1->GetRequestedRegion() );

	RGBIteratorType pix_buf1( rgb_input, rgb_input->GetRequestedRegion() );
	iterator1.GoToBegin();
	iterator2.GoToBegin();
	iterator3.GoToBegin();
	for ( pix_buf1.GoToBegin();
	      !pix_buf1.IsAtEnd() && !iterator1.IsAtEnd() && !iterator2.IsAtEnd() && !iterator3.IsAtEnd() ;
	      ++pix_buf1, ++iterator1, ++iterator2, ++iterator3 )
	{ 
		RGBPixelType pix_vals;
		pix_vals = pix_buf1.Get();

		const float red   = (float)pix_vals[0];
		const float green = (float)pix_vals[1];
		const float blue  = (float)pix_vals[2];

		float red1, lime, intensity;

		float total = red + green + blue;
		//intensity = total / (255*3);
		//TRY THIS
		intensity = (0.3*red+0.59*green+0.11*blue)/255;	 

		float r = (float)red   / total;
		float g = (float)green / total;
		float b = (float)blue  / total;

		float s = 1 - 3 * ((r<g) ? ((r<b) ? r : b ) : ((g<b) ? g : b ));

		float h = 0.5 * ( (r-g) + (r-b) ) / (sqrt((r-g)*(r-g)+(r-b)*(g-b)) );
		if(b > g)
			h = 2*M_PI - h;
		h =  h+1;			//************************************************needed???
		float cr =  cos( h );
		float cl =  sin( h );

		lime = s * cl;
		red1 = s * cr;

		iterator1.Set( red1 * 127 );
		iterator2.Set( lime * 127 );
		iterator3.Set( intensity * 127 );
	}
	ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
	ThreshFilterType::Pointer threshfilter = ThreshFilterType::New();
	CalculatorType::Pointer calculator = CalculatorType::New();
	scalarImageToHistogramGenerator->SetNumberOfBins( 128 );
	calculator->SetNumberOfThresholds( number_of_bins );
	threshfilter->SetOutsideValue( 0 );
	threshfilter->SetInsideValue( (int)itk::NumericTraits<UcharPixelType>::max() );

	scalarImageToHistogramGenerator->SetInput( intensity_image1 );
	scalarImageToHistogramGenerator->Compute();
	calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
	threshfilter->SetInput( intensity_image1 );
	calculator->Update();
	const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();

	int lowerThreshold,upperThreshold;
	if( number_in_foreground > number_of_bins ){
		std::cerr<<"Number of thresholds included in the foreground is greater than the number of thresholds computed\n";
		return;
	}

	//testing
	CalculatorType::OutputType::const_iterator itNum1 = thresholdVector.begin();
	for( int i=0; i<number_of_bins; ++i ) {
		std::cout<<(int)static_cast<UcharPixelType>(*itNum1)<<std::endl;
		++itNum1;
	}//

	if( foreground_dark )
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.end();
 		for( int i=0; i<(number_of_bins-number_in_foreground+1); ++i ) --itNum;
		upperThreshold = static_cast<UcharPixelType>(*itNum);
		lowerThreshold = itk::NumericTraits<UcharPixelType>::min();
	} 
	else 
	{
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
		for( int i=0; i<(number_of_bins-number_in_foreground+1); ++i ) ++itNum;
		lowerThreshold = static_cast<UcharPixelType>(*itNum);
		upperThreshold = itk::NumericTraits<UcharPixelType>::max();
	}

	threshfilter->SetLowerThreshold( lowerThreshold );
	threshfilter->SetUpperThreshold( upperThreshold );
	threshfilter->Update();

	intial_binarization = UcharImageType::New();
	intial_binarization = threshfilter->GetOutput();

	typedef itk::RescaleIntensityImageFilter< FloatImageType, UcharImageType > RescaleFlUcType;
	RescaleFlUcType::Pointer RescaleIntIO1 = RescaleFlUcType::New();
	RescaleIntIO1->SetOutputMaximum( 127 );
	RescaleIntIO1->SetOutputMinimum( 0 );
	RescaleIntIO1->SetInput( intensity_image1 );
	RescaleIntIO1->Update();
	RescaleFlUcType::Pointer RescaleIntIO2 = RescaleFlUcType::New();
	RescaleIntIO2->SetOutputMaximum( 127 );
	RescaleIntIO2->SetOutputMinimum( 0 );
	RescaleIntIO2->SetInput( lime_image1 );
	RescaleIntIO2->Update();
	RescaleFlUcType::Pointer RescaleIntIO3 = RescaleFlUcType::New();
	RescaleIntIO3->SetOutputMaximum( 127 );
	RescaleIntIO3->SetOutputMinimum( 0 );
	RescaleIntIO3->SetInput( red_image1 );
	RescaleIntIO3->Update();
	red_image = RescaleIntIO3->GetOutput();
	lime_image = RescaleIntIO2->GetOutput();
	intensity_image = RescaleIntIO1->GetOutput();

	float max=(float)LLONG_MIN;
	float min=(float)LLONG_MAX;
	iterator2.GoToBegin();
	for ( pix_buf1.GoToBegin();
	      !iterator2.IsAtEnd();
	      ++iterator2 ){
		if( max < iterator2.Get() ) max = iterator2.Get();
		if( min > iterator2.Get() ) min = iterator2.Get();
	}
	std::cout<<"MAX:"<<max<<"\nMIN:"<<min<<std::endl;
	typedef itk::ImageFileWriter< UcharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "intensity.tif" );
	writer->SetInput( intensity_image );
	writer->Update();
	WriterType::Pointer writer1 = WriterType::New();
	writer1->SetFileName( "lime.tif" );
	writer1->SetInput( lime_image );
	writer1->Update();
	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetFileName( "red.tif" );
	writer2->SetInput( red_image );
	writer2->Update();

	bin_done = 1;
}

//Run Pixel classification
void ColorSegmentation::FindArchetypalColors()
{
	dh::RGB_Atype at;

	typedef itk::ImageRegionConstIterator< UcharImageType > UcharIteratorType;

	UcharIteratorType pix_buf3( red_image, red_image->GetRequestedRegion() );
	UcharIteratorType pix_buf2( lime_image, lime_image->GetRequestedRegion() );
	UcharIteratorType pix_buf1( intensity_image, intensity_image->GetRequestedRegion() );

	// ========== Create 3D Histogram =========
	if(hist)
		delete hist;
	hist = new dh::RGBHistogram(2);

	pix_buf1.GoToBegin();	pix_buf2.GoToBegin();	pix_buf3.GoToBegin();
	for ( ;
	      !pix_buf1.IsAtEnd() && !pix_buf2.IsAtEnd() && !pix_buf3.IsAtEnd();
	      ++pix_buf1, ++pix_buf2, ++pix_buf3 )
	{
		dh::_RGB P = dh::_RGB( pix_buf1.Get(), pix_buf2.Get(), pix_buf3.Get() );
		hist->inc( P );
	}

	hist->find_bounding_box();
	hist->smooth();
	hist->delete_secondary_blobs();

	std::cout << "Mode \t" << ((int)(hist->mode.R)) * 2 << ", " << ((int)(hist->mode.G)) * 2<< ", " << ((int)(hist->mode.B)) * 2
		<<" has max frequency " << hist->max_freq<< "." << std::endl;

	//hist->dump();

	//================== Create Histogram Projections ==================
	std::cout << "Creating histogram images..." << std::endl;

	dh::Col_Slice* RB_Hist;
	dh::Col_Slice* RG_Hist;
	dh::Col_Slice* GB_Hist;
	
	if (RLI_MODE)
	{ 
		RB_Hist = new dh::RI_Slice();
		RG_Hist = new dh::RL_Slice();
		GB_Hist = new dh::LI_Slice();
	}
	else
	{ 
		dh::RGB_Slice RB_Hist( 1, dh::_RGB(0, 0, 0), dh::aG, 128, dh::aB, dh::POS, dh::aR, dh::POS );
		dh::RGB_Slice RG_Hist( 1, dh::_RGB(0, 0, 0), dh::aB, 128, dh::aG, dh::NEG, dh::aR, dh::POS );
		dh::RGB_Slice GB_Hist( 1, dh::_RGB(0, 0, 0), dh::aR, 128, dh::aB, dh::POS, dh::aG, dh::NEG ); 
	}	

	dh::Slice3D_Space hist_imgs( RB_Hist, RG_Hist, GB_Hist );

	//===================================================================
	// =========== Find Color Archetypes ================================
	at.bkgrnd = hist->mode;

	const int max_res = 15;   // MAGIC NUMBER !!! 
	double res[max_res+1];

	// ======= Find Rough Estimate (initial seeds) ===========	

	const int resolution = 8;	// MAGIC NUMBER !!! ************************ used to be 4
								// Tune for accuracy (hi) vs. Speed (lo)
								// Must use > 1.

	dh::SeedGrid seed_grid( hist, resolution, (bool)foreground_dark );

	dh::_RGB a1, a2;

	seed_grid.find_seeds( a1, a2, dh::eval_state );

	std::cout << "   Seed1 = " << a1.R * 2 << " "<< a1.G * 2 << " " <<  a1.B * 2 << " " << std::endl;
	std::cout << "   Seed2 = " << a2.R * 2 << " "<< a2.G * 2 << " " <<  a2.B * 2 << " " << std::endl;

	std::cout << "Finding most distinct colors..." << std::endl;
	
	bool a1_moved = false;
	bool a2_moved = false;
	
	do
	{ 
		do
	    { 
			for ( int iter = 0; iter <= max_res; iter++ )
 			{ 
				go_best_dir( a1, a1_moved, a2 );
				go_best_dir( a2, a2_moved, a1 );
				res[iter] = dh::eval_state ( at.bkgrnd, a1, a2 );
		    }
	    } while ( ( a1_moved || a2_moved ) &&
		           fabs( res[max_res] - res[0] ) >
		           fabs( .02 * res[max_res]) );    // MAGIC NUMBER !!!
      
		int big_jump_size = 8;       // MAGIC NUMBER !!!

		go_best_dir( a1, a1_moved, a2, big_jump_size, 1 );
		go_best_dir( a2, a2_moved, a1, big_jump_size, 1 );

	} while ( a1_moved || a2_moved );

	std::cout << "   Color1 = " << a1.R * 2 << " "<< a1.G * 2 << " " <<  a1.B * 2 << " " << std::endl;
	std::cout << "   Color2 = " << a2.R * 2 << " "<< a2.G * 2 << " " <<  a2.B * 2 << " " << std::endl;

	if( a1.B > a2.B )
	{
		arch_typ1= a1; 
		arch_typ2= a2;
	}
	else
	{
		arch_typ1= a2; 
		arch_typ2= a1;
	}
	
	bkgrnd_typ = at.bkgrnd;
}

void ColorSegmentation::ComputeClassWeights()
{ 
	dh::RGB_Atype atype;
	atype.red = arch_typ1;
	atype.blue = arch_typ2;
	atype.bkgrnd = bkgrnd_typ;

	float R_B_axis_len = dh::RGB_Classifier::RGB_euclidean_dist (atype.red, atype.blue );
	float Y_B_axis_len = dh::RGB_Classifier::RGB_euclidean_dist (atype.bkgrnd, atype.blue );

	// Set slide_wt to deal with differing spread
	//	slide_wt = Y_R_axis_len / ( Y_B_axis_len + Y_R_axis_len );

 	// Enable this line and disable previous to
			// totally disable slide_wt

	float slide_wt = .7;
	Pixel_Class pixel_class;

	//----------- Find decision planes ---------------
	dh::XYZ r = atype.red;
	dh::XYZ b = atype.blue;
	dh::XYZ w = atype.bkgrnd;
	const double spw = 0.50;
	dh::XYZ d =  (b *(1-spw)+ r*spw);
	static dh::XYZ split_point = d;
	//std::cout << "SP" << split_point << std::endl;
	dh::XYZ wb = b - w;
	dh::XYZ wr = r - w;
	dh::XYZ wd = d - w;
	//DMP(wb);
	//DMP(wr);
	//DMP(wd);
	dh::XYZ p = cross ( cross ( wb , wr ), wd );
	//DMP(p);
	//cout << "magp = " << magnitude ( p ) << endl;
	static dh::XYZ decision_plane = p / magnitude ( p );
	//DMP( RLI_WV_Classifier_2::decision_plane );
	//cout << RLI_WV_Classifier_2::decision_plane << endl;
	p = dh::cross ( dh::cross ( wb , wr ), wr );
	static dh::XYZ red_plane = p / magnitude ( p );
	static double red_sp_dist = dh::dot ( (d - r), red_plane );

	p = dh::cross ( dh::cross ( wb , wr ), wb );
	static dh::XYZ blue_plane = p / magnitude ( p );
	static double blue_sp_dist = dh::dot ( (d - b), blue_plane );
	
	//DMP( RLI_WV_Classifier_2::decision_plane );
	//DMP( RLI_WV_Classifier_2::red_plane );
	//DMP( RLI_WV_Classifier_2::blue_plane );
	//DMP( RLI_WV_Classifier_2::red_sp_dist );
	//DMP( RLI_WV_Classifier_2::blue_sp_dist );

	//------------------------------------------------

//**************************************
	typedef itk::ImageRegionIteratorWithIndex< FloatImageType > IteratorType;
	typedef itk::ImageRegionConstIterator< UcharImageType > UcharIteratorType;

	int size1 = rgb_input->GetLargestPossibleRegion().GetSize()[0];
	int size2 = rgb_input->GetLargestPossibleRegion().GetSize()[1];
	int size3 = rgb_input->GetLargestPossibleRegion().GetSize()[2];

	FloatImageType::Pointer red_weights = FloatImageType::New();
	FloatImageType::Pointer blue_weights = FloatImageType::New();

	FloatImageType::PointType origin1, origin2;
	origin1[0] = 0; origin2[0] = 0;
	origin1[1] = 0; origin2[1] = 0;
	origin1[2] = 0; origin2[2] = 0;
	red_weights->SetOrigin( origin1 ); blue_weights->SetOrigin( origin2 );

	FloatImageType::IndexType start1, start2;
	start1[0] = 0; start2[0] = 0; // first index on X
	start1[1] = 0; start2[1] = 0; // first index on Y
	start1[2] = 0; start2[2] = 0; // first index on Z
	FloatImageType::SizeType  size11, size22;
	size11[0] = size1; size22[0] = size1; // size along X
	size11[1] = size2; size22[1] = size2; // size along Y
	size11[2] = size3; size22[2] = size3; // size along Z
	FloatImageType::RegionType region1, region2;
	region1.SetSize( size11 );  region2.SetSize( size22 );
	region1.SetIndex( start1 ); region2.SetIndex( start2 );
	red_weights->SetRegions( region1 ); blue_weights->SetRegions( region2 );
	red_weights->Allocate();            blue_weights->Allocate();
	red_weights->FillBuffer(0);         blue_weights->FillBuffer(0);
	red_weights->Update();              blue_weights->Update();

	IteratorType iterator1 ( red_weights,      red_weights      ->GetRequestedRegion() );
	IteratorType iterator2 ( blue_weights,     blue_weights     ->GetRequestedRegion() );

	iterator1.GoToBegin();
	iterator2.GoToBegin();
	
	UcharIteratorType pix_buf3( red_image, red_image->GetRequestedRegion() );
	UcharIteratorType pix_buf2( lime_image, lime_image->GetRequestedRegion() );
	UcharIteratorType pix_buf1( intensity_image, intensity_image->GetRequestedRegion() );

	// ========== Create 3D Histogram =========
	pix_buf1.GoToBegin();	pix_buf2.GoToBegin();	pix_buf3.GoToBegin();
	for ( ;
	      !pix_buf1.IsAtEnd() && !pix_buf2.IsAtEnd() && !pix_buf3.IsAtEnd() && !iterator1.IsAtEnd() && !iterator2.IsAtEnd();
	      ++pix_buf1, ++pix_buf2, ++pix_buf3, ++iterator1, ++iterator2 )
	{
		dh::_RGB pixel = dh::_RGB( pix_buf1.Get(), pix_buf2.Get(), pix_buf3.Get() );

		float red_dist = 2.0 * ( 1.0 - slide_wt ) *
			dh::RGB_Classifier::RGB_euclidean_dist( pixel, atype.red, 1, 1, 1);
		//DMP(red_dist);

		float blue_dist = 2.0 * slide_wt *
			dh::RGB_Classifier::RGB_euclidean_dist( pixel, atype.blue, 1, 1, 1);
		//DMP(blue_dist);

		float bkgrnd_dist = dh::RGB_Classifier::RGB_euclidean_dist( pixel, atype.bkgrnd, 1, 1, 1);
		//DMP(bkgrnd_dist);

		//Cell_Class pixel_class;
		double s_plane_dist;

		float certainty;  	// 1 = certain; 0 = unknown

		s_plane_dist = dot ( ((dh::XYZ)pixel - split_point), decision_plane );
		if ( s_plane_dist >= 0 ) { // Pixel is "red"
			pixel_class = RED_CELL;
			double r_plane_dist = dh::dot ( ((dh::XYZ)pixel - (dh::XYZ)atype.red), red_plane );
			// This could be done with signum functions to make it faster...
			if ( r_plane_dist * red_sp_dist < 0 ) {
				certainty = 1.0;
			}
			else {
				certainty = s_plane_dist / ( s_plane_dist + fabs(r_plane_dist) );
			}
		}
		else
		{ // Pixel is "blue"
			pixel_class = BLUE_CELL;
			double b_plane_dist = dh::dot ( ((dh::XYZ)pixel - (dh::XYZ)atype.blue), blue_plane );
			// This could be done with signum functions to make it faster...
			if ( b_plane_dist * blue_sp_dist < 0 ){
				certainty = 1.0;
			}
			else{
				s_plane_dist = - s_plane_dist;
				certainty = s_plane_dist / ( s_plane_dist + fabs(b_plane_dist) );
			 }
		}
		//DMP(pixel_class);

		float bkgrnd_wt = 1 - certainty;


		if(  ( pixel_class == RED_CELL && (red_dist * (1 - bkgrnd_wt) < bkgrnd_dist * bkgrnd_wt))
			     ||(pixel_class == BLUE_CELL
			     && (blue_dist * (1 - bkgrnd_wt) < bkgrnd_dist * bkgrnd_wt))
			  ) 
		{
			certainty = pow(certainty, 3);
		}

		if( pixel_class != RED_CELL && pixel_class != BLUE_CELL ){
			iterator1.Set( blue_dist/(blue_dist+red_dist)*certainty );  //Function1
			iterator2.Set( red_dist/(blue_dist+red_dist) *certainty );  //Function2
		}
		else{
			iterator1.Set( (blue_dist+bkgrnd_dist)/(blue_dist+red_dist+bkgrnd_dist)*certainty );  //Function3
			iterator2.Set( (red_dist+bkgrnd_dist) /(blue_dist+red_dist+bkgrnd_dist)*certainty );  //Function4
		}
	}

	typedef itk::RescaleIntensityImageFilter< FloatImageType, UcharImageType > RescaleFlUcType;
	RescaleFlUcType::Pointer RescaleIntIO1 = RescaleFlUcType::New();
	RescaleIntIO1->SetOutputMaximum( 255 );
	RescaleIntIO1->SetOutputMinimum( 0 );
	RescaleIntIO1->SetInput( red_weights );
	RescaleIntIO1->Update();
	RescaleFlUcType::Pointer RescaleIntIO2 = RescaleFlUcType::New();
	RescaleIntIO2->SetOutputMaximum( 255 );
	RescaleIntIO2->SetOutputMinimum( 0 );
	RescaleIntIO2->SetInput( blue_weights );
	RescaleIntIO2->Update();
	red_wts  = RescaleIntIO1->GetOutput();
	blue_wts = RescaleIntIO2->GetOutput();

	typedef itk::ImageFileWriter< UcharImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "red_wts.tif" );
	writer->SetInput( red_wts );
	writer->Update();
	WriterType::Pointer writer1 = WriterType::New();
	writer1->SetFileName( "blue_wts.tif" );
	writer1->SetInput( blue_wts );
	writer1->Update();
}


void ColorSegmentation::go_best_dir( dh::_RGB& ma, bool& moved, const dh::_RGB& sa, const int r, int res )
{
	double new_val;
	double best_val = dh::NEG_HUGE;
	int best_dR, best_dG, best_dB;
	int x, y, z;
	
	if ( r % res != 0 )
	{ 
		std::cerr<<"go_best_dir: r % res != 0 may cause instability!"; 
	}
	if ( r <= 0 || res <=0 || res > r )
	{ 
		std::cerr<<"go_best_dir: Invalid parameters: r=" << r << " res=" << res;
	}

	for ( x = -r; x <= r; x+= res )
	{ 
		for ( y = -r; y <= r; y+= res )
		{ 
			for ( z = -r; z <= r; z+= res )
			{
				if ( // Point is inside blob                 v Blob Edge Threshold
			        hist->v(ma.R+x, ma.G+y, ma.B+z, true ) > 0
					  // and point is darker than background on each axis
				/*	  && ma.R+x < at.bkgrnd.R
					  && ma.G+y < at.bkgrnd.G */
				/*	  && ma.B+z < at.bkgrnd.B */   // Intensity, for RLI
				// 4/16/98- this is wrong place for this, with new seeding algorithm
				// add as an option?	 Or make flexible, to work on
				// light on dark OR dark on light?
				   )
			                     
				{ 
					new_val = dh::eval_state ( at.bkgrnd, dh::_RGB(ma.R+x, ma.G+y, ma.B+z), sa );
					if ( new_val > best_val )
					{ 
						best_val = new_val;
						best_dR = x;
						best_dG = y;
						best_dB = z;
					}
				}
			}
		}
	}
	if ( best_val != DBL_MIN )
	{ 
		ma.R += best_dR;
		ma.G += best_dG;
		ma.B += best_dB;
	}
	else
	{ 
		std::cerr<<"Optimization Algorithm Failure: point is outside blob!!!\n"; 
	}

	moved = ( best_dR == 0 || best_dG == 0 || best_dB == 0 ) ? false : true;
 }