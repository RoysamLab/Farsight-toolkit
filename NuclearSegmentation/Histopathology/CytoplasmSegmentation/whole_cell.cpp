#ifndef _CYTOPLASM_SEG_CXX_
#define _CYTOPLASM_SEG_CXX_

#define BIN_Image(a,b) bin_Image[(a)+(b)*size1]
#define inp_im_2D(a,b) INP_IM_2D[(a)+(b)*size1]
#define inp_im_2D1(a,b) INP_IM_2D1[(a)+(b)*size1]
#define nuc_im(a,b) NUC_IM[(a)+(b)*size1]
#define grad_imw(a,b) GRAD_IMW[(a)+(b)*(size1+2)]
/*
#define inp_im_2D(a,b) INP_IM_2D[(a)+(b)*size1]
#define inp_im_2D1(a,b) INP_IM_2D1[(a)+(b)*size1]
#define NUC(a,b) Nuc[(a)+(b)*size1]
#define NUCO(a,b) NucO[(a)+(b)*size1]
#define cytim(a,b) Cyt_Image[(a)+(b)*size1]
#define OP(a,b) Op[(a)+(b)*size1]
*/

#include "whole_cell.h"
//Constructor
WholeCellSeg::WholeCellSeg(){
	shift_bin = 0;
	num_levels = 1;
	num_levels_incl = 1;
	scaling = 255;
	mem_scaling = 1;
	use_mem_img = 0;
	bin_done = 0;
	seg_done = 0;
	nuc_im_set = 0;
	cyt_im_set = 0;
	mem_im_set = 0;
}

//Destructor
WholeCellSeg::~WholeCellSeg(){
}

//Set Parameters
void WholeCellSeg::set_parameters ( int *parameters ){
	shift_bin = *parameters; parameters++;
	num_levels = *parameters; parameters++;
	num_levels_incl = *parameters; parameters++;
	scaling = *parameters; parameters++;
	mem_scaling = *parameters; parameters++;
	use_mem_img = *parameters;
}

void WholeCellSeg::RunBinarization(){
	if( !nuc_im_set || !cyt_im_set ){
		std::cerr<<"Complete segmenting nuclei and set input imge before starting segmentation\n";
		return;
	}

	int size1=cyt_im_inp->GetLargestPossibleRegion().GetSize()[0];
	int size2=cyt_im_inp->GetLargestPossibleRegion().GetSize()[1];

	if( ( size1 != nuclab_inp->GetLargestPossibleRegion().GetSize()[0] ) || ( size2 != nuclab_inp->GetLargestPossibleRegion().GetSize()[1] ) ){
		std::cerr<<"The input images must be of the same size\n";
		return;
	}

	typedef unsigned short int UShortPixelType;

	typedef itk::ImageRegionIteratorWithIndex< UShortImageType > IteratorType;
	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;

	typedef itk::Statistics::ScalarImageToHistogramGenerator< IntImageType > ScalarImageToHistogramGeneratorType;
	typedef ScalarImageToHistogramGeneratorType::HistogramType HistogramType;
	typedef itk::OtsuMultipleThresholdsCalculator< HistogramType > CalculatorType;
	typedef itk::RescaleIntensityImageFilter< UShortImageType, IntImageType > RescaleUsIntType;
	typedef itk::RescaleIntensityImageFilter< UShortImageType, UShortImageType > RescaleUsUsType;
	typedef itk::BinaryThresholdImageFilter< IntImageType, UShortImageType >  ThreshFilterType;
	typedef itk::BinaryThresholdImageFilter< UShortImageType, UShortImageType >  ThresholdFilterType;
	typedef itk::OrImageFilter< UShortImageType, UShortImageType, UShortImageType > OrFilterType;
 	typedef itk::BinaryBallStructuringElement< UShortPixelType, 2 > StructuringElementType;
	typedef itk::BinaryErodeImageFilter< UShortImageType, UShortImageType, StructuringElementType > ErodeFilterType;
	typedef itk::BinaryDilateImageFilter< UShortImageType, UShortImageType, StructuringElementType > DilateFilterType;

	UShortImageType::Pointer intermediate_bin_im_out;

	unsigned char *in_Image;
	unsigned long int ind=0;

//Call Yousef's binarization method if the number of bin levels is < 2
	if( num_levels < 2 ){
		bin_Image = (unsigned short *) malloc (size1*size2*sizeof(unsigned short));
		for(int j=0; j<size2; ++j)
			for(int i=0; i<size1; ++i)
				BIN_Image(i,j)=255;
		in_Image = (unsigned char *) malloc (size1*size2);
		if( ( in_Image == NULL ) || ( bin_Image == NULL ) ){
			std::cerr << "Memory allocation for binarization of image failed\n";
			return;
		}

		RescaleUsUsType::Pointer rescaleususfilter = RescaleUsUsType::New();
		rescaleususfilter->SetInput( cyt_im_inp );
		rescaleususfilter->SetOutputMaximum( itk::NumericTraits<unsigned char>::max() );
		rescaleususfilter->SetOutputMinimum( 0 );
		rescaleususfilter->Update();
		UShortImageType::Pointer resc_cyt_im = UShortImageType::New();
		resc_cyt_im = rescaleususfilter->GetOutput();
		ConstIteratorType pix_buf1( resc_cyt_im, resc_cyt_im->GetRequestedRegion() );
		for ( pix_buf1.GoToBegin(); !pix_buf1.IsAtEnd(); ++pix_buf1, ++ind )
		in_Image[ind]=(unsigned char)(pix_buf1.Get());

		int ok = 0;
		ok = Cell_Binarization_2D(in_Image,bin_Image, size1, size2, shift_bin);
		free( in_Image );

		if( !ok ){
			std::cerr<<"Binarization Failed\n";
			return;
		}

//copy the output binary image into the ITK image
		intermediate_bin_im_out = UShortImageType::New();
		UShortImageType::PointType origin;
		origin[0] = 0;
		origin[1] = 0;
		intermediate_bin_im_out->SetOrigin( origin );

		UShortImageType::IndexType start;
		start[0] = 0;  // first index on X
		start[1] = 0;  // first index on Y
		UShortImageType::SizeType  size;
		size[0] = size1;  // size along X
		size[1] = size2;  // size along Y

		UShortImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );

		intermediate_bin_im_out->SetRegions( region );
		intermediate_bin_im_out->Allocate();
		intermediate_bin_im_out->FillBuffer(0);
		intermediate_bin_im_out->Update();

		unsigned short int dum,dum1;
		dum = 0;
		dum1 = USHRT_MAX;

		unsigned int asd,asd1; asd=0; asd1=0;
		IteratorType iterator ( intermediate_bin_im_out, intermediate_bin_im_out->GetRequestedRegion() );
		for(unsigned long int i=0; i<(unsigned long int(size1*size2)); ++i){
			if( bin_Image[i] )
			iterator.Set( dum1 );
			else
			iterator.Set( dum );
			++iterator;
		}
	}
//Call multi level binarization method if the number of bin levels is >= 2
	else{
		RescaleUsIntType::Pointer rescaleusintfilter = RescaleUsIntType::New();
		ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();
		rescaleusintfilter->SetInput( cyt_im_inp );
		rescaleusintfilter->SetOutputMaximum( itk::NumericTraits<unsigned short>::max() );
		rescaleusintfilter->SetOutputMinimum( 0 );
		rescaleusintfilter->Update();

		ThreshFilterType::Pointer threshfilter = ThreshFilterType::New();
		CalculatorType::Pointer calculator = CalculatorType::New();
		scalarImageToHistogramGenerator->SetNumberOfBins( 255 );
		calculator->SetNumberOfThresholds( num_levels );
		threshfilter->SetOutsideValue( (int)0 );
		threshfilter->SetInsideValue( (int)USHRT_MAX );
		scalarImageToHistogramGenerator->SetInput( rescaleusintfilter->GetOutput() );
		scalarImageToHistogramGenerator->Compute();
		calculator->SetInputHistogram( scalarImageToHistogramGenerator->GetOutput() );
		threshfilter->SetInput( rescaleusintfilter->GetOutput() );
		calculator->Update();
		const CalculatorType::OutputType &thresholdVector = calculator->GetOutput();
		CalculatorType::OutputType::const_iterator itNum = thresholdVector.begin();
		int lowerThreshold,upperThreshold;

		for( int i=0; i<(num_levels-num_levels_incl); ++i ) ++itNum;
		lowerThreshold = static_cast<unsigned short>(*itNum);
		upperThreshold = itk::NumericTraits<unsigned short>::max();

		threshfilter->SetLowerThreshold( lowerThreshold );
		threshfilter->SetUpperThreshold( upperThreshold );
		threshfilter->Update();

		intermediate_bin_im_out = UShortImageType::New();
		intermediate_bin_im_out = threshfilter->GetOutput();
	}

//Fill holes left by the nuclei
	ThresholdFilterType::Pointer binarythreshfilter = ThresholdFilterType::New();
	binarythreshfilter->SetInsideValue( (unsigned short int)USHRT_MAX );
	binarythreshfilter->SetOutsideValue( 0 );
	binarythreshfilter->SetLowerThreshold( 1 );
	binarythreshfilter->SetUpperThreshold( (unsigned short int)USHRT_MAX );
	binarythreshfilter->SetInput( nuclab_inp );
	OrFilterType::Pointer orfilter = OrFilterType::New();
	orfilter->SetInput1( binarythreshfilter->GetOutput() );
	orfilter->SetInput2( intermediate_bin_im_out );
//dialate and erode
	ErodeFilterType::Pointer  binaryErode  = ErodeFilterType::New();
	DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
	StructuringElementType  structuringElement;
	structuringElement.SetRadius( 3 );  // 3x3 structuring element
	structuringElement.CreateStructuringElement();
	binaryErode->SetKernel( structuringElement );
	binaryErode->SetErodeValue( (unsigned short int)USHRT_MAX );
	binaryDilate->SetDilateValue( (unsigned short int)USHRT_MAX );
	binaryDilate->SetKernel( structuringElement );
	binaryErode->SetInput( binaryDilate->GetOutput() );
	binaryDilate->SetInput( orfilter->GetOutput() );
//erode and dialate
	ErodeFilterType::Pointer  binaryErode1  = ErodeFilterType::New();
	DilateFilterType::Pointer binaryDilate1 = DilateFilterType::New();
	binaryErode1->SetKernel(  structuringElement );
	binaryErode1->SetErodeValue( (unsigned short int)USHRT_MAX );
	binaryDilate1->SetDilateValue( (unsigned short int)USHRT_MAX );
	binaryDilate1->SetKernel( structuringElement );
	binaryErode1->SetInput( binaryErode->GetOutput() );
	binaryDilate1->SetInput( binaryErode1->GetOutput() );
	binaryDilate1->Update();
//Get pointer to the final binary image and return it to calling function
	UShortImageType::Pointer image_bin = UShortImageType::New();
	image_bin = binaryDilate1->GetOutput();
	bin_im_out = image_bin;
	bin_done = 1;

//Update bin array
	ind=0;
	ConstIteratorType pix_buf3( bin_im_out, bin_im_out->GetRequestedRegion() );
	for ( pix_buf3.GoToBegin(); !pix_buf3.IsAtEnd(); ++pix_buf3, ++ind )
	bin_Image[ind]=(pix_buf3.Get());

	return;
}





void WholeCellSeg::RunSegmentation(){
	int size1=cyt_im_inp->GetLargestPossibleRegion().GetSize()[0];
	int size2=cyt_im_inp->GetLargestPossibleRegion().GetSize()[1];

	typedef itk::SmoothingRecursiveGaussianImageFilter< UShortImageType, UShortImageType > GaussianFilterType;
	typedef itk::RescaleIntensityImageFilter< UShortImageType, FltImageType > RescaleUSFltType;
	typedef itk::CastImageFilter< UShortImageType, IntImageType > CastUSIntType;
	typedef itk::CastImageFilter< IntImageType, UShortImageType > CastIntUSType;
	typedef itk::ImageRegionIteratorWithIndex< IntImageType > IteratorType1;
	typedef itk::ImageRegionConstIterator< FltImageType > ConstIteratorType1;
	typedef itk::ImageRegionConstIterator< UShortImageType > ConstIteratorType;
	typedef itk::MorphologicalWatershedFromMarkersImageFilter< IntImageType, IntImageType > WatershedFilterType;

//Get gradient image from cytoplasm image
	GaussianFilterType::Pointer  gaussianfilter = GaussianFilterType::New();
	gaussianfilter->SetSigma( 1.25 );
	gaussianfilter->SetInput( cyt_im_inp );
	gaussianfilter->Update();

//Rescale image
	RescaleUSFltType::Pointer rescaleusflt = RescaleUSFltType::New();
	rescaleusflt->SetOutputMaximum( scaling );
	rescaleusflt->SetOutputMinimum( 1 );
	rescaleusflt->SetInput( gaussianfilter->GetOutput() );
	rescaleusflt->Update();

//Get the rescaled gradient image from ITK into an array of known size and indexing system
	float *INP_IM_2D;
	FltImageType::Pointer grad_img = FltImageType::New();
	grad_img = rescaleusflt->GetOutput();
	INP_IM_2D = (float *) malloc (size1*size2*sizeof(float));
	if ( INP_IM_2D == NULL ){
		std::cerr << "Unable to allocate memory for 2D Gradient Image";
		return;
	}
	ConstIteratorType1 pix_buf( grad_img, grad_img->GetRequestedRegion() );
	long ind=0;
	for ( pix_buf.GoToBegin(); !pix_buf.IsAtEnd(); ++pix_buf, ++ind )
		INP_IM_2D[ind] = ( pix_buf.Get() );
	int testing=0;

	if( use_mem_img ){
		GaussianFilterType::Pointer  gaussianfilter1 = GaussianFilterType::New();
		gaussianfilter1->SetSigma( 1.25 );
		gaussianfilter1->SetInput( mem_im_inp );
		gaussianfilter1->Update();
		RescaleUSFltType::Pointer rescaleusflt1 = RescaleUSFltType::New();
		rescaleusflt1->SetOutputMaximum( scaling );
		rescaleusflt1->SetOutputMinimum( 1 );
		rescaleusflt1->SetInput( gaussianfilter1->GetOutput() );
		rescaleusflt1->Update();
		FltImageType::Pointer grad_img1 = FltImageType::New();
		grad_img1 = rescaleusflt1->GetOutput();
		float *INP_IM_2D1;
		INP_IM_2D1 = (float *) malloc (size1*size2*sizeof(float));
		if ( INP_IM_2D1 == NULL ){
			std::cerr << "Unable to allocate memory for 2D Membrane Image";
			return;
		}
		ConstIteratorType1 pix_buf2( grad_img1, grad_img1->GetRequestedRegion() );
		ind=0;
		for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2, ++ind )
			INP_IM_2D1[ind]=(pix_buf2.Get());
		for(int j=0; j<size2; j++)
			for(int i=0; i<size1; i++)
				inp_im_2D(i,j) = (inp_im_2D(i,j)/mem_scaling)*(inp_im_2D1(i,j)*mem_scaling);
		free( INP_IM_2D1 );
	}


//Get the nucleus label image from ITK into an array of known size and indexing system
	unsigned short *NUC_IM;
	NUC_IM = (unsigned short *) malloc (size1*size2*sizeof(unsigned short));
	if ( NUC_IM == NULL ){
		std::cerr << "Unable to allocate memory for Nucleus Label Image";
		return;
	}
	ConstIteratorType pix_buf2( nuclab_inp, nuclab_inp->GetRequestedRegion() );
	ind=0;
	for ( pix_buf2.GoToBegin(); !pix_buf2.IsAtEnd(); ++pix_buf2, ++ind )
		NUC_IM[ind]=(pix_buf2.Get());

//allocate memory for the gradient weighted distance map
	float *GRAD_IMW;
	GRAD_IMW = (float *) malloc ((size1+2)*(size2+2)*sizeof(float));
	if ( GRAD_IMW == NULL ){
		std::cerr << "Unable to allocate memory for Gradient Weighted Distance Image";
		return;
	}

//Create Gradient Weighted Distance Map
	float flt_mini = -1*FLT_MAX;
	for(int i=0; i<size1; i++)
		for(int j=0; j<size2; j++){
			if(!nuc_im(i,j)) grad_imw(i+1,j+1) = FLT_MAX;
			else grad_imw(i+1,j+1)=0;
		}

	for(int i=0; i<size1; i++)
		for(int j=0; j<size2; j++)
			if(!BIN_Image(i,j)) grad_imw(i+1,j+1)=flt_mini;

	free( NUC_IM );
	free( bin_Image );

	for(int i=0; i<(size1+2); i++){
		grad_imw(i,0)=flt_mini;
		grad_imw(i,size2+1)=flt_mini;
	}

	for(int i=0; i<(size2+2); i++){
		grad_imw(0,i)=flt_mini;
		grad_imw(size1+1,i)=flt_mini;
	}

	int ok;
 	ok = gradient_enhanced_distance_map( INP_IM_2D, GRAD_IMW, size1, size2);

	free( INP_IM_2D );

//Getting gradient weighted distance map into ITK array
	IntImageType::Pointer image2;
	image2 = IntImageType::New();
	IntImageType::PointType origint;
	origint[0] = 0;
	origint[1] = 0;
	image2->SetOrigin( origint );

	IntImageType::IndexType startt;
	startt[0] = 0;  // first index on X
	startt[1] = 0;  // first index on Y
	IntImageType::SizeType  sizet;
	sizet[0] = size1;  // size along X
	sizet[1] = size2;  // size along Y
	IntImageType::RegionType regiont;
	regiont.SetSize( sizet );
	regiont.SetIndex( startt );
	image2->SetRegions( regiont );
	image2->Allocate();
	image2->FillBuffer(0);
	image2->Update();
	//copy the output image into the ITK image
	IteratorType1 iteratort(image2,image2->GetRequestedRegion());
	for(int j=0; j<size2; j++){
		for(int i=0; i<size1; i++){
			iteratort.Set(grad_imw(i+1,j+1));
			++iteratort;
		}
	}
	free( GRAD_IMW );

	CastUSIntType::Pointer castUSIntfilter = CastUSIntType::New();
	castUSIntfilter->SetInput( nuclab_inp );
	castUSIntfilter->Update();
	IntImageType::Pointer nuclab_inp_int = IntImageType::New();
	nuclab_inp_int = castUSIntfilter->GetOutput();

//Seeded watershed to get the cytoplasm regions
	WatershedFilterType::Pointer watershedfilter = WatershedFilterType::New();
	watershedfilter->SetInput1( image2 );
	watershedfilter->SetInput2( nuclab_inp_int );
	watershedfilter->SetMarkWatershedLine( 1 );
	watershedfilter->Update();

	CastIntUSType::Pointer castIntUSfilter = CastIntUSType::New();
	castIntUSfilter->SetInput( watershedfilter->GetOutput() );
	castIntUSfilter->Update();
	seg_im_out = castIntUSfilter->GetOutput();

//Write the output for testing
/*	RescaleIntIOType::Pointer RescaleIntIO1 = RescaleIntIOType::New();
	RescaleIntIO1->SetOutputMaximum( 255 );
	RescaleIntIO1->SetOutputMinimum( 0 );
	RescaleIntIO1->SetInput( rescaleiointfilter->GetOutput() ); //watershedfilter->GetOutput() image1
	RescaleIntIO1->Update();
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "bin_info.tif" );
	writer->SetInput( RescaleIntIO1->GetOutput() );//RescaleIntIO1--finalO/P
	writer->Update();*/

//Get Array into IDL
/*	IntImageType::Pointer image_fin = IntImageType::New();
	image_fin = watershedfilter->GetOutput();
	ConstIteratorType3 pix_bufed3( image_fin, image_fin->GetRequestedRegion() );
	pix_bufed3.GoToBegin();

	for(int j=size2-1; j>=0; j--)
		for(int i=0; i<size1; i++){
			OP(i,j)=(pix_bufed3.Get());
			++pix_bufed3;
		}
*/
	seg_done = 1;
}

UShortImageType::Pointer WholeCellSeg::getBinPointer(){
	if( bin_done ) return bin_im_out;
	else return NULL;
}

UShortImageType::Pointer WholeCellSeg::getSegPointer(){
	if( seg_done ) return seg_im_out;
	else return NULL;
}

#endif