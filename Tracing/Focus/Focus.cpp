#include "Focus.h"

//Constructor
Focus::Focus(UCharImageType::Pointer input)
	: radius(5)
{
	imgIn = input;
	imgInColor = NULL;
	varImg = NULL;
}
Focus::Focus(RGBImageType::Pointer input)
	: radius(5)
{
	imgIn = NULL;
	imgInColor = input;
	varImg = NULL;
}

void Focus::MakeVarianceImage()
{

	UCharImageType::Pointer img;
	if(imgIn)
		img = imgIn;
	else
	{
		typedef itk::RGBToLuminanceImageFilter< RGBImageType, UCharImageType > ConvertFilterType;
		ConvertFilterType::Pointer convert = ConvertFilterType::New();
		convert->SetInput( imgInColor );
		convert->Update();
		img = convert->GetOutput();
	}

	typedef itk::VarianceImageFunction< UCharImageType2D > VarianceFunctionType;
	typedef itk::ExtractImageFilter< UCharImageType, UShortImageType2D > ExtractFilterType;

	std::cout << "Making Variance Image";

	int size1 = (int)img->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)img->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)img->GetLargestPossibleRegion().GetSize()[2];

	varImg = FloatImageType::New();

	FloatImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	origin[2] = 0;
	varImg->SetOrigin( origin );

	FloatImageType::IndexType start = {{ 0,0,0 }};
	FloatImageType::SizeType  size = {{ size1, size2, size3 }};
	FloatImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	varImg->SetRegions( region ); 
	varImg->Allocate();

	for(int z = 0; z < size3; ++z)
	{
		std::cout << ".";
		UCharImageType::IndexType index = { {0,0,z} };
		UCharImageType::SizeType size = { {size1, size2, 0} };
		UCharImageType::RegionType region;
		region.SetSize(size);
		region.SetIndex(index);

		ExtractFilterType::Pointer extract = ExtractFilterType::New();
		extract->SetInput( img );
		extract->SetExtractionRegion( region );
		extract->Update();

		UShortImageType2D::Pointer im2 = extract->GetOutput();

		typedef itk::MeanImageFilter<UShortImageType2D,UShortImageType2D> MeanFilterType;
		typedef itk::SquareImageFilter<UShortImageType2D,UShortImageType2D> SquareFilterType;
		typedef itk::SubtractImageFilter<UShortImageType2D,UShortImageType2D> SubtractFilterType;

		UShortImageType2D::SizeType radiussize = {{radius, radius}};
		MeanFilterType::Pointer mean1 = MeanFilterType::New();
		mean1->SetInput(im2);
		mean1->SetRadius(radiussize);

		SquareFilterType::Pointer sq1 = SquareFilterType::New();
		sq1->SetInput(mean1->GetOutput());

		SquareFilterType::Pointer sq2 = SquareFilterType::New();
		sq2->SetInput(im2);

		MeanFilterType::Pointer mean2 = MeanFilterType::New();
		mean2->SetInput(sq2->GetOutput());
		mean2->SetRadius(radiussize);

		SubtractFilterType::Pointer sb1 = SubtractFilterType::New();
		sb1->SetInput1(mean2->GetOutput());
		sb1->SetInput2(sq1->GetOutput());

		sb1->Update();

		typedef itk::ImageRegionIterator<FloatImageType> FloatIteratorType;
		typedef itk::ImageRegionIterator<UShortImageType2D> UShortIteratorType;

		FloatIteratorType fit(varImg,region);
		UShortIteratorType uit(sb1->GetOutput(),sb1->GetOutput()->GetLargestPossibleRegion());
		for(fit.GoToBegin(),uit.GoToBegin();!uit.IsAtEnd(); ++fit, ++uit)
		{
			//printf("%d ", int(uit.Get()));
			fit.Set(uit.Get());
		}


		/* // DONT USE THIS BECAUSE IT IS VERY SLOW
		VarianceFunctionType::Pointer varFunction = VarianceFunctionType::New();
		varFunction->SetInputImage( im2 );
		varFunction->SetNeighborhoodRadius( radius );

		for(int x=0; x<size1; ++x)
		{
		for(int y=0; y<size2; ++y)
		{
		UCharImageType2D::IndexType ind;
		ind[0] = x;
		ind[1] = y;
		float var = (float)varFunction->EvaluateAtIndex( ind );

		FloatImageType::IndexType ind3;
		ind3[0] = x;
		ind3[1] = y;
		ind3[2] = z;
		varImg->SetPixel(ind3, var);
		}
		}
		*/

	}

	typedef itk::ImageFileWriter<FloatImageType> FloatWriter;
	FloatWriter::Pointer fwriter = FloatWriter::New();
	fwriter->SetInput(varImg);
	fwriter->SetFileName("test_fwriter.mhd");
	fwriter->Update();
	/*
	//Rescale:
	typedef itk::RescaleIntensityImageFilter< FloatImageType, UCharImageType > RescaleType;
	RescaleType::Pointer rescale = RescaleType::New();
	rescale->SetOutputMaximum( 255 );
	rescale->SetOutputMinimum( 0 );
	rescale->SetInput( varImg );
	rescale->Update();
	*/
	std::cout << "done\n";
}

UCharImageType2D::Pointer Focus::MakeProjection()
{
	if(!imgIn)
		return NULL;

	int size1 = (int)imgIn->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)imgIn->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)imgIn->GetLargestPossibleRegion().GetSize()[2];

	UCharImageType2D::Pointer outImg = UCharImageType2D::New();
	UCharImageType2D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	outImg->SetOrigin( origin );

	UCharImageType2D::IndexType start = {{ 0,0 }};
	UCharImageType2D::SizeType  size = {{ size1, size2 }};
	UCharImageType2D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	outImg->SetRegions( region ); 
	outImg->Allocate();

	typedef itk::ImageRegionIteratorWithIndex< UCharImageType2D > IteratorType;
	IteratorType it( outImg, outImg->GetLargestPossibleRegion() );
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		UCharImageType2D::IndexType ind = it.GetIndex();
		int max = -1000;
		int mSlice = 0;
		for(int z=0; z<size3; ++z)
		{
			FloatImageType::IndexType ind3f;
			ind3f[0] = ind[0];
			ind3f[1] = ind[1];
			ind3f[2] = z;
			FloatImageType::PixelType pix = varImg->GetPixel(ind3f);
			if((int)pix > max)
			{
				max = (int)pix;
				mSlice = z;
			}
		}
		UCharImageType::IndexType ind3;
		ind3[0] = ind[0];
		ind3[1] = ind[1];
		ind3[2] = mSlice;
		it.Set( imgIn->GetPixel(ind3) );
	}

	return outImg;
}

std::pair<UCharImageType2D::Pointer,UShortImageType2D::Pointer> Focus::MakeProjection2()
{
	if(!imgIn)
		return std::pair<UCharImageType2D::Pointer,UShortImageType2D::Pointer>(NULL,NULL);

	int size1 = (int)imgIn->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)imgIn->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)imgIn->GetLargestPossibleRegion().GetSize()[2];

	UCharImageType2D::Pointer outImg = UCharImageType2D::New();
	UShortImageType2D::Pointer outDepth = UShortImageType2D::New();
	UCharImageType2D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	outImg->SetOrigin( origin );
	outDepth->SetOrigin( origin );

	UCharImageType2D::IndexType start = {{ 0,0 }};
	UCharImageType2D::SizeType  size = {{ size1, size2 }};
	UCharImageType2D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	outImg->SetRegions( region ); 
	outDepth->SetRegions ( region );
	outImg->Allocate();
	outDepth->Allocate();

	typedef itk::ImageRegionIteratorWithIndex< UCharImageType2D > IteratorType;
	IteratorType it( outImg, outImg->GetLargestPossibleRegion() );
	typedef itk::ImageRegionIteratorWithIndex < UShortImageType2D> IteratorType2;
	IteratorType2 it2( outDepth, outDepth->GetLargestPossibleRegion() );
	for(it.GoToBegin(),it2.GoToBegin(); !it.IsAtEnd(); ++it,++it2)
	{
		UCharImageType2D::IndexType ind = it.GetIndex();
		int max = -1000;
		int mSlice = 0;
		for(int z=0; z<size3; ++z)
		{
			FloatImageType::IndexType ind3f;
			ind3f[0] = ind[0];
			ind3f[1] = ind[1];
			ind3f[2] = z;
			FloatImageType::PixelType pix = varImg->GetPixel(ind3f);
			if((int)pix > max)
			{
				max = (int)pix;
				mSlice = z;
			}
		}
		UCharImageType::IndexType ind3;
		ind3[0] = ind[0];
		ind3[1] = ind[1];
		ind3[2] = mSlice;
		it.Set( imgIn->GetPixel(ind3) );
		it2.Set ( mSlice );
	}

	std::pair< UCharImageType2D::Pointer, UShortImageType2D::Pointer> retpair;
	retpair.first = outImg;
	retpair.second = outDepth;
	return retpair;
}


RGBImageType2D::Pointer Focus::MakeProjectionColor()
{
	if(!imgInColor)
		return NULL;

	int size1 = (int)imgInColor->GetLargestPossibleRegion().GetSize()[0];
	int size2 = (int)imgInColor->GetLargestPossibleRegion().GetSize()[1];
	int size3 = (int)imgInColor->GetLargestPossibleRegion().GetSize()[2];

	RGBImageType2D::Pointer outImg = RGBImageType2D::New();
	RGBImageType2D::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	outImg->SetOrigin( origin );

	RGBImageType2D::IndexType start = {{ 0,0 }};
	RGBImageType2D::SizeType  size = {{ size1, size2 }};
	RGBImageType2D::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	outImg->SetRegions( region ); 
	outImg->Allocate();

	typedef itk::ImageRegionIteratorWithIndex< RGBImageType2D > IteratorType;
	IteratorType it( outImg, outImg->GetLargestPossibleRegion() );
	for(it.GoToBegin(); !it.IsAtEnd(); ++it)
	{
		RGBImageType2D::IndexType ind = it.GetIndex();
		int max = -1000;
		int mSlice = 0;
		for(int z=0; z<size3; ++z)
		{
			FloatImageType::IndexType ind3f;
			ind3f[0] = ind[0];
			ind3f[1] = ind[1];
			ind3f[2] = z;
			FloatImageType::PixelType pix = varImg->GetPixel(ind3f);
			if((int)pix > max)
			{
				max = (int)pix;
				mSlice = z;
			}
		}
		RGBImageType::IndexType ind3;
		ind3[0] = ind[0];
		ind3[1] = ind[1];
		ind3[2] = mSlice;
		it.Set( imgInColor->GetPixel(ind3) );
	}

	return outImg;
}

std::vector<float> Focus::FindVariance(float x, float y, float ppi)
{
	int ix = (int)(x*ppi);
	int iy = (int)(y*ppi);

	return FindVariance(ix,iy);
}

std::vector<float> Focus::FindVariance(int x, int y)
{
	std::vector<float> variances;

	if(!varImg)
		return variances;

	int sizeX = (int)varImg->GetLargestPossibleRegion().GetSize()[0];
	int sizeY = (int)varImg->GetLargestPossibleRegion().GetSize()[1];
	int sizeZ = (int)varImg->GetLargestPossibleRegion().GetSize()[2];

	if(x >= sizeX || y >= sizeY)
		return variances;

	for(int z=0; z<sizeZ; ++z)
	{
		FloatImageType::IndexType ind;
		ind[0] = x;
		ind[1] = y;
		ind[2] = z;

		FloatImageType::PixelType val = varImg->GetPixel(ind);
		variances.push_back(val);
	}

	return variances;
}
/*
std::vector<float> Focus::FindVariance(int x, int y)
{
typedef itk::Image< UCharPixelType, 2 > UCharImageType2D;
typedef itk::VarianceImageFunction< UCharImageType2D > VarianceFunctionType;
typedef itk::ExtractImageFilter< UCharImageType, UCharImageType2D > ExtractFilterType;

int sizeX = (int)img->GetLargestPossibleRegion().GetSize()[0];
int sizeY = (int)img->GetLargestPossibleRegion().GetSize()[1];
int sizeZ = (int)img->GetLargestPossibleRegion().GetSize()[2];

std::vector<float> variances;

for(int z = 0; z < sizeZ; ++z)
{
UCharImageType::IndexType index = { 0,0,z };
UCharImageType::SizeType size = { sizeX, sizeY, 0 };
UCharImageType::RegionType region;
region.SetSize(size);
region.SetIndex(index);

ExtractFilterType::Pointer extract = ExtractFilterType::New();
extract->SetInput( img );
extract->SetExtractionRegion( region );
extract->Update();

UCharImageType2D::Pointer im2 = extract->GetOutput();

VarianceFunctionType::Pointer varFunction = VarianceFunctionType::New();
varFunction->SetInputImage( im2 );
varFunction->SetNeighborhoodRadius( radius );

UCharImageType2D::IndexType ind;
ind[0] = x;
ind[1] = y;

float var = (float)varFunction->EvaluateAtIndex( ind );
variances.push_back(var);
}

return variances;
}
*/