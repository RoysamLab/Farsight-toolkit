/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/

#include "ImageOperation.h"

ImageOperation::ImageOperation()
{
  SM = 0;
  SN = 0;
  SZ = 0;
  mask_set = false;
  display_set = false;
  u1 = 100;
  u2 = 0;
  sigma1 = sigma2 = 100;
  num_soma = 0;
  coding_method = 1;
}

void ImageOperation::ImComputeInitBackgroundModel()
{

  typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;

  SliceIteratorType it(I, I->GetRequestedRegion() );
  it.SetFirstDirection( 0 );
  it.SetSecondDirection( 1 );
  it.GoToBegin();

  u2 = 0;
  std::vector<float> intensity;
  while( !it.IsAtEnd() )
  {
   while ( !it.IsAtEndOfSlice() )
   {
     while ( !it.IsAtEndOfLine() )
     {
		intensity.push_back(it.Get());
		++it;
	 }
	 it.NextLine();
   }
   it.NextSlice();
  }

  for(std::vector<float>::iterator j=intensity.begin();j!=intensity.end();++j)     
  {
      u2 += *j;
  }
  u2 /= intensity.size();

  for(std::vector<float>::iterator j=intensity.begin();j!=intensity.end();++j)     
  {
      sigma2 += pow(*j - u2, 2);
  }

  sigma2 = sqrt(sigma2/intensity.size());

 
  std::cout<<"background model:"<<u2<<","<<sigma2<<std::endl;
}
void ImageOperation::ImComputeInitForegroundModel()
{
   typedef itk::NearestNeighborInterpolateImageFunction< 
                       ImageType, float>  InterpolatorType;
   InterpolatorType::Pointer interpolator = InterpolatorType::New();
   interpolator->SetInputImage(I);
   I_cu.clear();
	//consider all the seed points as initial foreground
   for( int j = 0; j < SeedPt.NP; j++ )
   {
	   ImageType::IndexType index;
	   index[0] = SeedPt.Pt[j].x;
	   index[1] = SeedPt.Pt[j].y;
	   index[2] = SeedPt.Pt[j].z;
	   float temp = interpolator->EvaluateAtIndex(index);
	   I_cu.push_back(temp);
	}
   u1 = 0;
   sigma1 = 0;
  for(std::vector<float>::iterator j=I_cu.begin();j!=I_cu.end();++j)     
  {
      u1 += *j;
  }
  u1 /= I_cu.size();
  for(std::vector<float>::iterator j=I_cu.begin();j!=I_cu.end();++j)     
  {
      sigma1 += pow(*j - u1, 2);
  }

  sigma1 = sqrt(sigma1/I_cu.size());
  std::cout<<"foreground model:"<<u1<<","<<sigma1<<std::endl;
  //I_cu.clear();

}
void ImageOperation::ImComputeForegroundModel(PointList3D Cu, std::vector<float> Ru)
{
   typedef itk::NearestNeighborInterpolateImageFunction< 
                       ImageType, float>  InterpolatorType;
   InterpolatorType::Pointer interpolator = InterpolatorType::New();
   interpolator->SetInputImage(I);

   Vector3D v1,v2,v3, vtemp;
   float pi = 3.1415926;
   int m = 8;
   Point3D temp_r_pt;

   for( int j = 0; j < Cu.NP; j++ )
   {
	   ImageType::IndexType index;
	   index[0] = Cu.Pt[j].x;
	   index[1] = Cu.Pt[j].y;
	   index[2] = Cu.Pt[j].z;
	   float temp = interpolator->EvaluateAtIndex(index);
	   I_cu.push_back(temp);


	   /*if( j == 0 )
	   {
	      v1.x = Cu.Pt[0].x - Cu.Pt[1].x;
          v1.y = Cu.Pt[0].y - Cu.Pt[1].y;
          v1.z = Cu.Pt[0].z - Cu.Pt[1].z;
          v1.ConvertUnit();
	   }
	   else
	   {
	   	  v1.x = Cu.Pt[j].x - Cu.Pt[j-1].x;
          v1.y = Cu.Pt[j].y - Cu.Pt[j-1].y;
          v1.z = Cu.Pt[j].z - Cu.Pt[j-1].z;
          v1.ConvertUnit();
	   }

       v2.x = -v1.z;
       v2.y = 0;
       v2.z = v1.x;
       v2.ConvertUnit();
       v3.x = 1;
       v3.y = -(pow(v1.x,2) + pow(v1.z,2))/(v1.x*v1.y + std::numeric_limits<float>::epsilon());
       v3.z = v1.z/(v1.x + std::numeric_limits<float>::epsilon());
       v3.ConvertUnit();

	    for( int k = 0; k < m; k++ )
	    {
	      float theta = (2 * pi * k)/m;
	      vtemp.x = v2.x * cos(theta) + v3.x * sin(theta);
	      vtemp.y = v2.y * cos(theta) + v3.y * sin(theta);
		  vtemp.z = v2.z * cos(theta) + v3.z * sin(theta);
		  vtemp.ConvertUnit();

		  vnl_vector<float> oj(3);
		  oj(0) = vtemp.x;
		  oj(1) = vtemp.y;
		  oj(2) = vtemp.z;

		  temp_r_pt.x = Cu.Pt[j].x + 1/2 * Ru[j] * vtemp.x;
		  temp_r_pt.y = Cu.Pt[j].y + 1/2 * Ru[j] * vtemp.y;
		  temp_r_pt.z = Cu.Pt[j].z + 1/2 * Ru[j] * vtemp.z;
		
		  if( temp_r_pt.check_out_of_range_3D(SM,SN,SZ) )
			continue;

          ProbImageType::IndexType temp_index; 
          temp_index[0] = temp_r_pt.x;
          temp_index[1] = temp_r_pt.y;
		  temp_index[2] = temp_r_pt.z;
		  float temp = interpolator->EvaluateAtIndex(temp_index);
		  I_cu.push_back(temp);
		}*/
   }

   u1 = 0;
   sigma1 = 0;
  for(std::vector<float>::iterator j=I_cu.begin();j!=I_cu.end();++j)     
  {
      u1 += *j;
  }
  u1 /= I_cu.size();
  for(std::vector<float>::iterator j=I_cu.begin();j!=I_cu.end();++j)     
  {
      sigma1 += pow(*j - u1, 2);
  }

  sigma1 = sqrt(sigma1/I_cu.size());

  std::cout<<"foreground model:"<<u1<<","<<sigma1<<std::endl;
   
}

IOImagePointer ImageOperation::ImBkSub(IOImagePointer In)
{
	std::cout<<"subtracting background..."<<std::endl;
	int radius = 10;
   
    typedef itk::SubtractImageFilter<IOImageType, IOImageType, IOImageType>
			filterType;
	filterType::Pointer filter;
	filter = filterType::New();

    filter->SetInput1(In);
	filter->SetInput2(ImGaussian_XY(In, radius));
	filter->Update();
	return filter->GetOutput();
}

void ImageOperation::ImRemoveSlice(int in)
{
   if( in == 0 )
	  return;

   SZ = I->GetLargestPossibleRegion().GetSize()[2];

   typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
   SliceIteratorType It( I, I->GetRequestedRegion() );
  
   It.SetFirstDirection( 0 );
   It.SetSecondDirection( 1 );
   It.GoToBegin();

   while( !It.IsAtEnd() )
   {
    while ( !It.IsAtEndOfSlice() )
    {
     while ( !It.IsAtEndOfLine() )
     {
	   if( It.GetIndex()[2] < in || It.GetIndex()[2] >= SZ - in )
		 It.Set(0);
       ++It;
	 }
     It.NextLine();
    }
    It.NextSlice();
   }
}

void ImageOperation::ImDisplayRead(const char *filename, int shrink_factor)
{
	/*//IOImagePointer I; 
	typedef itk::ImageFileReader<IORGBImageType> ReaderType;		
	ReaderType::Pointer reader = ReaderType::New();	

	reader->SetFileName(filename);
	IORGBImagePointer IIO = reader->GetOutput();
	reader->Update();

    typedef itk::CastImageFilter<IORGBImageType, RGBImageType> CasterType;

    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(IIO);

	caster->Update();

	RGBImagePointer IRGB3D = RGBImageType::New();
	IRGB3D = caster->GetOutput();

	typedef itk::ShrinkImageFilter< RGBImageType,RGBImageType > ShrinkFilterType;
    ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New(); 

    shrinkFilter->SetInput( IRGB3D );
    unsigned int dfactors[3] = { shrink_factor, shrink_factor, 1};
    shrinkFilter->SetShrinkFactors(dfactors);
    shrinkFilter->UpdateLargestPossibleRegion();
    IRGB3D = shrinkFilter->GetOutput();
    shrinkFilter->Update();
	IRGB = ImMaxProjection(IRGB3D);*/


	typedef itk::ImageFileReader<IOImageType> ReaderType1;		
	ReaderType1::Pointer reader1 = ReaderType1::New();	

	reader1->SetFileName(filename);
	IOImagePointer IIO1 = reader1->GetOutput();
	reader1->Update();

    typedef itk::CastImageFilter<IOImageType, ImageType> CasterType1;

    CasterType1::Pointer caster1 = CasterType1::New();
    caster1->SetInput(IIO1);

	caster1->Update();

	IDisplay = caster1->GetOutput();

	typedef itk::ShrinkImageFilter< ImageType, ImageType > ShrinkFilterType1;
    ShrinkFilterType1::Pointer shrinkFilter1 = ShrinkFilterType1::New(); 

    shrinkFilter1->SetInput( IDisplay );
    unsigned int dfactors1[3] = { shrink_factor, shrink_factor, 1};
    shrinkFilter1->SetShrinkFactors(dfactors1);
    shrinkFilter1->UpdateLargestPossibleRegion();
    IDisplay = shrinkFilter1->GetOutput();
    shrinkFilter1->Update();

	display_set = true;
}

void ImageOperation::ImSeriesReadWrite(std::vector< std::string > filenames, const char *save_name, int shrink_factor, bool sixteen_bit)
{
 if( sixteen_bit )
 {
  typedef itk::ImageSeriesReader< IOImageType1 >  ReaderType;
  typedef itk::ImageFileWriter< IOImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //reader->SetImageIO( itk::TIFFImageIO::New() );

  reader->SetFileNames( filenames  );
  reader->Update();

  //typedef itk::CastImageFilter<IOImageType1, IOImageType> CasterType;
  //CasterType::Pointer caster = CasterType::New();
  //caster->SetInput(reader->GetOutput());
  //caster->Update();

  typedef itk::RescaleIntensityImageFilter< IOImageType1, IOImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( reader->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  writer->SetFileName( save_name );
  writer->SetInput( ImShrink(rescale->GetOutput(), shrink_factor) );

  try 
  { 
    writer->Update(); 
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
  } 
 }
 else
 {
  typedef itk::ImageSeriesReader< IOImageType >  ReaderType;
  typedef itk::ImageFileWriter< IOImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  //reader->SetImageIO( itk::TIFFImageIO::New() );

  reader->SetFileNames( filenames  );
  reader->Update();

  //typedef itk::CastImageFilter<IOImageType1, IOImageType> CasterType;
  //CasterType::Pointer caster = CasterType::New();
  //caster->SetInput(reader->GetOutput());
  //caster->Update();

  typedef itk::RescaleIntensityImageFilter< IOImageType, IOImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( reader->GetOutput() );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();

  writer->SetFileName( save_name );
  writer->SetInput( ImShrink(rescale->GetOutput(), shrink_factor) );

  try 
  { 
    writer->Update(); 
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
  } 
 }
}

void ImageOperation::ImRead(const char *filename)
{
	//IOImagePointer I; 
	typedef itk::ImageFileReader<IOImageType> ReaderType;		
	ReaderType::Pointer reader = ReaderType::New();	

	reader->SetFileName(filename);
	IOImagePointer IIO = reader->GetOutput();
	reader->Update();

    typedef itk::CastImageFilter<IOImageType, ImageType> CasterType;

    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(IIO);

	caster->Update();
	I = caster->GetOutput();

	//ImSmoothing(1);
	//ImRemoveSlice(in);

    typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
    rescale->SetInput( I );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 255 );
	I = rescale->GetOutput();
    rescale->Update();
    
	SM = I->GetLargestPossibleRegion().GetSize()[0];
    SN = I->GetLargestPossibleRegion().GetSize()[1];
    SZ = I->GetLargestPossibleRegion().GetSize()[2];
}

void ImageOperation::ImRead_NoSmooth(const char *filename, int in)
{
	//IOImagePointer I; 
	typedef itk::ImageFileReader<IOImageType> ReaderType;		
	ReaderType::Pointer reader = ReaderType::New();	

	reader->SetFileName(filename);
	IOImagePointer IIO = reader->GetOutput();
	reader->Update();

    typedef itk::CastImageFilter<IOImageType, ImageType> CasterType;

    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(IIO);

	caster->Update();
	I = caster->GetOutput();

	//ImSmoothing();
	ImRemoveSlice(in);

    typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
    rescale->SetInput( I );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 255 );

    //I = ImGaussian(rescale->GetOutput(),1);
	I = rescale->GetOutput();
    rescale->Update();
}

void ImageOperation::ImSmoothing(int smoothing_scale)
{

	if( smoothing_scale != 0 )
	{
      I = ImGaussian(I,smoothing_scale);
	}

   /*typedef itk::LaplacianRecursiveGaussianImageFilter<
                        ImageType, ImageType >  FilterType;

	FilterType::Pointer laplacian = FilterType::New();
    laplacian->SetNormalizeAcrossScale( true );
    laplacian->SetInput( I );
    laplacian->SetSigma( 1 );
	laplacian->Update();
    I = laplacian->GetOutput();

   typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
   SliceIteratorType It( I, I->GetRequestedRegion() );
   It.SetFirstDirection( 0 );
   It.SetSecondDirection( 1 );
   It.GoToBegin();

   while( !It.IsAtEnd() )
   {
    while ( !It.IsAtEndOfSlice() )
    {
     while ( !It.IsAtEndOfLine() )
     {
	   if( It.Get() > 0)
		   It.Set(0);
	   if( It.Get() < 0)
		   It.Set(It.Get()*-1);

       ++It;
	 }
     It.NextLine();
    }
    It.NextSlice();
   }*/

    /*typedef itk::CastImageFilter<ImageType, ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(I);
    caster->Update();

    typedef itk::GradientAnisotropicDiffusionImageFilter<
               ProbImageType, ProbImageType >  FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput( ImGaussian(caster->GetOutput(), 1) );
	filter->SetNumberOfIterations( 5 );
    filter->SetTimeStep( 0.0625 );

    filter->Update();

	typedef itk::RescaleIntensityImageFilter< ProbImageType, ImageType> RescaleFilterType;
    RescaleFilterType::Pointer rescale = RescaleFilterType::New();
    rescale->SetInput( filter->GetOutput() );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 255 );

    I = rescale->GetOutput();
    rescale->Update();  */

}


void ImageOperation::ConvertReadImage()
{
   /*//ImagePointer IOut;
   typedef itk::CastImageFilter<IOImageType, ImageType> CasterType;
   CasterType::Pointer caster = CasterType::New();
   caster->SetInput(IIO);

   typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> RescaleFilterType;
   RescaleFilterType::Pointer rescale = RescaleFilterType::New();

   rescale->SetInput( caster->GetOutput() );
   rescale->SetOutputMinimum( 0 );
   rescale->SetOutputMaximum( 255 );

   I = rescale->GetOutput();
   rescale->Update();

   //return IOut; */
}

LabelImagePointer2D ImageOperation::ImMaxProjection(LabelImagePointer IInput)
{
   int dim = 2;

   typedef itk::MaximumProjectionImageFilter< LabelImageType, LabelImageType > FilterType;
   FilterType::Pointer filter = FilterType::New();
   filter->SetInput( IInput );
   filter->SetProjectionDimension( dim );
   // to be sure that the result is ok with several threads, even on a single
   // proc computer
   filter->SetNumberOfThreads( 2 );
   filter->Update();
   LabelImageType::SizeType inputSize = filter->GetOutput()->GetLargestPossibleRegion().GetSize();
   typedef itk::ExtractImageFilter< LabelImageType, LabelImageType2D > ExtractType;
   ExtractType::Pointer extract = ExtractType::New();
   extract->SetInput( filter->GetOutput() );
   LabelImageType::SizeType size;
   for(int i=0; i<=3; i++) 
   {
    if(i == dim) 
    {
     size[i] = 0;
    } 
    else 
    { 
     size[i] = inputSize[i];
    }
   }
   LabelImageType::IndexType idx;
   idx.Fill(0);
   LabelImageType::RegionType region;
   region.SetSize( size );
   region.SetIndex( idx );
   extract->SetExtractionRegion( region );
   extract->Update();
   LabelImagePointer2D I2D = extract->GetOutput();

   return I2D;
  /* typedef itk::ImageLinearIteratorWithIndex< LabelImageType2D > LinearIteratorType;
   typedef itk::ImageSliceConstIteratorWithIndex< LabelImageType > SliceIteratorType;
   unsigned int projectionDirection = 2;
   unsigned int i, j;
   unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
  {
   if (i != projectionDirection)
  {
   direction[j] = i;
   j++;
   }
  }

  LabelImageType2D::RegionType region;
  LabelImageType2D::RegionType::SizeType size;
  LabelImageType2D::RegionType::IndexType index;
  LabelImageType::RegionType requestedRegion = IInput->GetRequestedRegion();

  index[ direction[0] ] = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ] = requestedRegion.GetSize()[ direction[0] ];
  size[ 1- direction[0] ] = requestedRegion.GetSize()[ direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );
  LabelImagePointer2D I2D = LabelImageType2D::New();
  I2D->SetRegions( region );
  I2D->Allocate();

  SliceIteratorType inputIt( IInput, IInput->GetRequestedRegion() );
  LinearIteratorType outputIt( I2D, I2D->GetRequestedRegion() );
  inputIt.SetFirstDirection( direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
  {
   while ( ! outputIt.IsAtEndOfLine() )
  {
   outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
   ++outputIt;
  }
   outputIt.NextLine();
 }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
     outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
     ++inputIt;
     ++outputIt;
    }
   outputIt.NextLine();
   inputIt.NextLine();
   }

   outputIt.GoToBegin();
   inputIt.NextSlice();
  }
  
  return I2D; */
}


ImagePointer2D ImageOperation::ImMaxProjection(ImagePointer IInput)
{
   int dim = 2;

   typedef itk::MaximumProjectionImageFilter< ImageType, ImageType > FilterType;
   FilterType::Pointer filter = FilterType::New();
   filter->SetInput( IInput );
   filter->SetProjectionDimension( dim );
   // to be sure that the result is ok with several threads, even on a single
   // proc computer
   filter->SetNumberOfThreads( 2 );
   filter->Update();
   ImageType::SizeType inputSize = filter->GetOutput()->GetLargestPossibleRegion().GetSize();
   typedef itk::ExtractImageFilter< ImageType, ImageType2D > ExtractType;
   ExtractType::Pointer extract = ExtractType::New();
   extract->SetInput( filter->GetOutput() );
   ImageType::SizeType size;
   for(int i=0; i<=3; i++) 
   {
    if(i == dim) 
    {
     size[i] = 0;
    } 
    else 
    { 
     size[i] = inputSize[i];
    }
   }
   ImageType::IndexType idx;
   idx.Fill(0);
   ImageType::RegionType region;
   region.SetSize( size );
   region.SetIndex( idx );
   extract->SetExtractionRegion( region );
   extract->Update();
   ImagePointer2D I2D = extract->GetOutput();

   return I2D;

  /* typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
   typedef itk::ImageSliceConstIteratorWithIndex< ImageType > SliceIteratorType;
   unsigned int projectionDirection = 2;
   unsigned int i, j;
   unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
  {
   if (i != projectionDirection)
  {
   direction[j] = i;
   j++;
   }
  }

  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size;
  ImageType2D::RegionType::IndexType index;
  ImageType::RegionType requestedRegion = IInput->GetRequestedRegion();

  index[ direction[0] ] = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ] = requestedRegion.GetSize()[ direction[0] ];
  size[ 1- direction[0] ] = requestedRegion.GetSize()[ direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );
  ImagePointer2D I2D = ImageType2D::New();
  I2D->SetRegions( region );
  I2D->Allocate();

  SliceIteratorType inputIt( IInput, IInput->GetRequestedRegion() );
  LinearIteratorType outputIt( I2D, I2D->GetRequestedRegion() );
  inputIt.SetFirstDirection( direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
  {
   while ( ! outputIt.IsAtEndOfLine() )
  {
   outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
   ++outputIt;
  }
   outputIt.NextLine();
 }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
     outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
     ++inputIt;
     ++outputIt;
    }
   outputIt.NextLine();
   inputIt.NextLine();
   }

   outputIt.GoToBegin();
   inputIt.NextSlice();
  }
  
  return I2D; */
}

RGBImagePointer2D ImageOperation::ImMaxProjection1(ImagePointer IInput)
{
   typedef itk::ImageLinearIteratorWithIndex< ImageType2D > LinearIteratorType;
   typedef itk::ImageSliceConstIteratorWithIndex< ImageType > SliceIteratorType;
   unsigned int projectionDirection = 2;
   unsigned int i, j;
   unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
  {
   if (i != projectionDirection)
  {
   direction[j] = i;
   j++;
   }
  }

  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size;
  ImageType2D::RegionType::IndexType index;
  ImageType::RegionType requestedRegion = IInput->GetRequestedRegion();

  index[ direction[0] ] = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ] = requestedRegion.GetSize()[ direction[0] ];
  size[ 1- direction[0] ] = requestedRegion.GetSize()[ direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );
  ImagePointer2D I2D = ImageType2D::New();
  I2D->SetRegions( region );
  I2D->Allocate();

  SliceIteratorType inputIt( IInput, IInput->GetRequestedRegion() );
  LinearIteratorType outputIt( I2D, I2D->GetRequestedRegion() );
  inputIt.SetFirstDirection( direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
  {
   while ( ! outputIt.IsAtEndOfLine() )
  {
   outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
   ++outputIt;
  }
   outputIt.NextLine();
 }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
     outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
     ++inputIt;
     ++outputIt;
    }
   outputIt.NextLine();
   inputIt.NextLine();
   }

   outputIt.GoToBegin();
   inputIt.NextSlice();
  }
  
  typedef itk::CastImageFilter<ImageType2D, RGBImageType2D> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(I2D);
  caster->Update();

  return caster->GetOutput();
}


ProbImagePointer2D ImageOperation::ImMaxProjection(ProbImagePointer IInput)
{

   int dim = 2;

   typedef itk::MaximumProjectionImageFilter< ProbImageType, ProbImageType > FilterType;
   FilterType::Pointer filter = FilterType::New();
   filter->SetInput( IInput );
   filter->SetProjectionDimension( dim );
   // to be sure that the result is ok with several threads, even on a single
   // proc computer
   filter->SetNumberOfThreads( 2 );
   filter->Update();
   ProbImageType::SizeType inputSize = filter->GetOutput()->GetLargestPossibleRegion().GetSize();
   typedef itk::ExtractImageFilter< ProbImageType, ProbImageType2D > ExtractType;
   ExtractType::Pointer extract = ExtractType::New();
   extract->SetInput( filter->GetOutput() );
   ProbImageType::SizeType size;
   for(int i=0; i<=3; i++) 
   {
    if(i == dim) 
    {
     size[i] = 0;
    } 
    else 
    { 
     size[i] = inputSize[i];
    }
   }
   ProbImageType::IndexType idx;
   idx.Fill(0);
   ProbImageType::RegionType region;
   region.SetSize( size );
   region.SetIndex( idx );
   extract->SetExtractionRegion( region );
   extract->Update();
   ProbImagePointer2D I2D = extract->GetOutput();

   return I2D;

   /*typedef itk::ImageLinearIteratorWithIndex< ProbImageType2D > LinearIteratorType;
   typedef itk::ImageSliceConstIteratorWithIndex< ProbImageType > SliceIteratorType;
   unsigned int projectionDirection = 2;
   unsigned int i, j;
   unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
  {
   if (i != projectionDirection)
  {
   direction[j] = i;
   j++;
   }
  }

  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size;
  ImageType2D::RegionType::IndexType index;
  ImageType::RegionType requestedRegion = IInput->GetRequestedRegion();

  index[ direction[0] ] = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ] = requestedRegion.GetSize()[ direction[0] ];
  size[ 1- direction[0] ] = requestedRegion.GetSize()[ direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );
  ProbImagePointer2D I2D = ProbImageType2D::New();
  I2D->SetRegions( region );
  I2D->Allocate();

  SliceIteratorType inputIt( IInput, IInput->GetRequestedRegion() );
  LinearIteratorType outputIt( I2D, I2D->GetRequestedRegion() );
  inputIt.SetFirstDirection( direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
  {
   while ( ! outputIt.IsAtEndOfLine() )
  {
   outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
   ++outputIt;
  }
   outputIt.NextLine();
 }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
     outputIt.Set( vnl_math_max( outputIt.Get(), inputIt.Get() ));
     ++inputIt;
     ++outputIt;
    }
   outputIt.NextLine();
   inputIt.NextLine();
   }

   outputIt.GoToBegin();
   inputIt.NextSlice();
  }
  
  return I2D;*/
}

RGBImagePointer2D ImageOperation::ImMaxProjection(RGBImagePointer IInput)
{
   typedef itk::ImageLinearIteratorWithIndex< RGBImageType2D > LinearIteratorType;
   typedef itk::ImageSliceConstIteratorWithIndex< RGBImageType > SliceIteratorType;
   unsigned int projectionDirection = 2;
   unsigned int i, j;
   unsigned int direction[2];
  for (i = 0, j = 0; i < 3; ++i )
  {
   if (i != projectionDirection)
  {
   direction[j] = i;
   j++;
   }
  }

  ImageType2D::RegionType region;
  ImageType2D::RegionType::SizeType size;
  ImageType2D::RegionType::IndexType index;
  ImageType::RegionType requestedRegion = IInput->GetRequestedRegion();

  index[ direction[0] ] = requestedRegion.GetIndex()[ direction[0] ];
  index[ 1- direction[0] ] = requestedRegion.GetIndex()[ direction[1] ];
  size[ direction[0] ] = requestedRegion.GetSize()[ direction[0] ];
  size[ 1- direction[0] ] = requestedRegion.GetSize()[ direction[1] ];

  region.SetSize( size );
  region.SetIndex( index );
  RGBImagePointer2D I2D = RGBImageType2D::New();
  I2D->SetRegions( region );
  I2D->Allocate();

  SliceIteratorType inputIt( IInput, IInput->GetRequestedRegion() );
  LinearIteratorType outputIt( I2D, I2D->GetRequestedRegion() );
  inputIt.SetFirstDirection( direction[1] );
  inputIt.SetSecondDirection( direction[0] );
  outputIt.SetDirection( 1 - direction[0] );

  outputIt.GoToBegin();
  while ( ! outputIt.IsAtEnd() )
  {
   while ( ! outputIt.IsAtEndOfLine() )
  {
   outputIt.Set( itk::NumericTraits<unsigned short>::NonpositiveMin() );
   ++outputIt;
  }
   outputIt.NextLine();
 }

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
	 RGBImageType2D::PixelType temp;
	 temp[0] = vnl_math_max( outputIt.Get()[0], inputIt.Get()[0] );
	 temp[1] = vnl_math_max( outputIt.Get()[1], inputIt.Get()[1] );
	 temp[2] = vnl_math_max( outputIt.Get()[2], inputIt.Get()[2] );
     outputIt.Set( temp );
     ++inputIt;
     ++outputIt;
    }
   outputIt.NextLine();
   inputIt.NextLine();
   }

   outputIt.GoToBegin();
   inputIt.NextSlice();
  }
  
  return I2D;

}

ImagePointer ImageOperation::ImCopy()
{
    typedef itk::CastImageFilter<ImageType, ImageType> CasterType;

    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(I);
	caster->Update();
	return caster->GetOutput();
}

void ImageOperation::ImWrite(const char *filename, int image_id)
{
	IOImagePointer1 IIO;
	typedef itk::ImageFileWriter<IOImageType1> WriterType;
    typedef itk::CastImageFilter<ProbImageType,IOImageType1> ProbCasterType;
	ProbCasterType::Pointer prob_caster = ProbCasterType::New();
	typedef itk::CastImageFilter<ImageType, IOImageType1> CasterType;
	CasterType::Pointer caster = CasterType::New(); 

	if( image_id == 1 )
	{
	 caster->SetInput(I);
	 IIO = caster->GetOutput();
	 caster->Update();
	}
	else if( image_id == 2 )
	{
	 prob_caster->SetInput(IVessel);
	 IIO = prob_caster->GetOutput();
	 prob_caster->Update();
	}

	typedef itk::ImageFileWriter<IOImageType1> WriterType;

	WriterType::Pointer writer;	
	writer = WriterType::New();	

	writer->SetFileName(filename);
	writer->SetInput(IIO);
	writer->Update();
}

void ImageOperation::ImWrite_Soma(const char *filename)
{
	IOImagePointer1 IIO;
	typedef itk::ImageFileWriter<IOImageType1> WriterType;
	typedef itk::CastImageFilter<ImageType, IOImageType1> CasterType;
	CasterType::Pointer caster = CasterType::New(); 

	caster->SetInput(IMask);
	IIO = caster->GetOutput();
	caster->Update();

	typedef itk::ImageFileWriter<IOImageType1> WriterType;

	WriterType::Pointer writer;	
	writer = WriterType::New();	

	writer->SetFileName(filename);
	writer->SetInput(IIO);
	writer->Update();
}

void ImageOperation::ConvertWriteImage()	
{
	/*//IOImagePointer IIO;

	typedef itk::CastImageFilter<ImageType,IOImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
		
	caster->SetInput(I);
	IIO = caster->GetOutput();
	caster->Update();

	//return IIO;*/
}


IOImagePointer ImageOperation::ImInvert(IOImagePointer In)
{
   std::cout<<"inverting image..."<<std::endl;
   typedef itk::ImageSliceIteratorWithIndex< IOImageType > SliceIteratorType;
   SliceIteratorType It( In, In->GetRequestedRegion() );
   It.SetFirstDirection( 0 );
   It.SetSecondDirection( 1 );
   It.GoToBegin();

   while( !It.IsAtEnd() )
   {
    while ( !It.IsAtEndOfSlice() )
    {
     while ( !It.IsAtEndOfLine() )
     {
	   It.Set( 255 - It.Get() );
       ++It;
	 }
     It.NextLine();
    }
    It.NextSlice();
   }

   return In;
}

void ImageOperation::ImShrink(int shrink_factor)
{
   typedef itk::ShrinkImageFilter< ImageType,ImageType > ShrinkFilterType;
   ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New(); 

   shrinkFilter->SetInput( I );
   unsigned int dfactors[3] = { shrink_factor, shrink_factor, 1};
   shrinkFilter->SetShrinkFactors(dfactors);
   shrinkFilter->UpdateLargestPossibleRegion();
   I = shrinkFilter->GetOutput();
   shrinkFilter->Update();

   SM = I->GetLargestPossibleRegion().GetSize()[0];
   SN = I->GetLargestPossibleRegion().GetSize()[1];
   SZ = I->GetLargestPossibleRegion().GetSize()[2];
}

IOImagePointer ImageOperation::ImShrink(IOImagePointer In, int shrink_factor)
{
   typedef itk::ShrinkImageFilter< IOImageType,IOImageType > ShrinkFilterType;
   ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New(); 

   shrinkFilter->SetInput( In );
   unsigned int dfactors[3] = { shrink_factor, shrink_factor, 1};
   shrinkFilter->SetShrinkFactors(dfactors);
   shrinkFilter->UpdateLargestPossibleRegion();
   In = shrinkFilter->GetOutput();
   shrinkFilter->Update();

   return In;
}

void ImageOperation::ImRemove_RedSeeds(PointList3D Cu, std::vector<float> Ru)
{
  int expand_ratio = 1;
  typedef itk::VesselTubeSpatialObject<3> VesselTubeType;
  typedef itk::VesselTubeSpatialObjectPoint<3> VesselTubePointType;
  VesselTubeType::Pointer VesselTube = VesselTubeType::New();
  VesselTubeType::PointListType list;
  for( int i=0; i<Cu.NP; i++)
  {
   VesselTubePointType p;
   p.SetPosition(Cu.Pt[i].x, Cu.Pt[i].y, Cu.Pt[i].z);
   p.SetRadius(Ru[i] * expand_ratio);
   list.push_back(p);
  }
  VesselTube->GetProperty()->SetName("VesselTube");
  VesselTube->SetId(1);
  VesselTube->SetPoints(list);

  typedef itk::Point< float, 3 > PointType;

  for( int i = 0; i < SeedPt.NP; i++ )
  {
    if( visit_label(i) == 0 )
	{
      PointType p;
	  p[0] = SeedPt.Pt[i].x;
	  p[1] = SeedPt.Pt[i].y;
	  p[2] = SeedPt.Pt[i].z;
	  if( VesselTube->IsInside(p) )
	  {
		  visit_label(i) = 1;
	  }
	}
  }
}

void ImageOperation::ImRefresh_TracingImage()
{
  typedef itk::ImageSliceIteratorWithIndex< LabelImageType > SliceIteratorType;
   SliceIteratorType It( IL_Tracing, IL_Tracing->GetRequestedRegion() );
   It.SetFirstDirection( 0 );
   It.SetSecondDirection( 1 );
   It.GoToBegin();

   while( !It.IsAtEnd() )
   {
    while ( !It.IsAtEndOfSlice() )
    {
     while ( !It.IsAtEndOfLine() )
     {
	   It.Set( 0 );
       ++It;
	 }
     It.NextLine();
    }
    It.NextSlice();
   }
}

void ImageOperation::ImRefresh_LabelImage()
{
   typedef itk::CastImageFilter<ImageType, LabelImageType> CasterType1;
   CasterType1::Pointer caster1 = CasterType1::New();
   CasterType1::Pointer caster2 = CasterType1::New();
   caster1->SetInput(I);
   caster2->SetInput(I);
   IL = caster1->GetOutput();
   caster1->Update();
   IL_Tracing = caster2->GetOutput();
   caster2->Update();

   typedef itk::ImageSliceIteratorWithIndex< LabelImageType > SliceIteratorType;
   SliceIteratorType It1( IL, IL->GetRequestedRegion() );
   SliceIteratorType It2( IL_Tracing, IL_Tracing->GetRequestedRegion() );
   It1.SetFirstDirection( 0 );
   It1.SetSecondDirection( 1 );
   It1.GoToBegin();
   It2.SetFirstDirection( 0 );
   It2.SetSecondDirection( 1 );
   It2.GoToBegin();
   while( !It1.IsAtEnd() )
   {
    while ( !It1.IsAtEndOfSlice() )
    {
     while ( !It1.IsAtEndOfLine() )
     {
	   It1.Set( 0 );
	   It2.Set( 0 );
       ++It1;
	   ++It2;
	 }
     It1.NextLine();
	 It2.NextLine();
    }
    It1.NextSlice();
    It2.NextSlice();
   }
}

ImagePointer ImageOperation::ImRescale(ProbImagePointer IInput)
{
  typedef itk::RescaleIntensityImageFilter< ProbImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();

  rescale->SetInput( IInput );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();
  ImagePointer output_image = rescale->GetOutput();

  return output_image;
}

ImagePointer ImageOperation::ImRescale(LabelImagePointer IInput)
{
  typedef itk::RescaleIntensityImageFilter< LabelImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();

  rescale->SetInput( IInput );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();
  ImagePointer output_image = rescale->GetOutput();

  return output_image;
}

ImagePointer ImageOperation::ImRescale(ImagePointer IInput)
{
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( IInput );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();
  return rescale->GetOutput();
}

void ImageOperation::computeGVF(int noise_level, int num_iteration, int smoothing_scale)
{ 
  if( smoothing_scale == 0 )
  {
   typedef itk::GradientImageFilter<ProbImageType, float, float> GradientImageFilterType;
   GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();

   typedef itk::RescaleIntensityImageFilter< ImageType, ProbImageType> RescaleFilterType;
   RescaleFilterType::Pointer rescale = RescaleFilterType::New();
   rescale->SetInput( I );
   rescale->SetOutputMinimum( 0 );
   rescale->SetOutputMaximum( 1 );
   rescale->Update();

   gradientFilter->SetInput(rescale->GetOutput());
   //gradientFilter->SetInput(I);

   try
   {
    gradientFilter->Update();
   }
   catch( itk::ExceptionObject & err )
   {
    std::cerr << "Exception caught: " << err << std::endl;
   }

   //IG = gradientFilter->GetOutput();

   typedef itk::GradientVectorFlowImageFilter<GradientImageType, GradientImageType> GradientVectorFlowFilterType;
   GradientVectorFlowFilterType::Pointer GVFFilter = GradientVectorFlowFilterType::New();
 
    //if( num_iteration == 0 )
    //{
    // IGVF = gradientFilter->GetOutput();
    //}
    //else
    //{
     GVFFilter->SetInput(gradientFilter->GetOutput());
     GVFFilter->SetNoiseLevel(noise_level);
     GVFFilter->SetIterationNum(num_iteration);

     try
     {
      GVFFilter->Update();
     }
     catch( itk::ExceptionObject & err )
     {
      std::cerr << "Exception caught: " << err << std::endl;
     }
     IGVF = GVFFilter->GetOutput();
    //}

  }
  else
  {
   typedef itk::GradientRecursiveGaussianImageFilter<ProbImageType, GradientImageType> GradientImageFilterType;
   GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();

   typedef itk::RescaleIntensityImageFilter< ImageType, ProbImageType> RescaleFilterType;
   RescaleFilterType::Pointer rescale = RescaleFilterType::New();
   rescale->SetInput( I );
   rescale->SetOutputMinimum( 0 );
   rescale->SetOutputMaximum( 1 );
   rescale->Update();

   gradientFilter->SetSigma(smoothing_scale);
   gradientFilter->SetInput(rescale->GetOutput());

   try
   {
    gradientFilter->Update();
   }
   catch( itk::ExceptionObject & err )
   {
    std::cerr << "Exception caught: " << err << std::endl;
   }

   //IG = gradientFilter->GetOutput();

   typedef itk::GradientVectorFlowImageFilter<GradientImageType, GradientImageType> GradientVectorFlowFilterType;
   GradientVectorFlowFilterType::Pointer GVFFilter = GradientVectorFlowFilterType::New();

   GVFFilter->SetInput(gradientFilter->GetOutput());
   GVFFilter->SetNoiseLevel(noise_level);
   GVFFilter->SetIterationNum(num_iteration);

   try
   {
    GVFFilter->Update();
   }
   catch( itk::ExceptionObject & err )
   {
    std::cerr << "Exception caught: " << err << std::endl;
   }
   IGVF = GVFFilter->GetOutput();
  }

}

void ImageOperation::normalizeGVF()
{
  typedef itk::ImageSliceIteratorWithIndex< GradientImageType > SliceIteratorType;
  SliceIteratorType inputIt( IGVF, IGVF->GetRequestedRegion() );
  inputIt.SetFirstDirection( 0 );
  inputIt.SetSecondDirection( 1 );
  inputIt.GoToBegin();
 while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
	 if(inputIt.Get().GetNorm() != 0)
	 {
	   inputIt.Set( inputIt.Get()/( inputIt.Get().GetNorm() + std::numeric_limits<float>::epsilon() ) );
	 }
     ++inputIt;
    }
   inputIt.NextLine();
   }
   inputIt.NextSlice();
  }

 /* typedef itk::ImageDuplicator< GradientImageType > DuplicatorType;
  DuplicatorType::Pointer Duplicator = DuplicatorType::New();
  Duplicator->SetInputImage(IGVF);
  Duplicator->Update();
  IGVF_Norm = Duplicator->GetOutput();


  typedef itk::CastImageFilter<ImageType,ProbImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();

  caster->SetInput(I);
  caster->Update();
  GVF_Mag = caster->GetOutput();


  typedef itk::ImageSliceIteratorWithIndex< GradientImageType > SliceIteratorType;
  typedef itk::ImageSliceIteratorWithIndex< ProbImageType > SliceIteratorType1;

  SliceIteratorType inputIt( IGVF, IGVF->GetRequestedRegion() );
  SliceIteratorType outputIt( IGVF_Norm, IGVF_Norm->GetRequestedRegion() );
  SliceIteratorType1 outputIt1( GVF_Mag, GVF_Mag->GetRequestedRegion() );
  inputIt.SetFirstDirection( 0 );
  inputIt.SetSecondDirection( 1 );
  outputIt.SetFirstDirection( 0 );
  outputIt.SetSecondDirection( 1 );
  outputIt1.SetFirstDirection( 0 );
  outputIt1.SetSecondDirection( 1 );

  inputIt.GoToBegin();
  outputIt.GoToBegin();
  outputIt1.GoToBegin();

  while( !inputIt.IsAtEnd() )
  {
   while ( !inputIt.IsAtEndOfSlice() )
  {
    while ( !inputIt.IsAtEndOfLine() )
   {
	 if(inputIt.Get().GetNorm() == 0)
	 {
	   outputIt.Set( inputIt.Get());
	 }
	 else
	 {
	   outputIt.Set( inputIt.Get()/inputIt.Get().GetNorm() );
	 }

     
     outputIt1.Set( inputIt.Get().GetNorm() );
     ++inputIt;
     ++outputIt;
	 ++outputIt1;
    }
   inputIt.NextLine();
   outputIt.NextLine();
   outputIt1.NextLine();
   }
   inputIt.NextSlice();
   outputIt.NextSlice();
   outputIt1.NextSlice();
  } */
}

void ImageOperation::SeedAdjustment(int iter_num)
{
   
   normalizeGVF();

   if( iter_num == 0 )
   {
      //sort seeds by their saliency
     PointList3D New_SeedPt;

     vnl_vector<double> saliency(SeedPt.NP);
     for( int i = 0; i < SeedPt.NP; i++)
     {
	  GradientImageType::IndexType index; 
	  SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
	  index[0] = ceil(SeedPt.Pt[i].x);
	  index[1] = ceil(SeedPt.Pt[i].y);
	  index[2] = ceil(SeedPt.Pt[i].z);
	  saliency(i) = IVessel->GetPixel(index);
     }

     for( unsigned int i = 0; i < saliency.size(); i++)
     {
      int index = saliency.arg_max();
	  New_SeedPt.AddPt(SeedPt.Pt[index]);
	  saliency(index) = -1;
     } 

      SeedPt = New_SeedPt;

	  visit_label.set_size(SeedPt.GetSize());
	  visit_label.fill(0);
	   return;
   }

   //int iter_num = 100;

   typedef itk::VectorLinearInterpolateImageFunction< 
                       GradientImageType, float >  GradientInterpolatorType;

   GradientInterpolatorType::Pointer interpolator = GradientInterpolatorType::New();
   interpolator->SetInputImage(IGVF);


   //move seeds along gradient vector flow
   int j = 0;

   Point3D temp_pt;

   while( j < iter_num )
   {
	   for( int i = 0; i < SeedPt.NP; i++ )
	   {
		   GradientImageType::IndexType index; 
		  // std::cout<<SeedPt.Pt[i].x<<","<<SeedPt.Pt[i].y<<","<<SeedPt.Pt[i].z<<std::endl;
		   SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
		   index[0] = (SeedPt.Pt[i].x);
		   index[1] = (SeedPt.Pt[i].y);
		   index[2] = (SeedPt.Pt[i].z);
		    //std::cout<<index[0]<<","<<index[1]<<","<<index[2]<<std::endl;
		   GradientPixelType gradient = interpolator->EvaluateAtIndex(index);
           //std::cout<<"check point 20"<<std::endl;
		   //GradientPixelType gradient = IGVF->GetPixel(index);
		   //std::cout<<gradient[0]<<","<<gradient[1]<<","<<gradient[2]<<std::endl;
		   //Point3D temp_pt(gradient[0],gradient[1],gradient[2]);
		   temp_pt.x = gradient[0];
           temp_pt.y = gradient[1];
		   temp_pt.z = gradient[2];
           SeedPt.Pt[i] = SeedPt.Pt[i] + temp_pt;
	   }
	   j++;
   }

   /*//filter out seeds in the background
   int i = 0;
   while( i < SeedPt.NP )
   {
	  ImageType::IndexType index; 

	  SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);

	  index[0] = ceil(SeedPt.Pt[i].x);
	  index[1] = ceil(SeedPt.Pt[i].y);
	  index[2] = ceil(SeedPt.Pt[i].z);
	  
      if( VBW->GetPixel(index) == 0)
	  {    
		  SeedPt.RemovePt(i);
		  continue;
      }
	  i++;
   }*/


   seed_centroid();

   //sort seeds by their saliency
   PointList3D New_SeedPt;

   vnl_vector<double> saliency(SeedPt.NP);
   for( int i = 0; i < SeedPt.NP; i++)
   {
	  GradientImageType::IndexType index; 
	  SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
	  index[0] = ceil(SeedPt.Pt[i].x);
	  index[1] = ceil(SeedPt.Pt[i].y);
	  index[2] = ceil(SeedPt.Pt[i].z);
	  saliency(i) = IVessel->GetPixel(index);
   }

   for( unsigned int i = 0; i < saliency.size(); i++)
   {
      int index = saliency.arg_max();
	  New_SeedPt.AddPt(SeedPt.Pt[index]);
	  saliency(index) = -1;
   } 

    SeedPt = New_SeedPt;

	visit_label.set_size(SeedPt.GetSize());
	visit_label.fill(0);

} 

void ImageOperation::SeedDetection(float th, int detection_method, int seed_radius)
{
   
   SM = I->GetLargestPossibleRegion().GetSize()[0];
   SN = I->GetLargestPossibleRegion().GetSize()[1];
   SZ = I->GetLargestPossibleRegion().GetSize()[2];

   SeedType seed;
   SamplePointer seeds_candidate = SampleType::New();
   seeds_candidate->SetMeasurementVectorSize( 3 );
   
   //int max_num_seed = 100000;
   //SeedPt.SetN(max_num_seed);

   //int remove_radius = 5;

   Point3D temp_pt;

   int rad = 1;

   typedef itk::CastImageFilter<ImageType,LabelImageType> CasterType;
   CasterType::Pointer caster = CasterType::New();
   caster->SetInput(I);
   caster->Update();
   LabelImagePointer I_Seed = caster->GetOutput();
 
  typedef itk::NeighborhoodIterator< ProbImageType > NeighborhoodIteratorType;
  typedef itk::NeighborhoodIterator< ImageType > NeighborhoodIteratorType1;
  typedef itk::NeighborhoodIterator< GradientImageType > NeighborhoodIteratorType2;

  NeighborhoodIteratorType::RadiusType radius;
  NeighborhoodIteratorType1::RadiusType radius1;
  NeighborhoodIteratorType2::RadiusType radius2;
  radius.Fill(rad);
  radius1.Fill(rad);
  radius2.Fill(rad);

  NeighborhoodIteratorType it( radius, IVessel, IVessel->GetRequestedRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
   if( detection_method == 0 ) //Vesselness and Ridgeness
   {
     if(it.GetCenterPixel() >= th )
	 {
		  seed[0] = it.GetIndex()[0];
		  seed[1] = it.GetIndex()[1];
		  seed[2] = it.GetIndex()[2];
	    if( seed_radius != 0 )
		{
		  if( !SeedSparsify( seeds_candidate, seed, seed_radius ) )
		   seeds_candidate->PushBack( seed );
		}
		else
		{
		   seeds_candidate->PushBack( seed );
		}
	 } 
   }
   else if( detection_method == 1) //Local Maxima Seeds
   {
   
     ProbImageType::PixelType max = it.GetCenterPixel();
     for(unsigned i = 0; i < it.Size(); i++)
     {
      if( it.GetPixel(i) > max )
      {
       max = it.GetPixel(i);
      }
	 }
	 if(max == it.GetCenterPixel()) 
	 {
	    //SeedPt.AddPt(it.GetIndex()[0],it.GetIndex()[1],it.GetIndex()[2]);
		  seed[0] = it.GetIndex()[0];
		  seed[1] = it.GetIndex()[1];
		  seed[2] = it.GetIndex()[2];
	      if( !SeedSparsify( seeds_candidate, seed, seed_radius ) )
		  seeds_candidate->PushBack( seed );
	 }
   }
  
  }

    //SeedPt.SetN( seeds_candidate->Size() );
    SeedPt.RemoveAllPts();
    for( unsigned int i = 0; i < seeds_candidate->Size(); i++ )
	{
		SeedPt.AddPt( seeds_candidate->GetMeasurementVector(i)[0], seeds_candidate->GetMeasurementVector(i)[1], seeds_candidate->GetMeasurementVector(i)[2] );	
	}

	visit_label.set_size(SeedPt.GetSize());
	visit_label.fill(0);
  
}

bool ImageOperation::SeedSparsify(SamplePointer seeds_candidate, SeedType query_point, int radius)
{
    bool removal = false;
    Point3D temp_pt1, temp_pt2;
	temp_pt1.x = query_point[0];
    temp_pt1.y = query_point[1];
	temp_pt1.z = query_point[2];

	for( unsigned int i = 0; i < seeds_candidate->Size(); i++ )
	{
	   SeedType seed;
	   seed = seeds_candidate->GetMeasurementVector(i);
	   temp_pt2.x = seed[0];
	   temp_pt2.y = seed[1];
	   temp_pt2.z = seed[2];
	   if( temp_pt1.GetDistTo(temp_pt2) < radius )
	   {
	     removal = true;
		 break;
	   }
	}

	/*typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

    treeGenerator->SetSample( seeds_candidate );
    treeGenerator->SetBucketSize( 16 );
    treeGenerator->Update();

	typedef TreeGeneratorType::KdTreeType TreeType;
    TreeType::Pointer tree = treeGenerator->GetOutput();

	TreeType::InstanceIdentifierVectorType neighbors;
    tree->Search( query_point, (double)radius, neighbors );
      
    if( neighbors.size() > 0 )
	{
	   removal = true;
	} */
	return removal;
}

void ImageOperation::seed_centroid()
{
   typedef itk::CastImageFilter<ImageType,LabelImageType> CasterType;
   CasterType::Pointer caster = CasterType::New();
   caster->SetInput(I);
   caster->Update();
   LabelImagePointer I_Seed = caster->GetOutput();

   typedef itk::ImageSliceIteratorWithIndex< LabelImageType > SliceIteratorType;
   SliceIteratorType inputIt( I_Seed, I_Seed->GetRequestedRegion() );
   inputIt.SetFirstDirection( 0 );
   inputIt.SetSecondDirection( 1 );
   inputIt.GoToBegin();
   while( !inputIt.IsAtEnd() )
   {
    while ( !inputIt.IsAtEndOfSlice() )
   {
     while ( !inputIt.IsAtEndOfLine() )
    {
      inputIt.Set(0);
      ++inputIt;
     }
    inputIt.NextLine();
    }
    inputIt.NextSlice();
   }

   for( int i = 0; i < SeedPt.GetSize(); i++ )
   {
	   LabelImageType::IndexType index;
	   SeedPt.Pt[i].check_out_of_range_3D(SM,SN,SZ);
      /*for( int j = -1; j <= 1; j++ )
	  {
	    for( int k = -1; k <= 1; k++ )
		{
		  for( int z = -1; z <= 1; z++)
		  {
		  	index[0] = SeedPt.Pt[i].x + j;
	        index[1] = SeedPt.Pt[i].y + k;
	        index[2] = SeedPt.Pt[i].z + z;

			if( index[0] < 0 || index[1] < 0 || index[2] < 0 || 
	        index[0] >= SM || index[1] >= SN || index[2] >= SZ )
			   continue;

		     I_Seed->SetPixel(index, 1);
		  }
		}
	  }*/

	  index[0] = SeedPt.Pt[i].x;
	  index[1] = SeedPt.Pt[i].y;
	  index[2] = SeedPt.Pt[i].z;
	  I_Seed->SetPixel(index,1);
   }

  // Set up a connected components filter to label the binary objects.
  typedef itk::ConnectedComponentImageFilter< LabelImageType, LabelImageType > ConnectedComponentType;
  ConnectedComponentType::Pointer connectedComponentFilter = ConnectedComponentType::New();
  connectedComponentFilter->SetInput( I_Seed );

  // Relabel the components in order of size.
  typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelType;
  RelabelType::Pointer relabeler = RelabelType::New();
  relabeler->SetInput( connectedComponentFilter->GetOutput() );

  typedef itk::LabelGeometryImageFilter< LabelImageType > LabelGeometryType;
  LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();
  labelGeometryFilter->SetInput( relabeler->GetOutput() );

  labelGeometryFilter->Update();


  LabelGeometryType::LabelPointType index_temp;
  int labelValue = labelGeometryFilter->GetNumberOfLabels()-1;

  Point3D temp;
  //PointList3D new_SeedPt(labelValue+1);
  SeedPt.RemoveAllPts();
  for( int i = 1; i <= labelValue; i++)
  {
    index_temp = labelGeometryFilter->GetCentroid(i);
	temp.x = index_temp[0];
	temp.y = index_temp[1];
	temp.z = index_temp[2];
	SeedPt.AddPt(temp);
  }

}
void ImageOperation::SetCodingMethod(int in)
{
   coding_method = in + 1;
}

void ImageOperation::ImCoding( PointList3D Cu, std::vector<float> Ru, int snake_id, bool tracing )
{

   int coding_radius = 1;
   //int coding_method = 2;
   //std::cout<<"coding_method"<<coding_method<<std::endl;

   typedef itk::VesselTubeSpatialObject<3> VesselTubeType;
   typedef itk::VesselTubeSpatialObjectPoint<3> VesselTubePointType;
   VesselTubeType::Pointer VesselTube = VesselTubeType::New();
   VesselTubeType::PointListType list;
  

   if( coding_method == 2 )
   {
     if( tracing )
	 {
       for( int i=3; i<Cu.NP-3; i++)
       {
        VesselTubePointType p;
        p.SetPosition(Cu.Pt[i].x, Cu.Pt[i].y, Cu.Pt[i].z);
        p.SetRadius(Ru[i]);
        list.push_back(p);
       }
       VesselTube->GetProperty()->SetName("VesselTube");
       VesselTube->SetId(1);
       VesselTube->SetPoints(list);
	 }
	 else
	 {
	   for( int i=0; i<Cu.NP; i++)
       {
        VesselTubePointType p;
        p.SetPosition(Cu.Pt[i].x, Cu.Pt[i].y, Cu.Pt[i].z);
        p.SetRadius(Ru[i]);
        list.push_back(p);
       }
       VesselTube->GetProperty()->SetName("VesselTube");
       VesselTube->SetId(1);
       VesselTube->SetPoints(list);
	 }
   }


   if( tracing )
   {
	 if( Cu.NP <= 5 )
		 return;

    if( 1 )
    {
	 for( int i = 3; i < Cu.GetSize() - 3; i++ )
	 {
		 LabelImageType::IndexType index;
		 index[0] = Cu.Pt[i].x;
		 index[1] = Cu.Pt[i].y;
		 index[2] = Cu.Pt[i].z;
		 //IL->SetPixel( index, 1 );
		 for( int ix = -coding_radius; ix <= coding_radius; ix++ )
		 {
		   for( int iy = -coding_radius; iy <= coding_radius; iy++ )
		   {
		     for( int iz = -coding_radius; iz <= coding_radius; iz++ )
			 { 
               LabelImageType::IndexType new_index;
			   Point3D temp_pt;
			   temp_pt.x = Cu.Pt[i].x + ix;
               temp_pt.y = Cu.Pt[i].y + iy;
			   temp_pt.z = Cu.Pt[i].z + iz;
			   temp_pt.check_out_of_range_3D(SM,SN,SZ);
			   new_index[0] = temp_pt.x;
		       new_index[1] = temp_pt.y;
		       new_index[2] = temp_pt.z;
			   IL_Tracing->SetPixel( new_index, 1 );
			 }
		   }
		 }
	 }
	}
	/*else
	{
		    typedef itk::Point< float, 3 > PointType;
            typedef itk::ImageSliceIteratorWithIndex< LabelImageType > SliceIteratorType;
            SliceIteratorType It( IL_Tracing, IL_Tracing->GetRequestedRegion() );
			It.SetFirstDirection( 0 );
            It.SetSecondDirection( 1 );
            It.GoToBegin();
            while( !It.IsAtEnd() )
            {
              while ( !It.IsAtEndOfSlice() )
              {
                while ( !It.IsAtEndOfLine() )
                {
	             PointType p;

				 p[0] = It.GetIndex()[0];
                 p[1] = It.GetIndex()[1];
				 p[2] = It.GetIndex()[2];
				 if( VesselTube->IsInside(p) )
				    It.Set(1);

                 ++It;
				}
                 It.NextLine();
			  }
               It.NextSlice();
			}
	}*/
   }
   else
   {
    if( coding_method == 1 )
	{
	 for( int i = 0; i < Cu.GetSize(); i++ )
	 {
		 LabelImageType::IndexType index;
		 index[0] = Cu.Pt[i].x;
		 index[1] = Cu.Pt[i].y;
		 index[2] = Cu.Pt[i].z;
		 //IL->SetPixel( index, 1 );
		 for( int ix = -coding_radius; ix <=coding_radius; ix++ )
		 {
		   for( int iy = -coding_radius; iy <=coding_radius; iy++ )
		   {
		     for( int iz = -coding_radius; iz <=coding_radius; iz++ )
			 { 
               LabelImageType::IndexType new_index;
			   Point3D temp_pt;
			   temp_pt.x = Cu.Pt[i].x + ix;
               temp_pt.y = Cu.Pt[i].y + iy;
			   temp_pt.z = Cu.Pt[i].z + iz;
			   temp_pt.check_out_of_range_3D(SM,SN,SZ);
			   new_index[0] = temp_pt.x;
		       new_index[1] = temp_pt.y;
		       new_index[2] = temp_pt.z;
			   IL->SetPixel( new_index, snake_id );
			 }
		   }
		 }
	 }
	}
	else
	{
		    typedef itk::Point< float, 3 > PointType;
	        typedef itk::ImageSliceIteratorWithIndex< LabelImageType > SliceIteratorType;
			typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType1;
            SliceIteratorType It( IL, IL->GetRequestedRegion() );
			SliceIteratorType1 It1( I, I->GetRequestedRegion() );
			It.SetFirstDirection( 0 );
            It.SetSecondDirection( 1 );
            It.GoToBegin();
			It1.SetFirstDirection( 0 );
            It1.SetSecondDirection( 1 );
            It1.GoToBegin();
			std::vector<float> b_int;
			std::vector<float> f_int;
            while( !It.IsAtEnd() )
            {
              while ( !It.IsAtEndOfSlice() )
              {
                while ( !It.IsAtEndOfLine() )
                {
	             PointType p;
				 p[0] = It.GetIndex()[0];
                 p[1] = It.GetIndex()[1];
				 p[2] = It.GetIndex()[2];
                 
				 if( VesselTube->IsInside(p) )
				 {
				   if( It.Get() == 0 ) 
				    It.Set( snake_id );
				   if( snake_id == 0 )
				    It.Set( snake_id );
				 }
                 int new_int = It.Get();
				 if( new_int == 0 )
					 b_int.push_back(It1.Get());
				 else
				     f_int.push_back(It1.Get());

                 ++It;
				 ++It1;
				}
                 It.NextLine();
				 It1.NextLine();
			  }
               It.NextSlice();
			   It1.NextSlice();
			}
			//update the background and foreground models
			u1 = sigma1 = u2 = sigma2 = 0;
			  for(std::vector<float>::iterator j=b_int.begin();j!=b_int.end();++j)     
              {
               u2 += *j;
              }
			  for(std::vector<float>::iterator k=f_int.begin();k!=f_int.end();++k)     
              {
			   u1 += *k;
              }
              u2 /= b_int.size();
			  u1 /= f_int.size();

             for(std::vector<float>::iterator j=b_int.begin();j!=b_int.end();++j)       
             {
              sigma2 += pow(*j - u2, 2);
             }
			 for(std::vector<float>::iterator k=f_int.begin();k!=f_int.end();++k)       
             {
			  sigma1 += pow(*k - u1, 2);
             }
             sigma2 = sqrt(sigma2/b_int.size());
			 sigma1 = sqrt(sigma1/f_int.size());
			 std::cout<<"models:"<<u1<<","<<sigma1<<","<<u2<<","<<sigma2<<std::endl;
	}
   }
}


ImagePointer2D ImageOperation::extract_one_slice_yz(ImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ImageType, ImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ImageType::SizeType size = inputRegion.GetSize();

   size[0] = 0;
   ImageType::IndexType start = inputRegion.GetIndex();
   ImageType::RegionType desiredRegion;

   start[0] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   //std::cout<<"index:"<<index<<std::endl;
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ImagePointer2D I2D = filter->GetOutput();

   typedef itk::PermuteAxesImageFilter<ImageType2D> PermuterType;
   PermuterType::PermuteOrderArrayType order;
   order[0] = 1;
   order[1] = 0;
   PermuterType::Pointer Permuter = PermuterType::New();
   Permuter->SetInput(I2D);
   Permuter->SetOrder(order);
   Permuter->Update();
   I2D = Permuter->GetOutput();
	 
   return I2D;
}

ProbImagePointer2D ImageOperation::extract_one_slice_yz(ProbImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ProbImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ProbImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ProbImageType, ProbImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ProbImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ProbImageType::SizeType size = inputRegion.GetSize();

   size[0] = 0;
   ProbImageType::IndexType start = inputRegion.GetIndex();
   ProbImageType::RegionType desiredRegion;

   start[0] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   //std::cout<<"index:"<<index<<std::endl;
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ProbImagePointer2D I2D = filter->GetOutput();

   typedef itk::PermuteAxesImageFilter<ProbImageType2D> PermuterType;
   PermuterType::PermuteOrderArrayType order;
   order[0] = 1;
   order[1] = 0;
   PermuterType::Pointer Permuter = PermuterType::New();
   Permuter->SetInput(I2D);
   Permuter->SetOrder(order);
   Permuter->Update();
   I2D = Permuter->GetOutput();

   return I2D;
}


ImagePointer2D ImageOperation::extract_one_slice_xz(ImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ImageType, ImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ImageType::SizeType size = inputRegion.GetSize();
   size[1] = 0;
   ImageType::IndexType start = inputRegion.GetIndex();
   ImageType::RegionType desiredRegion;

   start[1] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ImagePointer2D I2D = filter->GetOutput();
   return I2D;
}

ProbImagePointer2D ImageOperation::extract_one_slice_xz(ProbImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ProbImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ProbImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ProbImageType, ProbImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ProbImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ProbImageType::SizeType size = inputRegion.GetSize();
   size[1] = 0;
   ProbImageType::IndexType start = inputRegion.GetIndex();
   ProbImageType::RegionType desiredRegion;

   start[1] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   //std::cout<<"index:"<<index<<std::endl;
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ProbImagePointer2D I2D = filter->GetOutput();
   return I2D;
}


ImagePointer2D ImageOperation::extract_one_slice(ImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ImageType, ImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ImageType::SizeType size = inputRegion.GetSize();
   size[2] = 0;
   ImageType::IndexType start = inputRegion.GetIndex();
   ImageType::RegionType desiredRegion;

   start[2] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ImagePointer2D I2D = filter->GetOutput();
   return I2D;
}

ProbImagePointer2D ImageOperation::extract_one_slice(ProbImagePointer I_input, int index)
{
   typedef itk::ImageDuplicator< ProbImageType > DuplicatorType;
   DuplicatorType::Pointer Duplicator = DuplicatorType::New();
   Duplicator->SetInputImage(I_input);
   Duplicator->Update();
   ProbImagePointer I_temp = Duplicator->GetOutput();

   typedef itk::ExtractImageFilter< ProbImageType, ProbImageType2D > FilterType;
   FilterType::Pointer filter = FilterType::New();
   ProbImageType::RegionType inputRegion = I_temp->GetLargestPossibleRegion();
   ProbImageType::SizeType size = inputRegion.GetSize();
   size[2] = 0;
   ProbImageType::IndexType start = inputRegion.GetIndex();
   ProbImageType::RegionType desiredRegion;

   start[2] = index;
   desiredRegion.SetSize( size );
   desiredRegion.SetIndex( start );
	
   //std::cout<<"index:"<<index<<std::endl;
   filter->SetExtractionRegion(desiredRegion);
   filter->SetInput(I_temp);

   filter->Update();
   ProbImagePointer2D I2D = filter->GetOutput();
   return I2D;
}

void ImageOperation::ImFastMarchingAnimation(int threshold)
{

	  /*typedef itk::BinaryThresholdImageFilter< ProbImageType, 
                        ImageType>    ThresholdingFilterType;
      ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
      thresholder->SetLowerThreshold( 0.0 );
      thresholder->SetUpperThreshold( threshold  );
      thresholder->SetOutsideValue(  0  );
      thresholder->SetInsideValue(  1 );
	  thresholder->SetInput( IDist );
      thresholder->Update();

     ISeg = thresholder->GetOutput();*/


	 /*typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
     SliceIteratorType it1( BW_temp, BW_temp->GetRequestedRegion() );
	 SliceIteratorType it2( ISeg, ISeg->GetRequestedRegion() );

	 it1.SetFirstDirection( 0 );
     it1.SetSecondDirection( 1 );
	 it2.SetFirstDirection( 0 );
     it2.SetSecondDirection( 1 );

	 it1.GoToBegin();
     it2.GoToBegin();

    while( !it1.IsAtEnd() )
    {
	  while ( !it1.IsAtEndOfSlice() )
      {
	 	while( !it1.IsAtEndOfLine())
	 	{
	       if( it1.Get() == 1 && it2.Get() == 0 )
	 		 it2.Set(1);
 
   	      ++it1;
          ++it2;
 		}
		 it1.NextLine();
		 it2.NextLine();
	  }
	 	 it1.NextSlice();
		 it2.NextSlice();
    }*/
}

void ImageOperation::clear_IMask()
{
     typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
	 SliceIteratorType it2( IMask, IMask->GetRequestedRegion() );

	 it2.SetFirstDirection( 0 );
     it2.SetSecondDirection( 1 );
     it2.GoToBegin();

    while( !it2.IsAtEnd() )
    {
	  while ( !it2.IsAtEndOfSlice() )
      {
	 	while( !it2.IsAtEndOfLine())
	 	{
	      it2.Set(0);
          ++it2;
 		}
		 it2.NextLine();
	  }
		 it2.NextSlice();
    } 
}

void ImageOperation::ImMasking(int shrink_factor)
{	
   typedef itk::ShrinkImageFilter< ImageType,ImageType > ShrinkFilterType;
   ShrinkFilterType::Pointer shrinkFilter = ShrinkFilterType::New(); 

   shrinkFilter->SetInput( IMask );
   unsigned int dfactors[3] = { shrink_factor, shrink_factor, 1};
   shrinkFilter->SetShrinkFactors(dfactors);
   shrinkFilter->UpdateLargestPossibleRegion();
   IMask = shrinkFilter->GetOutput();
   shrinkFilter->Update();


   typedef itk::ImageSliceIteratorWithIndex< ImageType > SliceIteratorType;
   SliceIteratorType It( I, I->GetRequestedRegion() );
   SliceIteratorType It1( IMask, IMask->GetRequestedRegion() );
   It.SetFirstDirection( 0 );
   It.SetSecondDirection( 1 );
   It.GoToBegin();
   It1.SetFirstDirection( 0 );
   It1.SetSecondDirection( 1 );
   It1.GoToBegin();

   while( !It.IsAtEnd() )
   {
    while ( !It.IsAtEndOfSlice() )
    {
      while ( !It.IsAtEndOfLine() )
      {
	   if( It1.Get() != 0 )
		   It.Set(0);
       ++It;
	   ++It1;
	  }
      It.NextLine();
	  It1.NextLine();
	}
    It.NextSlice();
	It1.NextSlice();
   }

   mask_set = true;

   //extract the centroids
  typedef itk::ConnectedComponentImageFilter< ImageType, LabelImageType > ConnectedComponentType;
  ConnectedComponentType::Pointer connectedComponentFilter = ConnectedComponentType::New();
  connectedComponentFilter->SetInput( IMask );

  // Relabel the components in order of size.
  typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelType;
  RelabelType::Pointer relabeler = RelabelType::New();
  relabeler->SetInput( connectedComponentFilter->GetOutput() );
  relabeler->SetMinimumObjectSize(100);

  typedef itk::LabelGeometryImageFilter< LabelImageType > LabelGeometryType;
  LabelGeometryType::Pointer labelGeometryFilter = LabelGeometryType::New();
  labelGeometryFilter->SetInput( relabeler->GetOutput() );

  labelGeometryFilter->Update();

  ISoma = relabeler->GetOutput();

  LabelGeometryType::LabelPointType index_temp;
  int labelValue = labelGeometryFilter->GetNumberOfLabels()-1;

  num_soma = labelValue;

  //compute the centroids of soma
  Point3D temp;
  //PointList3D new_SeedPt(labelValue+1);
  Centroid.RemoveAllPts();
  for( int i = 1; i <= labelValue; i++)
  {
    index_temp = labelGeometryFilter->GetCentroid(i);
	temp.x = index_temp[0];
	temp.y = index_temp[1];
	temp.z = index_temp[2];
	Centroid.AddPt(temp);
  }

  //compute the voronoi map
  typedef itk::DanielssonDistanceMapImageFilter< LabelImageType, LabelImageType > DistanceMapFilterType;
  DistanceMapFilterType::Pointer DMfilter = DistanceMapFilterType::New();
  DMfilter->SetInput( ISoma );
  DMfilter->Update();
  IVoronoi = DMfilter->GetVoronoiMap();
  
}



ImagePointer ImageOperation::ImGaussian(ImagePointer I_In, int sigma)
{
	//int sigma = 1;

    typedef itk::RecursiveGaussianImageFilter <ImageType, ImageType> FilterType;

    FilterType::Pointer filterX = FilterType::New();
    FilterType::Pointer filterY = FilterType::New();
	FilterType::Pointer filterZ = FilterType::New();
    filterX->SetDirection(0); // x direction
    filterY->SetDirection(1); // y direction
	filterZ->SetDirection(2); 
    filterX->SetOrder(FilterType::ZeroOrder);
    filterY->SetOrder(FilterType::ZeroOrder);
	filterZ->SetOrder(FilterType::ZeroOrder);
    filterX->SetNormalizeAcrossScale(false);
    filterY->SetNormalizeAcrossScale(false);
	filterZ->SetNormalizeAcrossScale(false);
    filterX->SetSigma(sigma);
    filterY->SetSigma(sigma);
	filterZ->SetSigma(sigma);

    filterX->SetInput( I_In );
    filterY->SetInput(filterX->GetOutput() );
	filterZ->SetInput(filterY->GetOutput() );
    filterZ->Update();
    ImagePointer Out = filterZ->GetOutput();
	Out->DisconnectPipeline();
    return Out;
}

ProbImagePointer ImageOperation::ImGaussian(ProbImagePointer I_In, int sigma)
{
    typedef itk::RecursiveGaussianImageFilter <ProbImageType, ProbImageType> FilterType;

    FilterType::Pointer filterX = FilterType::New();
    FilterType::Pointer filterY = FilterType::New();
	FilterType::Pointer filterZ = FilterType::New();
    filterX->SetDirection(0); // x direction
    filterY->SetDirection(1); // y direction
	filterZ->SetDirection(2); 
    filterX->SetOrder(FilterType::ZeroOrder);
    filterY->SetOrder(FilterType::ZeroOrder);
	filterZ->SetOrder(FilterType::ZeroOrder);
    filterX->SetNormalizeAcrossScale(false);
    filterY->SetNormalizeAcrossScale(false);
	filterZ->SetNormalizeAcrossScale(false);
    filterX->SetSigma(sigma);
    filterY->SetSigma(sigma);
	filterZ->SetSigma(sigma);

    filterX->SetInput( I_In );
    filterY->SetInput(filterX->GetOutput() );
	filterZ->SetInput(filterY->GetOutput() );
    filterZ->Update();
    ProbImagePointer Out = filterZ->GetOutput();
	Out->DisconnectPipeline();
    return Out;
}

IOImagePointer ImageOperation::ImGaussian_XY(IOImagePointer I_In, int sigma)
{
    typedef itk::RecursiveGaussianImageFilter <IOImageType, IOImageType> FilterType;

    FilterType::Pointer filterX = FilterType::New();
    FilterType::Pointer filterY = FilterType::New();
    filterX->SetDirection(0); // x direction
    filterY->SetDirection(1); // y direction
    filterX->SetOrder(FilterType::ZeroOrder);
    filterY->SetOrder(FilterType::ZeroOrder);
    filterX->SetNormalizeAcrossScale(false);
    filterY->SetNormalizeAcrossScale(false);
    filterX->SetSigma(sigma);
    filterY->SetSigma(sigma);

    filterX->SetInput( I_In );
    filterY->SetInput(filterX->GetOutput() );
    filterY->Update();
    IOImagePointer Out = filterY->GetOutput();
	Out->DisconnectPipeline();
    return Out;
}


void ImageOperation::ComputeGVFVesselness()
{
  double FrangiAlpha = 0.5;
  double FrangiBeta = 0.5;
  double FrangiC = 10;
  double A = 2 * pow(FrangiAlpha,2);
  double B = 2 * pow(FrangiBeta,2);
  double C = 2 * pow(FrangiC,2);

  typedef itk::CastImageFilter<ImageType,ProbImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(I);
  caster->Update();
  IVessel = caster->GetOutput();

  typedef itk::CastImageFilter<GradientImageType,GradientImageType> CasterType1;
  CasterType1::Pointer caster2 = CasterType1::New();
  caster2->SetInput(IGVF);
  caster2->Update();
  V1 = caster2->GetOutput();

  /*CasterType::Pointer caster1 = CasterType::New();
  caster1->SetInput(I);
  caster1->Update();
  Vx1 = caster1->GetOutput();


  CasterType::Pointer caster2 = CasterType::New();
  caster2->SetInput(I);
  caster2->Update();
  Vy1 = caster2->GetOutput();


  CasterType::Pointer caster3 = CasterType::New();
  caster3->SetInput(I);
  caster3->Update();
  Vz1 = caster3->GetOutput();*/


  typedef itk::RecursiveGaussianImageFilter<
  ProbImageType, ProbImageType > FilterType;
  //typedef itk::GradientImageFilter<ProbImageType, float, float> FilterType;
  ProbImagePointer Dx = extract_one_component(0,IGVF);
  ProbImagePointer Dy = extract_one_component(1,IGVF);
  ProbImagePointer Dz = extract_one_component(2,IGVF); 
  

  /*FilterType::Pointer filter_x = FilterType::New();
  filter_x->SetInput(Dx);
  filter_x->Update();
   Dxx = extract_one_component(0,filter_x->GetOutput());
   Dxy = extract_one_component(1,filter_x->GetOutput());
   Dxz = extract_one_component(2,filter_x->GetOutput());

  FilterType::Pointer filter_y = FilterType::New();
  filter_y->SetInput(Dy);
  filter_y->Update();
   Dyy = extract_one_component(1,filter_y->GetOutput());
   Dyz = extract_one_component(2,filter_y->GetOutput());

  FilterType::Pointer filter_z = FilterType::New();
  filter_z->SetInput(Dz);
  filter_z->Update();
   Dzz = extract_one_component(2,filter_z->GetOutput());*/

  //Dxx
  FilterType::Pointer filter1 = FilterType::New();
  filter1->SetDirection(0);
  filter1->SetOrder(FilterType::FirstOrder);
  filter1->SetInput(Dx);
  filter1->SetSigma(1);
  filter1->Update();
  ProbImagePointer Dxx = filter1->GetOutput();

  //Dxy
  FilterType::Pointer filter2 = FilterType::New();
  filter2->SetDirection(1);
  filter2->SetOrder(FilterType::FirstOrder);
  filter2->SetInput(Dx);
  filter2->SetSigma(1);
  filter2->Update();
  ProbImagePointer Dxy = filter2->GetOutput();

  //Dxz
  FilterType::Pointer filter3 = FilterType::New();
  filter3->SetDirection(2);
  filter3->SetOrder(FilterType::FirstOrder);
  filter3->SetInput(Dx);
  filter3->SetSigma(1);
  filter3->Update();
  ProbImagePointer Dxz = filter3->GetOutput();

  //Dyy
  FilterType::Pointer filter4 = FilterType::New();
  filter4->SetDirection(1);
  filter4->SetOrder(FilterType::FirstOrder);
  filter4->SetInput(Dy);
  filter4->SetSigma(1);
  filter4->Update();
  ProbImagePointer Dyy = filter4->GetOutput();

  //Dyz
  FilterType::Pointer filter5 = FilterType::New();
  filter5->SetDirection(2);
  filter5->SetOrder(FilterType::FirstOrder);
  filter5->SetInput(Dy);
  filter5->SetSigma(1);
  filter5->Update();
  ProbImagePointer Dyz = filter5->GetOutput();

  //Dzz
  FilterType::Pointer filter6 = FilterType::New();
  filter6->SetDirection(2);
  filter6->SetOrder(FilterType::FirstOrder);
  filter6->SetInput(Dz);
  filter6->SetSigma(1);
  filter6->Update();
  ProbImagePointer Dzz = filter6->GetOutput(); 


  typedef itk::ImageSliceIteratorWithIndex< ProbImageType > SliceIteratorType;
  typedef itk::ImageSliceIteratorWithIndex< GradientImageType > GradientSliceIteratorType;
  typedef itk::ImageSliceIteratorWithIndex< ImageType > ISliceIteratorType;

  SliceIteratorType inputIt1( Dxx, Dxx->GetRequestedRegion() );
  inputIt1.SetFirstDirection( 0 );
  inputIt1.SetSecondDirection( 1 );
  SliceIteratorType inputIt2( Dxy, Dxy->GetRequestedRegion() );
  inputIt2.SetFirstDirection( 0 );
  inputIt2.SetSecondDirection( 1 );
  SliceIteratorType inputIt3( Dxz, Dxz->GetRequestedRegion() );
  inputIt3.SetFirstDirection( 0 );
  inputIt3.SetSecondDirection( 1 );
  SliceIteratorType inputIt4( Dyy, Dyy->GetRequestedRegion() );
  inputIt4.SetFirstDirection( 0 );
  inputIt4.SetSecondDirection( 1 );
  SliceIteratorType inputIt5( Dyz, Dyz->GetRequestedRegion() );
  inputIt5.SetFirstDirection( 0 );
  inputIt5.SetSecondDirection( 1 );
  SliceIteratorType inputIt6( Dzz, Dzz->GetRequestedRegion() );
  inputIt6.SetFirstDirection( 0 );
  inputIt6.SetSecondDirection( 1 );

  ISliceIteratorType inputIt7(I, I->GetRequestedRegion() );
  inputIt7.SetFirstDirection( 0 );
  inputIt7.SetSecondDirection( 1 );

  /*SliceIteratorType outputIt1(Vx1, Vx1->GetRequestedRegion() );
  outputIt1.SetFirstDirection( 2 );
  outputIt1.SetSecondDirection( 1 );
  SliceIteratorType outputIt2(Vy1, Vy1->GetRequestedRegion() );
  outputIt2.SetFirstDirection( 2 );
  outputIt2.SetSecondDirection( 1 );
  SliceIteratorType outputIt3(Vz1, Vz1->GetRequestedRegion() );
  outputIt3.SetFirstDirection( 2 );
  outputIt3.SetSecondDirection( 1 );*/

  GradientSliceIteratorType inputIt0( IGVF, IGVF->GetRequestedRegion() );
  inputIt0.SetFirstDirection( 0 );
  inputIt0.SetSecondDirection( 1 );

  GradientSliceIteratorType outputIt( V1, V1->GetRequestedRegion() );
  outputIt.SetFirstDirection( 0 );
  outputIt.SetSecondDirection( 1 );

  SliceIteratorType outputIt4( IVessel, IVessel->GetRequestedRegion() );
  outputIt4.SetFirstDirection( 0 );
  outputIt4.SetSecondDirection( 1 );

  inputIt0.GoToBegin();
  inputIt1.GoToBegin();
  inputIt2.GoToBegin();
  inputIt3.GoToBegin();
  inputIt4.GoToBegin();
  inputIt5.GoToBegin();
  inputIt6.GoToBegin();
  inputIt7.GoToBegin();

  outputIt.GoToBegin();
  outputIt4.GoToBegin();
  
  double Ma[3][3];
  double V[3][3];
  double d[3];
  
  while( !inputIt1.IsAtEnd() )
  {
   while ( !inputIt1.IsAtEndOfSlice() )
  {
    while ( !inputIt1.IsAtEndOfLine() )
   {
     if( inputIt7.Get() == 0 )
     //if( 0 )
	 {
	   outputIt4.Set(0);
	   GradientPixelType temp;
	   temp[0] = 0;
	   temp[1] = 0;
	   temp[2] = 0;
	   outputIt.Set(temp);
	 }
	 else
	 {
     Ma[0][0] = inputIt1.Get();
     Ma[0][1] = inputIt2.Get();
     Ma[0][2] = inputIt3.Get();
	 Ma[1][0] = inputIt2.Get();
     Ma[1][1] = inputIt4.Get();
     Ma[1][2] = inputIt5.Get();
	 Ma[2][0] = inputIt3.Get();
	 Ma[2][1] = inputIt5.Get();
     Ma[2][2] = inputIt6.Get();
     eigen_decomposition(Ma,V,d);

     double Lambda1 = d[0];
     double Lambda2 = d[1];
	 double Lambda3 = d[2];
	 if(Lambda2>=0.0 || Lambda3>=0.0)
	 {
	   outputIt4.Set(0);
	 }
	 else
	 {
      double Ra  = Lambda2 / Lambda3; 
      double Rb  = Lambda1 / vcl_sqrt ( vnl_math_abs( Lambda2 * Lambda3 )); 
      double S  = vcl_sqrt( pow(Lambda1,2) + pow(Lambda2,2) + pow(Lambda3,2) );
      double vesMeasure_1  = 
         ( 1 - vcl_exp(-1.0*(( vnl_math_sqr( Ra ) ) / ( A ))));
	  double vesMeasure_2  = 
         vcl_exp ( -1.0 * ((vnl_math_sqr( Rb )) /  ( B )));
	  double vesMeasure_3  = 
         ( 1 - vcl_exp( -1.0 * (( vnl_math_sqr( S )) / ( C ))));
	  
	  float V_Saliency = vesMeasure_1 * vesMeasure_2 * vesMeasure_3;
	  //float V_Saliency = fabs(inputIt0.Get()[0] * (float)V[0][0] + inputIt0.Get()[1] * (float)V[0][1] + inputIt0.Get()[2] * (float)V[0][2]);
	  //std::cout<<Lambda1<<","<<Lambda2<<","<<Lambda3<<std::endl;
	  //std::cout<<V_Saliency<<std::endl;
	  outputIt4.Set(V_Saliency);
	 }

	 //outputIt1.Set(V[0][0]);
	 //outputIt2.Set(V[1][0]);
	 //outputIt3.Set(V[2][0]);
	 GradientPixelType temp;
	 temp[0] = (float)V[0][0];
	 temp[1] = (float)V[1][0];
	 temp[2] = (float)V[2][0];
	 outputIt.Set(temp);
	 }
	 ++inputIt0;
     ++inputIt1;
     ++inputIt2;
	 ++inputIt3;
	 ++inputIt4;
     ++inputIt5;
	 ++inputIt6;
	 ++inputIt7;
	 ++outputIt;
	 ++outputIt4;

    }
   inputIt0.NextLine();
   inputIt1.NextLine();
   inputIt2.NextLine();
   inputIt3.NextLine();
   inputIt4.NextLine();
   inputIt5.NextLine();
   inputIt6.NextLine();
   inputIt7.NextLine();
   outputIt.NextLine();
   outputIt4.NextLine();
   }
   inputIt0.NextSlice();
   inputIt1.NextSlice();
   inputIt2.NextSlice();
   inputIt3.NextSlice();
   inputIt4.NextSlice();
   inputIt5.NextSlice();
   inputIt6.NextSlice();
   inputIt7.NextSlice();
   outputIt.NextSlice();
   outputIt4.NextSlice();
  }

  //rescale the vesselness to [0,255]
  typedef itk::RescaleIntensityImageFilter< ProbImageType, ProbImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescale = RescaleFilterType::New();
  rescale->SetInput( IVessel );
  rescale->SetOutputMinimum( 0 );
  rescale->SetOutputMaximum( 255 );
  rescale->Update();
  IVessel = rescale->GetOutput();

}


ProbImagePointer extract_one_component(int index, GradientImagePointer G)
{
    typedef itk::ImageAdaptor<  GradientImageType, 
                              VectorPixelAccessor > ImageAdaptorType;
	ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
	VectorPixelAccessor  accessor;
    accessor.SetIndex( index );
    adaptor->SetPixelAccessor( accessor );
    adaptor->SetImage( G );
    typedef itk::CastImageFilter<ImageAdaptorType,ProbImageType> CasterType;
    CasterType::Pointer caster = CasterType::New();
    caster->SetInput(adaptor);
    caster->Update();
	ProbImagePointer output = caster->GetOutput();

    return output;
}

void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
	int n = 3;
    double e[3];
    double da[3];
    double dt, dat;
    double vet[3];
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            V[i][j] = A[i][j];
        }
    }
    tred2(V, d, e);
    tql2(V, d, e);
    
    /* Sort the eigen values and vectors by abs eigen value */
    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);
    if((da[0]>=da[1])&&(da[0]>da[2]))
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }
    else if((da[1]>=da[0])&&(da[1]>da[2]))  
    {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
        d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2]; 
    }
    if(da[0]>da[1])
    {
        dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
        d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2]; 
    }
}

static void tred2(double V[3][3], double d[3], double e[3]) {
    
	int n = 3;
/*  This is derived from the Algol procedures tred2 by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */
    int i, j, k;
    double scale;
    double f, g, h;
    double hh;
    for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }
    
    /* Householder reduction to tridiagonal form. */
    
    for (i = n-1; i > 0; i--) {
        /* Scale to avoid under/overflow. */
        scale = 0.0;
        h = 0.0;
        for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
        } else {
            
            /* Generate Householder vector. */
            
            for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
            f = d[i-1];
            g = sqrt(h);
            if (f > 0) { g = -g; }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++) { e[j] = 0.0; }
            
            /* Apply similarity transformation to remaining columns. */
            
            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
            hh = f / (h + h);
            for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
            for (j = 0; j < i; j++) {
                f = d[j]; g = e[j];
                for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    
    /* Accumulate transformations. */
    
    for (i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
            for (j = 0; j <= i; j++) {
                g = 0.0;
                for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
            }
        }
        for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
    }
    for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
}

/* Symmetric tridiagonal QL algorithm. */
static void tql2(double V[3][3], double d[3], double e[3]) {
    
	int n = 3;
/*  This is derived from the Algol procedures tql2, by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */
    
    int i, j, k, l, m;
    double f;
    double tst1;
    double eps;
    int iter;
    double g, p, r;
    double dl1, h, c, c2, c3, el1, s, s2;
    
    for (i = 1; i < n; i++) { e[i-1] = e[i]; }
    e[n-1] = 0.0;
    
    f = 0.0;
    tst1 = 0.0;
    eps = pow(2.0, -52.0);
    for (l = 0; l < n; l++) {
        
        /* Find small subdiagonal element */
        
        tst1 = MAX(tst1, fabs(d[l]) + fabs(e[l]));
        m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) { break; }
            m++;
        }
        
        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */
        
        if (m > l) {
            iter = 0;
            do {
                iter = iter + 1;  /* (Could check iteration count here.) */
                /* Compute implicit shift */
                g = d[l];
                p = (d[l+1] - g) / (2.0 * e[l]);
                r = hypot2(p, 1.0);
                if (p < 0) { r = -r; }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                dl1 = d[l+1];
                h = g - d[l];
                for (i = l+2; i < n; i++) { d[i] -= h; }
                f = f + h;
                /* Implicit QL transformation. */
                p = d[m]; c = 1.0; c2 = c; c3 = c;
                el1 = e[l+1]; s = 0.0; s2 = 0.0;
                for (i = m-1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);
                    /* Accumulate transformation. */
                    for (k = 0; k < n; k++) {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                
                /* Check for convergence. */
            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
    
    /* Sort eigenvalues and corresponding vectors. */
    for (i = 0; i < n-1; i++) {
        k = i;
        p = d[i];
        for (j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

double MAX(double a, double b)
{
  return ((a)>(b)?(a):(b));
}

static double hypot2(double x, double y) { return sqrt(x*x+y*y); }

double absd(double val){ if(val>0){ return val;} else { return -val;} };

PGfloat ***pgDmatrix(int M, int N, int K)
{
PGfloat ***Image = NULL;
int i,j;

Image = (PGfloat ***)malloc(M*sizeof(PGfloat **)); /* may need another * in sizeof */
for(i=0; i<M; i++)
{
  Image[i] = (PGfloat **)malloc(N*sizeof(PGfloat*));
  for(j=0; j<N; j++)
  {
    Image[i][j] = (PGfloat *)malloc(K*sizeof(PGfloat));
  }
 }

 return Image;
}

/* Frees a matrix allocated by dmatrix */
PGvoid pgFreeDmatrix(PGfloat ***Image, int M, int N, int K)
{
  int i,j;
   for(i=0; i<M; i++)
   {
     for(j=0; j<N; j++)
     {
       free(Image[i][j]);
     }
     free(Image[i]);
   }
   free(Image);
}

void GVFC(int XN, int YN, int ZN, float *f, float *ou, float *ov, float *oo, 
	  float mu, int ITER)
{
   double mag2, temp, tempx, tempy, tempz, fmax, fmin;
   int count, x, y, z, XN_1, XN_2, YN_1, YN_2, ZN_1, ZN_2;
   
   PGfloat ***fx, ***fy, ***fz, ***u, ***v, ***o, ***Lu, ***Lv, ***Lo, ***g, ***c1, ***c2, ***c3, ***b;
   
   /* define constants and create row-major double arrays */
   XN_1 = XN - 1;
   XN_2 = XN - 2;
   YN_1 = YN - 1;
   YN_2 = YN - 2;
   ZN_1 = ZN - 1;
   ZN_2 = ZN - 2;
   
   fx = pgDmatrix(YN,XN,ZN);
   fy = pgDmatrix(YN,XN,ZN);
   fz = pgDmatrix(YN,XN,ZN);
   u = pgDmatrix(YN,XN,ZN);
   v = pgDmatrix(YN,XN,ZN);
   o = pgDmatrix(YN,XN,ZN);
   Lu = pgDmatrix(YN,XN,ZN);
   Lv = pgDmatrix(YN,XN,ZN);
   Lo = pgDmatrix(YN,XN,ZN);
   g = pgDmatrix(YN,XN,ZN);
   c1 = pgDmatrix(YN,XN,ZN);
   c2 = pgDmatrix(YN,XN,ZN);
   c3 = pgDmatrix(YN,XN,ZN);
   b = pgDmatrix(YN,XN,ZN);
   
   
   /************** I: Normalize the edge map to [0,1] **************/
   fmax = -1e10;
   fmin = 1e10;
   for (x=0; x<=YN*XN*ZN-1; x++) {
      fmax = PG_MAX(fmax,f[x]);
      fmin = PG_MIN(fmin,f[x]);
   } 

   for (x=0; x<=YN*XN*ZN-1; x++) 
      f[x] = (f[x]-fmin)/(fmax-fmin);
  

   /**************** II: Compute edge map gradient *****************/
   /* I.1: Neumann boundary condition: 
    *      zero normal derivative at boundary 
    */
   /* Deal with corners */
  for(z=0; z<= ZN_1; z++)
   for (y=0; y<=YN_1; y++){
     fx[y][0][z] = fy[y][0][z] = fz[y][0][z] = 0;
     fx[y][XN_1][z] = fy[y][XN_1][z] = fz[y][XN_1][z] = 0;
   }
  for(x=0; x<= XN_1; x++)
   for (y=0; y<=YN_1; y++){
     fx[y][x][0] = fy[y][x][0] = fz[y][x][0] = 0;
     fx[y][x][ZN_1] = fy[y][x][ZN_1] = fz[y][x][ZN_1] = 0;
   }
  for(x=0; x<= XN_1; x++)
   for (z=0; z<=ZN_1; z++){
     fx[0][x][z] = fy[0][x][z] = fz[0][x][z] = 0;
     fx[YN_1][x][z] = fy[YN_1][x][z] = fz[YN_1][x][z] = 0;
   }
   
   /* I.2: Compute interior derivative using central difference */
  for(z=1; z<= ZN_2; z++)
   for (y=1; y<=YN_2; y++)
     for (x=1; x<=XN_2; x++) {
	 /* NOTE: f is stored in column major */
	 fx[y][x][z] = 0.5 * (f[y + (x+1)*YN + z*YN*XN] - f[y + (x-1)*YN + z*YN*XN]); 	 
	 fy[y][x][z] = 0.5 * (f[y+1 + x*YN + z*YN*XN] - f[y-1 + x*YN + z*YN*XN]);
     fz[y][x][z] = 0.5 * (f[y + x*YN + (z+1)*YN*XN] - f[y + x*YN + (z-1)*YN*XN]);
     }
   
   /******* III: Compute parameters and initializing arrays **********/
   temp = -1.0/(mu*mu);
  for(z=0; z<=ZN_1; z++)
   for (y=0; y<=YN_1; y++)
      for (x=0; x<=XN_1; x++) {
	 tempx = fx[y][x][z];
	 tempy = fy[y][x][z];
     tempz = fz[y][x][z];
	 /* initial GVF vector */
	 u[y][x][z] = tempx;
	 v[y][x][z] = tempy;
     o[y][x][z] = tempz;
	 /* gradient magnitude square */
	 mag2 = tempx*tempx + tempy*tempy + tempz*tempz; 
	 
	 g[y][x][z] = mu;
	 b[y][x][z] = mag2;

	 c1[y][x][z] = b[y][x][z] * tempx;
	 c2[y][x][z] = b[y][x][z] * tempy;
     c3[y][x][z] = b[y][x][z] * tempz;
      }
  
   /* free memory of fx, fy and fz */
   pgFreeDmatrix(fx,YN,XN,ZN);
   pgFreeDmatrix(fy,YN,XN,ZN);
   pgFreeDmatrix(fz,YN,XN,ZN);
   
   for(z=0; z<= ZN_1; z++)
   for (y=0; y<=YN_1; y++){
     Lu[y][0][z] = Lv[y][0][z] = Lo[y][0][z] = 0;
     Lu[y][XN_1][z] = Lv[y][XN_1][z] = Lo[y][XN_1][z] = 0;
   }
  for(x=0; x<= XN_1; x++)
   for (y=0; y<=YN_1; y++){
     Lu[y][x][0] = Lv[y][x][0] = Lo[y][x][0] = 0;
     Lu[y][x][ZN_1] = Lv[y][x][ZN_1] = Lo[y][x][ZN_1] = 0;
   }
  for(x=0; x<= XN_1; x++)
   for (z=0; z<=ZN_1; z++){
     Lu[0][x][z] = Lv[0][x][z] = Lo[0][x][z] = 0;
     Lu[YN_1][x][z] = Lv[YN_1][x][z] = Lo[YN_1][x][z] = 0;
   }

   /************* Solve GVF = (u,v) iteratively ***************/
   for (count=1; count<=ITER; count++) {
      /* IV: Compute Laplace operator using Neuman condition */
       
      /* IV.4: Compute interior */
   for(z=1; z<=ZN_2; z++)
    for (y=1; y<=YN_2; y++)
	 for (x=1; x<=XN_2; x++) {
	    Lu[y][x][z] = (u[y][x-1][z] + u[y][x+1][z] 
			+ u[y-1][x][z] + u[y+1][x][z] + u[y][x][z-1] + u[y][x][z+1])*1/6 - u[y][x][z];
	    Lv[y][x][z] = (v[y][x-1][z] + v[y][x+1][z] 
			+ v[y-1][x][z] + v[y+1][x][z] + v[y][x][z-1] + v[y][x][z+1])*1/6- v[y][x][z];
	    Lo[y][x][z] = (o[y][x-1][z] + o[y][x+1][z] 
			+ o[y-1][x][z] + o[y+1][x][z] + o[y][x][z-1] + o[y][x][z+1])*1/6 - o[y][x][z];
	 }
      
      /******** V: Update GVF ************/
   for(z=0; z<=ZN_1; z++)
    for (y=0; y<=YN_1; y++)
	 for (x=0; x<=XN_1; x++) {
	    u[y][x][z] = (1- b[y][x][z])*u[y][x][z] + g[y][x][z]*Lu[y][x][z]*4 + c1[y][x][z];
	    v[y][x][z] = (1- b[y][x][z])*v[y][x][z] + g[y][x][z]*Lv[y][x][z]*4 + c2[y][x][z];
        o[y][x][z] = (1- b[y][x][z])*o[y][x][z] + g[y][x][z]*Lo[y][x][z]*4 + c3[y][x][z];
	 }
      
      /* print iteration number */
      printf("%5d",count);
      if (count%15 == 0)
	 printf("\n");
   }
   printf("\n");
   
   /* copy u,v to the output in column major order */
  for(z=0; z<= ZN_1; z++)
   for (y=0; y<=YN_1; y++)
     for (x=0; x<=XN_1; x++) {
	 ou[z*XN*YN+x*YN+y] = u[y][x][z];
	 ov[z*XN*YN+x*YN+y] = v[y][x][z];
     oo[z*XN*YN+x*YN+y] = o[y][x][z];
     }
   
   /* free all the array memory */
   
   pgFreeDmatrix(u,YN,XN,ZN);
   pgFreeDmatrix(v,YN,XN,ZN);
   pgFreeDmatrix(o,YN,XN,ZN);
   pgFreeDmatrix(Lu,YN,XN,ZN);
   pgFreeDmatrix(Lv,YN,XN,ZN);
   pgFreeDmatrix(Lo,YN,XN,ZN);
   pgFreeDmatrix(g,YN,XN,ZN);
   pgFreeDmatrix(c1,YN,XN,ZN);
   pgFreeDmatrix(c2,YN,XN,ZN);
   pgFreeDmatrix(c3,YN,XN,ZN);
   pgFreeDmatrix(b,YN,XN,ZN);
   
   return;
}

float norm_density(float x, float mu, float sigma)
{
   float pi = 3.1415926;
   float p;
   p = 1/MAX(sigma * sqrt(2*pi),std::numeric_limits<float>::epsilon()) * exp( -pow(x-mu,2)/MAX(2*pow(sigma,2),std::numeric_limits<float>::epsilon()));
   return p;
}