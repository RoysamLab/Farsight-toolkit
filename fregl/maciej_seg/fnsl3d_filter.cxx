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

#include "fnsl3d_filter.h"

#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageSeriesReader.h"
#include "itkCastImageFilter.h"
#include "itkObject.h"

//: Constructor
fnsl3d_filter::
fnsl3d_filter()
{
}

//: Destructor
fnsl3d_filter::
~fnsl3d_filter()
{
}

/**	@biref	fitlers image using a median filter and then opens image
*	@param	img	ITK pointer to image
*	@param	radius	radius to use as the kernel for filtering
*	@return	ITK image pointer to filtered image
*/
// Code directly out of Maciej implemenation
fnsl3d_filter::ImageType::Pointer 
fnsl3d_filter::
median_open_image(ImageType::Pointer img,int radius)
{
  //set neighborhood
  ImageType::SizeType miIndexRadius;
  
  miIndexRadius[0] = radius; // radius along x
  miIndexRadius[1] = radius; // radius along y
  miIndexRadius[2] = 0;	// radius along z
  
  MedianFilterType::Pointer f_med = MedianFilterType::New();

  //set radius and input
  f_med->SetRadius( miIndexRadius );
  f_med->SetInput( img );

  OpenFilterType::Pointer f_open = OpenFilterType::New();

  StructuringElementType  structuringElement;
  
  structuringElement.SetRadius( radius );
  
  structuringElement.CreateStructuringElement();
  
  f_open->SetKernel( structuringElement );

  //connect open to medium filter
  f_open->SetInput( f_med->GetOutput() );

  try
    {
      f_open->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in OpenFilter: " << e << std::endl;
      exit(0);
    }


  return f_open->GetOutput();
}

/**	@brief	calculates the max,min,mean and std of the image
*	@param	image	ITK image pointer to input image
*	@param	max	maximum value of image
*	@param	min	minimum value of image
*	@param	mean mean value of image
*	@param	std standard deviation of image
*	@param	std_flag set true if std needs to be calculated
*/
// Code directly out of Maciej implemenation
void 
fnsl3d_filter::
get_min_max_mean_std(ImageType::Pointer image, double &max, double &min, 
                     double &mean, double &std, bool std_flag)
{
  //create iterator for image
  ConstIteratorType cit( image, image->GetRequestedRegion() );
  
  double count=0.0;
  max=0; 
  min=0; 
  bool init = true;
  mean=0.0;
  std=0.0;
  
  //iterate throught image and finds min max
  for ( cit.GoToBegin(); !cit.IsAtEnd(); ++cit )
    {
      double i=cit.Get();
      if(init)
        {
          max=i;
          min=i;
          init = false;
        }
      if(i>max)
        max=i;
      if(i<min)
        min=i;
      mean+=i;
      count++;
    }
  mean/=count;

  //calculate std
  if(std_flag)
    {
      for ( cit.GoToBegin(); !cit.IsAtEnd(); ++cit )
        {
          double i=cit.Get();
          std+=(i-mean)*(i-mean);
        }
      std/=count;
      std=sqrt(std);
    }
}

/**	@brief	filters the image using a morphological gradient filter
*	uses an open and erode filter to get morphological gradient
*	@param	img ITK image pointer to input image
*	@param	radius	size of filtering kernel
*	@return	ITK pointer to the morphological gradient of image
*/
fnsl3d_filter::ImageType::Pointer 
fnsl3d_filter::
morph_gradient(ImageType::Pointer img,int radius)
{
  //create filters
  ErodeFilterType::Pointer  f_erode  = ErodeFilterType::New();
  DilateFilterType::Pointer f_dilate = DilateFilterType::New();
  OpenFilterType::Pointer f_open = OpenFilterType::New();
  SubFilterType::Pointer f_sub = SubFilterType::New();

  StructuringElementType  structuringElement;
  structuringElement.SetRadius( radius );
  structuringElement.CreateStructuringElement();
  
  //set kernel 
  f_erode->SetKernel(  structuringElement );
  f_dilate->SetKernel( structuringElement );
  f_open->SetKernel( structuringElement );
  
  //connect open to medium filter
  f_open->SetInput( img );
  
  //connect open to erode and dilate
  f_erode->SetInput(  f_open->GetOutput() );
  f_dilate->SetInput( f_open->GetOutput() ); 
  
  //set inputs of substract filter for morphological gradient
  //(dilation-erosion=morphological gradient)
  f_sub->SetInput1( f_dilate->GetOutput() );
  f_sub->SetInput2( f_erode->GetOutput() );
  
  try
    {
      f_sub->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in SubFilter: " << e << std::endl;
      exit(0);
    }
  return f_sub->GetOutput();
}

/**	@brief	thresholds the image and then calculates the distance map
*	@param	img	ITK image pointer to input image
*	@param	hor	number of horizontal windows for thresholding
*	@param	ver	number of vertical windows for thresholding
*	@param	dep	number of depth windows for thresholding
*	@param	param	parameter used in calculating the threshold value
*	@return	ITK image pointer to filtered image
*/
fnsl3d_filter::ImageType::Pointer 
fnsl3d_filter::
threshold_distance_map(ImageType::Pointer img,int hor,int ver,int dep,
                       double param)
{
  //create temp image
  ImageType::Pointer im_thresh = ImageType::New();
  im_thresh->SetRegions(region_);
  im_thresh->CopyInformation( img );
  im_thresh->Allocate();
  
  int dsize=0;
  for(int k=0;k<dep;k++) {
    int vsize=0;
    for(int j=0;j<ver;j++) {
      int hsize=0;
      for(int i=0;i<hor;i++) {
        //divide up the regions into a grid for processing
        ImageType::IndexType start;
        start[0] = hsize;
        start[1] = vsize;
        start[2] = dsize;
        
        ImageType::SizeType size;
        size[0] =(image_size_[0] - hsize)/(hor - i);
        size[1] =(image_size_[1] - vsize)/(ver - j);
        size[2] =(image_size_[2] - dsize)/(dep - k);
        
        //set the desired region
        ImageType::RegionType desiredRegion;
        desiredRegion.SetSize(  size  );
        desiredRegion.SetIndex( start );
        
        //create region of interest filter
        RegionFilterType::Pointer f_region = RegionFilterType::New();

        //set the region of the filter and update filter
        f_region->SetRegionOfInterest( desiredRegion );
        f_region->SetInput( img );
        try
          {
            f_region->Update();
          }
        catch (itk::ExceptionObject & e)
          {
            std::cerr << "Exception in RegionFilter: " << e << std::endl;
            exit(0);
          }
        
        
        //create binary threshold filter
        ThreshFilterType::Pointer f_thresh = ThreshFilterType::New();

        //set the binary output
        const PixelType outsideValue = 0;
        const PixelType insideValue  = 255;
        
        f_thresh->SetOutsideValue( outsideValue );
        f_thresh->SetInsideValue(  insideValue  );
        
        double max, min, mean,std;
        //get the min max mean and std
        get_min_max_mean_std(f_region->GetOutput(),max,min,mean,std,true);
        
        //set the threshold
        PixelType lowerThreshold = 0;
        PixelType upperThreshold = mean + param*std;
        
        if( upperThreshold < 1 )
          upperThreshold = 1;
        else if( upperThreshold > 254 )
          upperThreshold = 254;
        
        f_thresh->SetLowerThreshold( lowerThreshold );
        f_thresh->SetUpperThreshold( upperThreshold );
        
        //update
        f_thresh->SetInput(f_region->GetOutput());
        
        try
          {
            f_thresh->Update();
          }
        catch (itk::ExceptionObject & e)
          {
            std::cerr << "Exception in ThreshFilter: " << e << std::endl;
            exit(0);
          }
        
        //create iterators
        IteratorType it( im_thresh,desiredRegion );
        ConstIndexIteratorType ciit( f_thresh->GetOutput(),
                                     f_thresh->GetOutput()->GetLargestPossibleRegion() );
        
        //save into image
        for ( ciit.GoToBegin(),it.GoToBegin(); !ciit.IsAtEnd(); 
              ++ciit, ++it)
          {
            it.Set(ciit.Get());
            //save miIndex if pixel is background
            if(ciit.Get() == 255)
              vIndex_.push_back(it.GetIndex());
          }
        hsize+=image_size_[0]/hor;
      }
      vsize+=image_size_[1]/ver;
    }
    dsize+=image_size_[2]/dep;
  }
  
  /*
  //save_image( im_thresh, "thred_image.tif");
  typedef itk::Image< unsigned char,  3 >   InputImageType;
  typedef itk::Image< float,  3 >   tempImageType;
  typedef itk::DanielssonDistanceMapImageFilter<
               InputImageType, tempImageType >  F_Type;
  typedef itk::ImageSeriesReader< InputImageType  > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->AddFileName( "thred_image.tif" );
  reader->AddFileName( "thred_image.tif" );
  reader->AddFileName( "thred_image.tif" );
  F_Type::Pointer f_daniel =  F_Type::New();
  f_daniel->InputIsBinaryOn();
  f_daniel->SetInput(reader->GetOutput());
  */
  /*
  typedef itk::Image< unsigned char,  3 >   InputImageType;
  typedef itk::CastImageFilter< ImageType, ImageType> CastFilterType;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput( im_thresh );
  typedef itk::DanielssonDistanceMapImageFilter<
               ImageType, ImageType >  F_Type;
  F_Type::Pointer f_daniel =  F_Type::New();
  f_daniel->InputIsBinaryOn();
  f_daniel->SetInput(castFilter->GetOutput());
  */
  
  // Please note: DanielFilter cannot handle image volume of one slice.
  DanielFilterType::Pointer f_daniel =  DanielFilterType::New();

  //update
  f_daniel->InputIsBinaryOn();
  f_daniel->SetInput(im_thresh);
  

  try
    {
      f_daniel->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in DanielFiltere: " << e << std::endl;
      exit(0);
    }
  return f_daniel->GetOutput();
}

/**	@brief	inverts the distance image and then combines with morphological gradient
*	@param	imgDaniel ITK image pointer to distance image
*	@param	imgMorphGrad ITK image pointer to morphological gradient image
*	@param	sigma	sigma for gaussian filter
*	@return	ITK image pointer to filtered image
*/
fnsl3d_filter::ImageType::Pointer 
fnsl3d_filter::
combine_invert_distance(ImageType::Pointer imgDaniel, 
                        ImageType::Pointer imgMorphGrad,
                        double sigma)
{
  double max, min, mean,std;
  //calculate the max and min
  get_min_max_mean_std(imgMorphGrad,max,min,mean,std,false);
  
  IteratorType it_dan(imgDaniel,imgDaniel->GetLargestPossibleRegion());
  ConstIteratorType cit_morph(imgMorphGrad,imgMorphGrad->GetLargestPossibleRegion());
  
  int i;
  //calculate normalized gradient and save the gradient into temp image
  for (i=0, it_dan.GoToBegin(),cit_morph.GoToBegin(); !cit_morph.IsAtEnd(); 
       ++it_dan, ++cit_morph,i++)
    {
      it_dan.Set(it_dan.Get() * exp(  (max-cit_morph.Get() ) / (max-min)  ));
    }

  GaussFilterType::Pointer f_gauss = GaussFilterType::New();

  //set the variance and update
  f_gauss->SetVariance(sigma);
  f_gauss->SetInput(imgDaniel);
  
  try
    {
      f_gauss->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in GaussFilter: " << e << std::endl;
      exit(0);
    }
  
  
  //get max to invert image
  get_min_max_mean_std(f_gauss->GetOutput(),max,min,mean,std,false);
  
  //invert distnace
  IteratorType it_gauss(f_gauss->GetOutput(),f_gauss->GetOutput()->GetLargestPossibleRegion());
  
  for(it_gauss.GoToBegin();!it_gauss.IsAtEnd();++it_gauss)
    {
      it_gauss.Set( max - it_gauss.Get() );
    }

  IndexIteratorType iit_gauss( f_gauss->GetOutput(), f_gauss->GetOutput()->GetLargestPossibleRegion() );
  //force background to be 255
  iit_gauss.GoToBegin();
  for (unsigned int i=0;i<vIndex_.size();i++)
    {
      iit_gauss.SetIndex(vIndex_[i]);
      iit_gauss.Set(255);
    }
 
  return f_gauss->GetOutput();
}

/**	@brief	prefilter the images and save all neccessary temp images
*	@param  img grayscale image 
*/
void
fnsl3d_filter::
run_filter(InputImageType::Pointer img, int rad_med, int rad_grad, 
           int gridX, int gridY, double paramT, double sigma )
{
  /*
  //load green channel
  mOutputImageType::Pointer point = load_channel(crsName,crsPath,crsType,
                                                 crStart,crEnd,crWidth);
  */

  region_ = img->GetLargestPossibleRegion();
  ImageType::SizeType im_size = region_.GetSize();
  image_size_[0] = im_size[0];
  image_size_[1] = im_size[1];
  image_size_[2] = im_size[2];
  
  typedef itk::CastImageFilter<fnsl3d_filter::InputImageType, 
    ImageType> castFilterType;
  castFilterType::Pointer caster = castFilterType::New();
  caster->SetInput(img);
  imgpIntenisty_ = caster->GetOutput();
 
  std::cout << "Median Filter" << std::endl;
  ImageType::Pointer im_thresh = median_open_image(imgpIntenisty_, rad_med);

  //create temp image
  /*
  ImageType::Pointer im_thresh = ImageType::New();
  im_thresh->SetRegions(region_);
  im_thresh->CopyInformation( point );
  im_thresh->Allocate();

  ConstIteratorType cit( point,point->GetLargestPossibleRegion() );

  IteratorType it( im_thresh,im_thresh->GetLargestPossibleRegion() );

  //save into temp image
  for ( cit.GoToBegin(),it.GoToBegin(); !cit.IsAtEnd(); ++cit,++it)
    {
      it.Set( cit.Get() );
    }
  */
  std::cout << "Threshold" << std::endl;
  ImageType::Pointer im_thrd_dst_map = 
    threshold_distance_map(im_thresh, gridX, gridY, 1 , paramT);
 
  std::cout << "Morphing" << std::endl;
  imgpGradient_ = morph_gradient(im_thresh, rad_grad);

  std::cout << "Inverting" << std::endl;
  imgpFilt_ = combine_invert_distance(im_thrd_dst_map, imgpGradient_ , sigma);

  // temp code for debugging
  /*
  save_image( imgpIntenisty_, "intensity_image.tif");
  save_image( im_thrd_dst_map, "thrd_dst_map.tif");
  save_image( imgpGradient_, "Gradient_image.tif");
  save_image( imgpFilt_, "Filt_image.tif");
  */
}

/**	@brief	saves the image to a file
*	@param img	ITK image pointer to input image
*	@param	file	file name of image
*/
void
fnsl3d_filter::
save_image(ImageType::Pointer img, std::string file)
{
  typedef float                              PixelType;
  typedef itk::Image<PixelType, 3>           ImageType;
  typedef unsigned char                      WritePixelType;
  typedef itk::Image< WritePixelType, 3 >    WriteImageType;
  typedef itk::RescaleIntensityImageFilter< ImageType, 
                                            WriteImageType > RescalerFilterType;
  typedef itk::ImageFileWriter< WriteImageType >  WriterType;

  RescalerFilterType::Pointer f_scaler_file = RescalerFilterType::New();
  
  
  WriterType::Pointer writer = WriterType::New();
  
  writer->SetFileName( file.c_str() );
  
  f_scaler_file->SetOutputMaximum( 255 );
  f_scaler_file->SetOutputMinimum(  0 );
  
  f_scaler_file->SetInput(img);
  
  writer->SetInput( f_scaler_file->GetOutput() );
  
  try
    {
      writer->Update();
    }
  catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Writer: " << e << std::endl;
      exit(0);
    }
  
}

int*
fnsl3d_filter::
get_size()
{
  return image_size_;
}
