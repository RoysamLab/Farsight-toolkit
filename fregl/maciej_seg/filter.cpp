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

/** @file filter.cpp
*   @brief class for filtering
*   This is the class that filter the image prior to merging
*
*   @author Maciej Wotjon
*/

//#include "FarSight.h"

#include "filter.h"
#include "itkVector.h"
#include "itkListSample.h"
#include "itkMatrix.h"
#include "itkMeanCalculator.h"
#include "itkCovarianceCalculator.h"

#include <iostream> // for cout
#include <sstream> // for ostringstream
#include <iomanip> // for setw
#include <algorithm> // for find_if
#include <vcl_fstream.h>
#include <fstream>

#ifndef M_PI
#  define M_PI 3.14159265358979323846f
#endif

/**	@brief ints filter class
*	@param gX size of grid in x drirection for thresholding
*	@param gY size of grid in y drirection for thresholding
*	@param rM radius of median filter
*	@param rG radius of morph opertations
*	@param pT parameter for thresholding
*	@param s parameter for gauss smooth
*/
filter::filter(int gX, int gY, int rM, int rG, double pT, double s,int channel)
: m_gridX(gX)
, m_gridY(gY)
, m_radMed(rM)
, m_radGrad(rG)
, m_paramT(pT)
, m_sigma(s)
, mChannel(channel)
{
	m_nImage_size[0] = m_nImage_size[1] = m_nImage_size[2] = 0;
}


filter::~filter(void)
{
	for(unsigned int i=0;i<mvPix.size();++i )	
	{
		delete mvPix[i];
	}
}

/**	@brief load the green channel from a series of images and save the rgp into a temp image
*	@return pointer ITK image
*/
filter::mOutputImageType::Pointer
filter::load_green_channel(const std::string &crsName,const std::string &crsPath,
						   const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth)
{

	//create reader
	mReaderType::Pointer reader = mReaderType::New();

	reader->SetFileNames( gen_names(crsPath.c_str(),
		crsName.c_str(),crsType.c_str(), crStart,crEnd,crWidth));

	//create adaptor
	mImageGreenAdaptorType::Pointer f_adaptor = mImageGreenAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerGreenAdaptorType::Pointer f_rescaler_adap = mRescalerGreenAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}
filter::mOutputImageType::Pointer
filter::load_green_channel(const std::string &filename)
{

	//create reader
	mReaderFileType::Pointer reader = mReaderFileType::New();

	reader->SetFileName( filename);

	//create adaptor
	mImageGreenAdaptorType::Pointer f_adaptor = mImageGreenAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerGreenAdaptorType::Pointer f_rescaler_adap = mRescalerGreenAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}


/**	@brief generates string for input file names
*	@param path path to files
*	@param file_first name of first file
*	the file needs to be a tiff and name scheme must have ascending numbers at end
*	@return	vector with generated string file names
*/
std::vector<std::string>
filter::gen_names(const std::string &path, const std::string &file_first,const std::string &type,int start,int end,int width)
{
	std::string temp,n;
	std::vector<std::string> output;
	for(int i=start;i<=end;i++)
	{
		std::ostringstream oss;
		oss << std::setw(width) << std::setfill('0') << i;

		temp = path + file_first + oss.str() + '.' + type;
		std::cerr << temp << std::endl;
		output.push_back(temp);
	}
	return output;
}


/**	@brief save intensity of the image and init pts and save image size
*	@param img ITK image pointer
*/
void 
filter::save_intensity(mOutputImageType::Pointer img)
{
	//get the whole image region
	mRegion = img->GetLargestPossibleRegion();
	//save the image size
	mOutputImageType::SizeType im_size = mRegion.GetSize();
	m_nImage_size[0] = im_size[0];
	m_nImage_size[1] = im_size[1];
	m_nImage_size[2] = im_size[2];


	//save intensity
	mConstIteratorType cit(img,img->GetLargestPossibleRegion());
	mImgpIntenisty = mOutputImageType::New();
	mImgpIntenisty->SetRegions( img->GetLargestPossibleRegion() );
	mImgpIntenisty->CopyInformation( img );
	mImgpIntenisty->Allocate();
	itk::ImageRegionIterator< mOutputImageType > it(mImgpIntenisty,mImgpIntenisty->GetLargestPossibleRegion());;
	for(it.GoToBegin(),cit.GoToBegin();!it.IsAtEnd();++it,++cit)
	{
		it.Set(cit.Get());
	}
}

/**	@biref	fitlers image using a median filter and then opens image
*	@param	img	ITK pointer to image
*	@param	radius	radius to use as the kernel for filtering
*	@return	ITK image pointer to filtered image
*/
filter::mOutputImageType::Pointer
filter::median_open_image(mOutputImageType::Pointer img, int radius)
{

	//set neighborhood
	mOutputImageType::SizeType miIndexRadius;

	miIndexRadius[0] = radius; // radius along x
	miIndexRadius[1] = radius; // radius along y
	miIndexRadius[2] = 0;	// radius aling z

	mMedianFilterType::Pointer f_med = mMedianFilterType::New();

	//set radius and input
	f_med->SetRadius( miIndexRadius );
	f_med->SetInput( img );

	mOpenFilterType::Pointer f_open = mOpenFilterType::New();

	mStructuringElementType  structuringElement;

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
void 
filter::get_min_max_mean_std(mOutputImageType::Pointer image, double &max, double &min, double &mean, double &std,bool std_flag){

	//create iterator for image
	mConstIteratorType cit( image, image->GetRequestedRegion() );

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
filter::mOutputImageType::Pointer 
filter::morph_gradient(mOutputImageType::Pointer img,int radius)
{
	//create filters
	mErodeFilterType::Pointer  f_erode  = mErodeFilterType::New();
	mDilateFilterType::Pointer f_dilate = mDilateFilterType::New();
	mOpenFilterType::Pointer f_open = mOpenFilterType::New();
	mSubFilterType::Pointer f_sub = mSubFilterType::New();

	mStructuringElementType  structuringElement;
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

	//set inputs of substract filter for morphological gradient (dilation-erosion=morphological gradient)
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

/**	@brief	saves the gradient
*	@param	img	ITK image pointer to input image
*/
void
filter::save_gradient(mOutputImageType::Pointer img)
{
	//save intensity
	mConstIteratorType cit(img,img->GetLargestPossibleRegion());
	mImgpGradient = mOutputImageType::New();
	mImgpGradient->SetRegions( img->GetLargestPossibleRegion() );
	mImgpGradient->CopyInformation( img );
	mImgpGradient->Allocate();
	itk::ImageRegionIterator< mOutputImageType > it(mImgpGradient,mImgpGradient->GetLargestPossibleRegion());;
	for(it.GoToBegin(),cit.GoToBegin();!it.IsAtEnd();++it,++cit)
	{
		it.Set(cit.Get());
	}
}

/**	@brief	thresholds the image and then calculates the distance map
*	@param	img	ITK image pointer to input image
*	@param	hor	number of horizontal windows for thresholding
*	@param	ver	number of vertical windows for thresholding
*	@param	dep	number of depth windows for thresholding
*	@param	param	parameter used in calculating the threshold value
*	@return	ITK image pointer to filtered image
*/
filter::mOutputImageType::Pointer 
filter::threshold_distance_map(mOutputImageType::Pointer img,int hor,int ver,int dep,double param)
{
	//create temp image
	mOutputImageType::Pointer im_thresh = mOutputImageType::New();
	im_thresh->SetRegions(mRegion);
	im_thresh->Allocate();

	int dsize=0;
	for(int k=0;k<dep;k++)
	{
		int vsize=0;
		for(int j=0;j<ver;j++)
		{
			int hsize=0;
			for(int i=0;i<hor;i++)
			{
				//divide up the regions into a grid for processing
				mOutputImageType::IndexType start;
				start[0] = hsize;
				start[1] = vsize;
				start[2] = dsize;

				mOutputImageType::SizeType size;
				size[0] =(m_nImage_size[0] - hsize)/(hor - i);
				size[1] =(m_nImage_size[1] - vsize)/(ver - j);
				size[2] =(m_nImage_size[2] - dsize)/(dep - k);

				//set the desired region
				mOutputImageType::RegionType desiredRegion;
				desiredRegion.SetSize(  size  );
				desiredRegion.SetIndex( start );

				//create region of interest filter
				mRegionFilterType::Pointer f_region = mRegionFilterType::New();

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
				mThreshFilterType::Pointer f_thresh = mThreshFilterType::New();

				//set the binary output
				const mOutputPixelType outsideValue = 0;
				const mOutputPixelType insideValue  = 255;

				f_thresh->SetOutsideValue( outsideValue );
				f_thresh->SetInsideValue(  insideValue  );

				double max, min, mean,std;
				//get the min max mean and std
				get_min_max_mean_std(f_region->GetOutput(),max,min,mean,std,true);

				//set the threshold
				mOutputPixelType lowerThreshold = 0;
				mOutputPixelType upperThreshold = mean + param*std;

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
				mIteratorType it( im_thresh,desiredRegion );
				mConstIndexIteratorType ciit( f_thresh->GetOutput(),f_thresh->GetOutput()->GetLargestPossibleRegion() );

				//save into image
				for ( ciit.GoToBegin(),it.GoToBegin(); !ciit.IsAtEnd(); 
					++ciit, ++it)
				{
					it.Set(ciit.Get());
					//save miIndex if pixel is background
					if(ciit.Get() == 255)
						mvIndex.push_back(it.GetIndex());
				}
				hsize+=m_nImage_size[0]/hor;
			}
			vsize+=m_nImage_size[1]/ver;
		}
		dsize+=m_nImage_size[2]/dep;
	}
	mDanielFilterType::Pointer f_daniel =  mDanielFilterType::New();

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

filter::mOutputImageType::Pointer 
filter::combine_invert_distance(mOutputImageType::Pointer imgDaniel,mOutputImageType::Pointer imgMorphGrad,double sigma)
{
	double max, min, mean,std;
	//calculate the max and min
	get_min_max_mean_std(imgMorphGrad,max,min,mean,std,false);

	mIteratorType it_dan(imgDaniel,imgDaniel->GetLargestPossibleRegion());
	mConstIteratorType cit_morph(imgMorphGrad,imgMorphGrad->GetLargestPossibleRegion());

	int i;
	//calculate normalized gradient and save the gradient into temp image
	for (i=0, it_dan.GoToBegin(),cit_morph.GoToBegin(); !cit_morph.IsAtEnd(); 
		++it_dan, ++cit_morph,i++)
	{
		it_dan.Set(it_dan.Get() * exp(  (max-cit_morph.Get() ) / (max-min)  ));
	}

	mGaussFilterType::Pointer f_gauss = mGaussFilterType::New();

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
	mIteratorType it_gauss(f_gauss->GetOutput(),f_gauss->GetOutput()->GetLargestPossibleRegion());

	for(it_gauss.GoToBegin();!it_gauss.IsAtEnd();++it_gauss)
	{
		it_gauss.Set( max - it_gauss.Get() );
	}

	mIndexIteratorType iit_gauss( f_gauss->GetOutput(), f_gauss->GetOutput()->GetLargestPossibleRegion() );
	//force background to be 255
	iit_gauss.GoToBegin();
	for (unsigned int i=0;i<mvIndex.size();i++)
	{
		iit_gauss.SetIndex(mvIndex[i]);
		iit_gauss.Set(255);
	}
	return f_gauss->GetOutput();
}

/**	@brief	saves the rgp image to a file
*	@param	file	file name of image
*/
filter::mInputImageFileType::Pointer
filter::save_image_rgb(std::string file)
{
	typedef itk::ImageFileWriter< mInputImageFileType >  tWriterType;
	tWriterType::Pointer writer = tWriterType::New();

	writer->SetFileName( file.c_str() );

	writer->SetInput( mImgpRGBimage );
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in Writer: " << e << std::endl;
		exit(0);
	}
        return mImgpRGBimage;
}

/**	@brief	saves the image to a file
*	@param img	ITK image pointer to input image
*	@param	file	file name of image
*/
void
filter::save_image(mOutputImageType::Pointer img,std::string file)
{

	mRescalerFilterType::Pointer f_scaler_file = mRescalerFilterType::New();


	mWriterType::Pointer writer = mWriterType::New();

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

/**	@brief	prefilter the images and save all neccessary temp images
*	@param	crsName name of file name of the image stacks, only the name without the numbering
*	@param	crsPath path to file
*	@param	crsType type of image, ei extension tif
*	@param	crStart starting number of images ie image_nameXX.tif XX->sstack
*	@param	crEnd ending stack number
*	@param	crWidth width of number ie image_nameXX.tif swidth -> 2 , image_nameXXX.tif swidth -> 3
*/
void
filter::run_filter(const std::string &crsName,const std::string &crsPath,const std::string &crsType,const int &crStart,const int &crEnd, const int &crWidth)
{
	//load green channel
	mOutputImageType::Pointer point = load_channel(crsName,crsPath,crsType,crStart,crEnd,
				   crWidth);

	save_intensity(point);

	std::cout << "Median Filter" << std::endl;
	point = median_open_image(point,m_radMed);

	//create temp image
	mOutputImageType::Pointer im_thresh = mOutputImageType::New();
	im_thresh->SetRegions(mRegion);
	im_thresh->CopyInformation( point );
	im_thresh->Allocate();

	mConstIteratorType cit( point,point->GetLargestPossibleRegion() );

	mIteratorType it( im_thresh,im_thresh->GetLargestPossibleRegion() );

	//save into temp image
	for ( cit.GoToBegin(),it.GoToBegin(); !cit.IsAtEnd(); ++cit,++it)
	{
		it.Set( cit.Get() );
	}

	std::cout << "Threshold" << std::endl;
	point = threshold_distance_map(point,m_gridX,m_gridY,1,m_paramT);

	std::cout << "Morphing" << std::endl;
	mOutputImageType::Pointer point2 = morph_gradient(im_thresh,m_radGrad);

	save_gradient(point2);

	std::cout << "Inverting" << std::endl;
	point = combine_invert_distance(point,point2,m_sigma);

	mImgpFilt = point;
}

void
filter::run_filter(const std::string& filename)
{
	//load green channel
	mOutputImageType::Pointer point = load_channel(filename);

	save_intensity(point);

        mOutputImageType::IndexType index;
        index[0] = 100;
        index[1] = 100;
        index[2] = 3;
        point->GetPixel(index);

	std::cout << "Median Filter" << std::endl;
	point = median_open_image(point,m_radMed);

	//create temp image
	mOutputImageType::Pointer im_thresh = mOutputImageType::New();
	im_thresh->SetRegions(mRegion);
	im_thresh->CopyInformation( point );
	im_thresh->Allocate();
	mConstIteratorType cit( point,point->GetLargestPossibleRegion() );

	mIteratorType it( im_thresh,im_thresh->GetLargestPossibleRegion() );

	//save into temp image
	for ( cit.GoToBegin(),it.GoToBegin(); !cit.IsAtEnd(); ++cit,++it)
	{
		it.Set( cit.Get() );
	}
 
	std::cout << "Threshold" << std::endl;
	point = threshold_distance_map(point,m_gridX,m_gridY,1,m_paramT);

	std::cout << "Morphing" << std::endl;
	mOutputImageType::Pointer point2 = morph_gradient(im_thresh,m_radGrad);

	save_gradient(point2);

	std::cout << "Inverting" << std::endl;
	point = combine_invert_distance(point,point2,m_sigma);

	mImgpFilt = point;
}

/**	@brief label image with cell labels
*	labels each cell with the number at the cell center
*/
void
filter::label_image(std::vector<cell*> cells)
{
	mInputGreenPixelFileType pt;
	pt.SetRed(255);
	pt.SetGreen(255);
	pt.SetBlue(255);

	mInputImageFileType::IndexType ind;
	//loop through all cells
	for(unsigned int i=0;i<cells.size();i++)
	{
		if(cells[i] == NULL)
			continue;
		std::string strmiIndex;
		std::ostringstream oss;
		oss << cells[i]->mLabel;
		strmiIndex = oss.str();
		
		//initialize the miIndex, make sure to update the center before calling
		mInputImageFileType::IndexType pixIndex;
		pixIndex[0] = cells[i]->mCenterX;
		pixIndex[1] = cells[i]->mCenterY;
		pixIndex[2] = cells[i]->mCenterZ;
		ind[2] = cells[i]->mCenterZ;

		//save the image size
		mInputImageFileType::SizeType size = mImgpRGBimage->GetLargestPossibleRegion().GetSize();

		//spacer to move numbers over
		// int nx_spacer=0;
		for(unsigned int j=0;j<strmiIndex.size();j++)
		{
			if( pixIndex[0] > size[0] + 3 || pixIndex[1] > size[1] + 5)
				break;

			//make pixel label for each number
			switch(strmiIndex[j])
			{
			case '1':
				//1
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '2':
				//1
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '3':
				//1
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '4':
				//1
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '5':
				//1
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//10
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//11
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//12
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '6':
				//1
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//10
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				//11
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//12
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '7':
				//1
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '8':
				//1
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '9':
				//1
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			case '0':
				//1
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 0;
				mImgpRGBimage->SetPixel(ind,pt);

				//2
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//3
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 1;
				mImgpRGBimage->SetPixel(ind,pt);

				//4
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//5
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 2;
				mImgpRGBimage->SetPixel(ind,pt);

				//6
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//7
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 3;
				mImgpRGBimage->SetPixel(ind,pt);

				//8
				ind[0] = pixIndex[0] + 0;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//9
				ind[0] = pixIndex[0] + 2;
				ind[1] = pixIndex[1] + 4;
				mImgpRGBimage->SetPixel(ind,pt);

				//10
				ind[0] = pixIndex[0] + 1;
				ind[1] = pixIndex[1] + 5;
				mImgpRGBimage->SetPixel(ind,pt);

				break;
			default:
				break;
			}
			//advance spacer
			pixIndex[0] += 4;
		}
	}
}
/**	@brief label image with cell boundaries
*	@param	cells vector of cells to be outlined
*/
void
filter::label_bound(std::vector<cell*> cells)
{
	//set color
	mInputGreenPixelFileType pt;
	pt.SetRed(255);
	pt.SetGreen(0);
	pt.SetBlue(0);

	mInputGreenPixelFileType ptG;
	mInputGreenPixelFileType ptN;
	if( mChannel == 0 )
	{
		std::cout << "Labeling Red Channel" << std::endl;
		ptG.SetRed(0);
		ptG.SetGreen(255);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(0);
		ptN.SetBlue(255);
	}
	else if( mChannel == 1 )
	{
		std::cout << "Labeling Green Channel" << std::endl;
		ptG.SetRed(255);
		ptG.SetGreen(0);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(0);
		ptN.SetBlue(255);
	}
	else if( mChannel == 2 )
	{
		std::cout << "Labeling Blue Channel" << std::endl;
		ptG.SetRed(255);
		ptG.SetGreen(0);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(255);
		ptN.SetBlue(0);
	}
	else
	{
		ptG.SetRed(255);
		ptG.SetGreen(255);
		ptG.SetBlue(255);

		ptN.SetRed(255);
		ptN.SetGreen(255);
		ptN.SetBlue(255);
	}


	mInputImageFileType::IndexType ind;
	//loop through all cells
	for(unsigned int i=0;i<cells.size();i++)
	{
		if(cells[i] == NULL)
			continue;
		//if( !cells[i]->mIsTrain )
			//continue;
		if( cells[i]->mClass == 2 )//&& cells[i]->mIsTrain)
		{
			pt = ptG;
		}
		else if( cells[i]->mClass == 1 )// && cells[i]->mIsTrain )
		{
			pt = ptN;
		}
		else
		{
			pt.SetRed(255);
			pt.SetGreen(0);
			pt.SetBlue(255);
		}

		for(unsigned int j=0;j<cells[i]->mvBoundPix2D.size();j++)
		{
			ind[0] = cells[i]->mvBoundPix2D[j]->x_;
			ind[1] = cells[i]->mvBoundPix2D[j]->y_;
			ind[2] = cells[i]->mvBoundPix2D[j]->z_;
			mImgpRGBimage->SetPixel(ind,pt);
		}
	}
}

/**	@brief get the intensity image
*	@return ITK pointer to image
*/
filter::mOutputImageType::Pointer
filter::getintensity(void)
{
	return mImgpIntenisty;
}

/**	@brief get the filtered image
*	@return ITK pointer to image
*/
filter::mOutputImageType::Pointer
filter::getfilt(void)
{
	return mImgpFilt;
}

/**	@brief get the gradient image
*	@return ITK pointer to image
*/
filter::mOutputImageType::Pointer
filter::getgradient(void)
{
	return mImgpGradient;
}

/** @breif get the image size
*	@return return array with cell size
*/
int*
filter::getsize(void)
{
	return m_nImage_size;
}

/**	@brief	load filtered images
*	@param	sFilt name of filtered image to load, can either be the name with extension or full path to file
*	@param  sGrad name of gradient image to load, can either be the name with extension or full path to file
*	@param  sInt name of intensity image to load, can either be the name with extension or full path to file
*	@param  sPath	path ot RGB image
*	@param  sRGB	name of begining of RGB file, the end is the numbered images which is generated
*	@param  sType	extension of RGB file
*	@param  start	staring index for RGB
*	@param  end	ending index for RGB
*	@param  width	width of numbered string at end of file name
*/
void
filter::load_images(std::string sFilt, std::string sGrad, std::string sInt, std::string sPath, std::string sRGB,
					std::string sType,const int &start,const int &end,const int &width,int &loadimages)
{
	typedef itk::ImageFileReader< mOutputImageType > ReaderType;

	ReaderType::Pointer reader1 = ReaderType::New();
	ReaderType::Pointer reader2 = ReaderType::New();
	ReaderType::Pointer reader3 = ReaderType::New();
	mReaderType::Pointer reader4 = mReaderType::New();

	//load filtered image
	reader1->SetFileName(sFilt.c_str());

	try
	{
		reader1->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		loadimages = 0;
		return;
	}

	mImgpFilt = mOutputImageType::New();
	mImgpFilt = reader1->GetOutput();

	//save gradient image
	reader2->SetFileName(sGrad.c_str());

	try
	{
		reader2->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		loadimages = 0;
		return;
	}
	save_gradient(reader2->GetOutput());

	//save intensity image
	reader3->SetFileName(sInt.c_str());

	try
	{
		reader3->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		loadimages = 0;
		return;
	}
	save_intensity(reader3->GetOutput());

	//load rgp image
	reader4->SetFileNames( gen_names(sPath.c_str(),
		sRGB.c_str(),sType.c_str(), start , end , width));
	try
	{
		reader4->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		exit(0);
	}
	
	//copy rgb image from reader to another image , needed to later get output
	mImgpRGBimage = mInputImageFileType::New();
	mImgpRGBimage->SetRegions( reader4->GetOutput()->GetLargestPossibleRegion() );
	mImgpRGBimage->CopyInformation( reader4->GetOutput() );
	mImgpRGBimage->Allocate();
	itk::ImageRegionConstIterator< mInputImageFileType > cit(reader4->GetOutput(),reader4->GetOutput()->GetLargestPossibleRegion());
	itk::ImageRegionIterator< mInputImageFileType > it(mImgpRGBimage,mImgpRGBimage->GetLargestPossibleRegion());;
	for(it.GoToBegin(),cit.GoToBegin();!it.IsAtEnd();++it,++cit)
	{
		it.Set(cit.Get());
	}

}

/**	@brief	clears pix vector
*/
void
filter::clear(void)
{
	for(unsigned int i=0;i<mvPix.size();++i )	
	{
		delete mvPix[i];
		mvPix[i] = NULL;
	}
}

/**	@brief	loads image specified by name
*	@param	sName name of file
*	@return	ITK image ptr to loaded image
*/
filter::mOutputImageType::Pointer
filter::load_image(std::string sName)
{
	typedef itk::ImageFileReader< mOutputImageType > ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	//load filtered image
	reader->SetFileName(sName.c_str());

	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		exit(0);
	}
	return reader->GetOutput();
}

/**	@brief load the green channel from a series of images and save the rgp into a temp image
*	@return pointer ITK image
*/
filter::mOutputImageType::Pointer
filter::load_channel(const std::string &crsName,const std::string &crsPath,const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth)
{

	//create reader
	mReaderType::Pointer reader = mReaderType::New();

	reader->SetFileNames( gen_names(crsPath.c_str(),
		crsName.c_str(),crsType.c_str(), crStart,crEnd,crWidth));

	//create adaptor
	mImageGreenAdaptorType::Pointer f_adaptor = mImageGreenAdaptorType::New();
	
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		exit(0);
	}

	mImgpRGBimage = reader->GetOutput();

	//mOutputImageType::Pointer greenChan = mOutputImageType::New();
	//mOutputImageType::Pointer blueChan = mOutputImageType::New();
	//mOutputImageType::Pointer redChan = mOutputImageType::New();
	//mOutputImageType::Pointer outputChan = mOutputImageType::New();

	//greenChan = load_green_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	//redChan = load_red_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	//blueChan = load_blue_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);

	//itk::ImageRegionConstIterator< mOutputImageType > cGreenIt(greenChan,greenChan->GetLargestPossibleRegion());

	//int greenSum = 0;
	//for(cGreenIt.GoToBegin();!cGreenIt.IsAtEnd();++cGreenIt)
	//{
	//	greenSum += cGreenIt.Get();
	//}

	//itk::ImageRegionConstIterator< mOutputImageType > cRedIt(redChan,redChan->GetLargestPossibleRegion());

	//int redSum = 0;
	//for(cRedIt.GoToBegin();!cRedIt.IsAtEnd();++cRedIt)
	//{
	//	redSum += cRedIt.Get();
	//}

	//itk::ImageRegionConstIterator< mOutputImageType > cBlueIt(blueChan,blueChan->GetLargestPossibleRegion());

	//int blueSum = 0;
	//for(cBlueIt.GoToBegin();!cBlueIt.IsAtEnd();++cBlueIt)
	//{
	//	blueSum += cBlueIt.Get();
	//}

	//if( redSum > greenSum && redSum > blueSum)
	//{
	//	std::cout << "Getting red channel" << std::endl;
	//	outputChan = redChan;
	//}
	//else if( blueSum > greenSum && blueSum > redSum)
	//{
	//	std::cout << "Getting blue channel" << std::endl;
	//	outputChan = blueChan;
	//}
	//else
	//{
	//	std::cout << "Getting green channel" << std::endl;
	//	outputChan = greenChan;
	//}
	mOutputImageType::Pointer outputChan = mOutputImageType::New();
	if( mChannel == 0 )
		outputChan = load_red_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	else if( mChannel == 1 )
		outputChan = load_green_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	else
		outputChan = load_blue_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	
	
	return outputChan;
}

filter::mOutputImageType::Pointer
filter::load_channel(const std::string & filename)
{

	//create reader
	mReaderFileType::Pointer reader = mReaderFileType::New();

	reader->SetFileName( filename );

	//create adaptor
	mImageGreenAdaptorType::Pointer f_adaptor = mImageGreenAdaptorType::New();
	
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in file reader: " << e << std::endl;
		exit(0);
	}

	mImgpRGBimage = reader->GetOutput();

	//mOutputImageType::Pointer greenChan = mOutputImageType::New();
	//mOutputImageType::Pointer blueChan = mOutputImageType::New();
	//mOutputImageType::Pointer redChan = mOutputImageType::New();
	//mOutputImageType::Pointer outputChan = mOutputImageType::New();

	//greenChan = load_green_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	//redChan = load_red_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);
	//blueChan = load_blue_channel(crsName,crsPath,crsType,crStart,crEnd,crWidth);

	//itk::ImageRegionConstIterator< mOutputImageType > cGreenIt(greenChan,greenChan->GetLargestPossibleRegion());

	//int greenSum = 0;
	//for(cGreenIt.GoToBegin();!cGreenIt.IsAtEnd();++cGreenIt)
	//{
	//	greenSum += cGreenIt.Get();
	//}

	//itk::ImageRegionConstIterator< mOutputImageType > cRedIt(redChan,redChan->GetLargestPossibleRegion());

	//int redSum = 0;
	//for(cRedIt.GoToBegin();!cRedIt.IsAtEnd();++cRedIt)
	//{
	//	redSum += cRedIt.Get();
	//}

	//itk::ImageRegionConstIterator< mOutputImageType > cBlueIt(blueChan,blueChan->GetLargestPossibleRegion());

	//int blueSum = 0;
	//for(cBlueIt.GoToBegin();!cBlueIt.IsAtEnd();++cBlueIt)
	//{
	//	blueSum += cBlueIt.Get();
	//}

	//if( redSum > greenSum && redSum > blueSum)
	//{
	//	std::cout << "Getting red channel" << std::endl;
	//	outputChan = redChan;
	//}
	//else if( blueSum > greenSum && blueSum > redSum)
	//{
	//	std::cout << "Getting blue channel" << std::endl;
	//	outputChan = blueChan;
	//}
	//else
	//{
	//	std::cout << "Getting green channel" << std::endl;
	//	outputChan = greenChan;
	//}
	mOutputImageType::Pointer outputChan = mOutputImageType::New();
	if( mChannel == 0 )
		outputChan = load_red_channel(filename);
	else if( mChannel == 1 )
		outputChan = load_green_channel(filename);
	else
		outputChan = load_blue_channel(filename);
	
	
	return outputChan;
}

/**	@brief load the green channel from a series of images and save the rgp into a temp image
*	@return pointer ITK image
*/
filter::mOutputImageType::Pointer
filter::load_blue_channel(const std::string &crsName,const std::string &crsPath,
						   const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth)
{

	//create reader
	mReaderType::Pointer reader = mReaderType::New();

	reader->SetFileNames( gen_names(crsPath.c_str(),
		crsName.c_str(),crsType.c_str(), crStart,crEnd,crWidth));

	//create adaptor
	mImageBlueAdaptorType::Pointer f_adaptor = mImageBlueAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerBlueAdaptorType::Pointer f_rescaler_adap = mRescalerBlueAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}
filter::mOutputImageType::Pointer
filter::load_blue_channel(const std::string &filename)
{

	//create reader
	mReaderFileType::Pointer reader = mReaderFileType::New();

	reader->SetFileName( filename );

	//create adaptor
	mImageBlueAdaptorType::Pointer f_adaptor = mImageBlueAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerBlueAdaptorType::Pointer f_rescaler_adap = mRescalerBlueAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}
/**	@brief load the green channel from a series of images and save the rgp into a temp image
*	@return pointer ITK image
*/
filter::mOutputImageType::Pointer
filter::load_red_channel(const std::string &crsName,const std::string &crsPath,
						   const std::string &crsType,const int &crStart,const int &crEnd,const int &crWidth)
{

	//create reader
	mReaderType::Pointer reader = mReaderType::New();

	reader->SetFileNames( gen_names(crsPath.c_str(),
		crsName.c_str(),crsType.c_str(), crStart,crEnd,crWidth));

	//create adaptor
	mImageRedAdaptorType::Pointer f_adaptor = mImageRedAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerRedAdaptorType::Pointer f_rescaler_adap = mRescalerRedAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}
filter::mOutputImageType::Pointer
filter::load_red_channel(const std::string &filename)
{

	//create reader
	mReaderFileType::Pointer reader = mReaderFileType::New();

	reader->SetFileName( filename);

	//create adaptor
	mImageRedAdaptorType::Pointer f_adaptor = mImageRedAdaptorType::New();
	
	f_adaptor->SetImage(reader->GetOutput());

	mRescalerRedAdaptorType::Pointer f_rescaler_adap = mRescalerRedAdaptorType::New();

	//set rescale values
	f_rescaler_adap->SetOutputMinimum(  0  );
	f_rescaler_adap->SetOutputMaximum( 255 );

	f_rescaler_adap->SetInput(f_adaptor);

	try
	{
		f_rescaler_adap->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in RescalerAdaptor: " << e << std::endl;
		exit(0);
	}
	return f_rescaler_adap->GetOutput();
}
/**	@brief project the segmentaion onto one plane
*	@param file file name to save image
*	@param cells vector of cell ptrs to project
*/
filter::mInputImageFileType::Pointer
filter::project_image(std::string &file,std::vector<cell*> cells)
{

	mInputGreenPixelFileType pt;
	pt.SetRed(0);
	pt.SetGreen(0);
	pt.SetBlue(0);

	mIteratorType cit(mImgpIntenisty,mImgpIntenisty->GetLargestPossibleRegion());
	mInputImageFileType::Pointer imProj = mInputImageFileType::New();

	imProj->SetRegions(mImgpIntenisty->GetLargestPossibleRegion());
	imProj->CopyInformation( mImgpIntenisty );
	imProj->Allocate();
	mOutputImageType::IndexType miIndex;
	mInputImageFileType::IndexType ind;


	for(int x=0; x<m_nImage_size[0]; x++ )
	{
		for(int y=0; y<m_nImage_size[1]; y++ )
		{
			int max = -1;
			for(int z=0; z<m_nImage_size[2]; z++ )
			{
				miIndex[0] = x;
				miIndex[1] = y;
				miIndex[2] = z;
				cit.SetIndex(miIndex);
				max = (cit.Get() > max) ? cit.Get() : max;
			}
			ind[0] = x;
			ind[1] = y;
			ind[2] = 0;
			if( mChannel == 0 )
				pt.SetRed(max);
			else if( mChannel == 1 )
				pt.SetGreen(max);
			else if( mChannel == 2 )
				pt.SetBlue(max);
			
			imProj->SetPixel(ind,pt);
		}
	}
	pt.SetRed(255);
	pt.SetGreen(0);
	pt.SetBlue(0);

	mInputGreenPixelFileType ptG;
	mInputGreenPixelFileType ptN;
	if( mChannel == 0 )
	{
		std::cout << "Labeling Red Channel" << std::endl;
		ptG.SetRed(0);
		ptG.SetGreen(255);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(0);
		ptN.SetBlue(255);
	}
	else if( mChannel == 1 )
	{
		std::cout << "Labeling Green Channel" << std::endl;
		ptG.SetRed(255);
		ptG.SetGreen(0);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(0);
		ptN.SetBlue(255);
	}
	else if( mChannel == 2 )
	{
		std::cout << "Labeling Blue Channel" << std::endl;
		ptG.SetRed(255);
		ptG.SetGreen(0);
		ptG.SetBlue(0);

		ptN.SetRed(0);
		ptN.SetGreen(255);
		ptN.SetBlue(0);
	}
	else
	{
		ptG.SetRed(255);
		ptG.SetGreen(255);
		ptG.SetBlue(255);

		ptN.SetRed(255);
		ptN.SetGreen(255);
		ptN.SetBlue(255);
	}
	//loop through all cells
	for(unsigned int i=0;i<cells.size();i++)
	{
		if(cells[i] == NULL)
			continue;
		if( cells[i]->mClass == 2 )//&& cells[i]->mIsTrain)
		{
			pt = ptG;
		}
		else if( cells[i]->mClass == 1 )// && cells[i]->mIsTrain )
		{
			pt = ptN;
		}
		else
		{
			continue;
		}

		for(unsigned int j=0;j<cells[i]->mvBoundPix2D.size();j++)
		{
			ind[0] = cells[i]->mvBoundPix2D[j]->x_;
			ind[1] = cells[i]->mvBoundPix2D[j]->y_;
			ind[2] = 0;
			imProj->SetPixel(ind,pt);
		}

	}
	typedef itk::ImageFileWriter< mInputImageFileType >  tWriterType;
	tWriterType::Pointer writer = tWriterType::New();

	writer->SetFileName( file.c_str() );

	writer->SetInput( imProj );
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & e)
	{
		std::cerr << "Exception in Writer: " << e << std::endl;
		exit(0);
	}
        return imProj;
}
