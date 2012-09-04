// ############################################################################################################################################################################
#ifndef _MaxClustering_h_
#define _MaxClustering_h_
// ############################################################################################################################################################################

#include <iostream>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <math.h>

// ITK INCLUDES
#include <itkImage.h>
#include <itkImageFileWriter.h>
// #include <itkSignedMaurerDistanceMapImageFilter.h>
// #include <itkAbsImageFilter.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
// #include <itkImageDuplicator.h>
// 
// #include "itklaplacianrecursivegaussianimagefilternew.h"


// OMP
#include "omp.h"


//#include "../ftkNuclearSegmentationNic.h"


namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
 	namespace nucSecNic{
		template < typename inputPixelType, typename binaryPixelType = unsigned short, typename seedDetecPixelType = unsigned short, typename loGResponsePixelType = double >
		class MaxClustering
		{
		public:
			
			typedef itk::Image< inputPixelType, 3 > inputImageType; // not sure about puting typename
			typedef itk::Image< binaryPixelType, 3 > binaryImageType;
			typedef itk::Image< seedDetecPixelType, 3 > seedDetectImageType;  	// For seed image
// 			
// 			typedef double loGResponsePixelType;
// 			typedef double distMapPixelType;
// 			
// // 			typedef itk::Image< float, 3 > distMapImageType;
// 			typedef itk::Image< distMapPixelType, 3 > distMapImageType;
			typedef itk::Image< loGResponsePixelType, 3 > loGResponseImageType;
// 			
			typedef double inputLoGIPixelType;
			typedef itk::Image< inputLoGIPixelType, 3 > inputLoGImageType;
			
// 			
			typedef itk::LaplacianRecursiveGaussianImageFilter< binaryImageType, loGResponseImageType >  LoGFilterType;
// 			typedef itk::LaplacianRecursiveGaussianImageFilterNew< inputLoGImageType, loGResponseImageType >  LoGFilterType;
// 			
			
			// FOR NOW I HARD CODE HERE THE MAXCLUSTERING RESULT, TO THE SEED DETECTING PIXEL TYPE
// 			typedef seedDetecPixelType maxClustPixelType;
			typedef itk::Image< seedDetecPixelType, 3 > maxClustImageType;
			
			
			MaxClustering(){
				std::cout << "Created Max Clustering";
			};
			~MaxClustering(){};
			void setParameters( bool getResultImg_MaxClustering = false );
			void setInput( const ftk::Image::Info* info, typename inputImageType::Pointer inputImage, typename binaryImageType::Pointer binaryImage, typename seedDetectImageType::Pointer seedDetecImage, typename loGResponseImageType::Pointer loGResponseImage );
			void runMaxClustering();
			typename itk::Image< seedDetecPixelType, 3 >::Pointer getMaxClustImage();

		private:
// 			void distMap_2D();
// 			void loGResponse( long long actualScale, typename loGResponseImageType::Pointer &loGResponse );
// 			
// // 			void runMinErrorThresholding();
// 			
			const Image::Info				*_info;						/*!< Information of the input image */
			long long					_numRows;					/*!< Number of rows in the input image */
			long long					_numColumns;					/*!< Number of colums in the input image */
			long long					_numStacks;					/*!< Number of stacks in the input image */
			long long					_totNumPixels;					/*!< Total number of pixels */
			double						_spacing_XY;
			double						_spacing_Z;
			unsigned int					_maxNumberIterations;
			inputPixelType					_maxValueInputPixelType;			/*!< Maximum value of the given input pixel type */
			binaryPixelType					_maxValueBinaryPixelType;
			seedDetecPixelType				_maxValueSeedDetectPixelType;
// 			
// 			long long 					_seedDetectMinScale;
// 			long long 					_seedDetectMaxScale;
			bool 						_getResultImg_MaxClustering;
// 			
// // 			unsigned int 					_numPoissonDist;				/*!< Number of poissson distribution */
// // 			double						_alpha_B; // alpha_1 ???
// // 			double						_alpha_1;
// // 			double						_alpha_F; // alpha_2 ??
// // 			double						_alpha_2; 
// // 			double						_P_I;
// 			
// // 			double						_alpha_C; // alpha_1 ???
// // 			double						_alpha_3; 
// // 			double						_P_I2;
// 			
			
			typename inputImageType::Pointer		_inputImage;
			typename inputImageType::PixelType 		*_inputImageArray; 
			typename binaryImageType::Pointer		_binaryImage;
			typename binaryImageType::PixelType		*_binaryImageArray;
			typename seedDetectImageType::Pointer		_seedDetectImage;
			typename seedDetectImageType::PixelType		*_seedDetectImageArray;
// 			
// 			typename distMapImageType::Pointer		_distMapImage;
// 			typename distMapImageType::PixelType		*_distMapImageArray;
// 			
			typename loGResponseImageType::Pointer		_loGResponseImage;
			typename loGResponseImageType::PixelType	*_loGResponseImageArray;
// 			
// 			typename loGResponseImageType::Pointer		_loGResponseImageMax;
// 			typename loGResponseImageType::PixelType	*_loGResponseImageMaxArray;
// 			
// 			typename inputLoGImageType::Pointer		_inputLoGImage;
// 			typename inputLoGImageType::PixelType		*_inputLoGImageArray;
// 			
// 			typename LoGFilterType::Pointer 		_laplacian;
// 			
// 			distMapPixelType			_maxMapDistance;
			
			typename maxClustImageType::Pointer		_maxClustImage;
			typename maxClustImageType::PixelType		*_maxClustImageArray;
			
			
		};
	};
};

#include "MaxClustering.hxx"

#endif



