// ############################################################################################################################################################################
#ifndef _SeedDetectionLoG_h_
#define _SeedDetectionLoG_h_
// ############################################################################################################################################################################

#include <iostream>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <math.h>

// ITK INCLUDES
#include <itkImage.h>
#include <itkImageFileWriter.h>


//#include "../ftkNuclearSegmentationNic.h"


namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
 	namespace nucSecNic{
		template < typename inputPixelType, typename binaryPixelType = unsigned short, typename seedDetecPixelType = unsigned short >
		class SeedDetectionLoG
		{
		public:
			
			typedef itk::Image< inputPixelType, 3 > inputImageType; // not sure about puting typename
			typedef itk::Image< binaryPixelType, 3 > binaryImageType;
			typedef itk::Image< seedDetecPixelType, 3 > seedDetectImageType;
			
			SeedDetectionLoG(){
				std::cout << "Created LoG";
			};
			~SeedDetectionLoG(){};
			void setParameters( long long seedDetectMinScale, long long seedDetectMaxScale, bool getResultImg_seedDetect );
			void setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage );
			void runSeedDetection();
// 			typename itk::Image< seedDetecPixelType, 3 >::Pointer getSeedDetecImage();

		private:
			
// 			void runMinErrorThresholding();
			
			const Image::Info				*_info;						/*!< Information of the input image */
			long long					_numRows;					/*!< Number of rows in the input image */
			long long					_numColumns;					/*!< Number of colums in the input image */
			long long					_numStacks;					/*!< Number of stacks in the input image */
			long long					_totNumPixels;					/*!< Total number of pixels */
			inputPixelType					_maxValueInputPixelType;			/*!< Maximum value of the given input pixel type */
			binaryPixelType					_maxValueBinaryPixelType;
			seedDetecPixelType				_maxValueSeedDetectPixelType;
			
			long long 					_seedDetectMinScale;
			long long 					_seedDetectMaxScale;
			bool 						_getResultImg_seedDetect;
			
// 			unsigned int 					_numPoissonDist;				/*!< Number of poissson distribution */
// 			double						_alpha_B; // alpha_1 ???
// 			double						_alpha_1;
// 			double						_alpha_F; // alpha_2 ??
// 			double						_alpha_2; 
// 			double						_P_I;
			
// 			double						_alpha_C; // alpha_1 ???
// 			double						_alpha_3; 
// 			double						_P_I2;
			
			
			typename inputImageType::Pointer		_inputImage;
			typename inputImageType::PixelType 		*_inputImageArray; 
			typename binaryImageType::Pointer		_binaryImage;
			typename binaryImageType::PixelType		*_binaryImageArray;
			typename seedDetectImageType::Pointer		_seedDetectImage;
			typename seedDetectImageType::PixelType		*_seedDetectImageArray;
			
			
		};
	};
};

#include "SeedDetectionLoG.hxx"

#endif



