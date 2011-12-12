// ############################################################################################################################################################################
#ifndef _BinarizeMixPoisson_h_
#define _BinarizeMixPoisson_h_
// ############################################################################################################################################################################

#include <iostream>

// ITK INCLUDES
#include <itkImage.h>

//#include "../ftkNuclearSegmentationNic.h"


namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
 	namespace nucSecNic{
		template < typename inputPixelType, typename binaryPixelType = unsigned short >
		class BinarizeMixPoisson
		{
		public:
			
			typedef itk::Image< inputPixelType, 3 > inputImageType; // not sure about puting typename
			typedef itk::Image< binaryPixelType, 3 > binaryImageType;
			
			BinarizeMixPoisson(){
				std::cout << "Created";
			};
			~BinarizeMixPoisson(){};
			void setParameters( unsigned int numberBins_mixPoisson, bool getResultImg_mixPoisson );
			void setInput( const ftk::Image::Info* info, typename itk::Image< inputPixelType, 3 >::Pointer inputImage );
			void runBinarization();
// 			{
// 				
// 				new binary;
// 			}
// 			itk::image<bool> getOutpu;
		private:
			
			void runMinErrorThresholding();
			
			
			
			const Image::Info				*_info;						/*!< Information of the input image */
			unsigned int					_numRows;					/*!< Number of rows in the input image */
			unsigned int					_numColumns;					/*!< Number of colums in the input image */
			unsigned int					_numStacks;					/*!< Number of stacks in the input image */
			long long					_totNumPixels;					/*!< Total number of pixels */
			
			unsigned int 					_numberBins_mixPoisson;
			bool 						_getResultImg_mixPoisson;
			
			double						_alpha_B;
			double						_alpha_F;
			double						_P_I;
			
			double						_alpha_C;
			double						_P_I2;
			
			
			typename inputImageType::Pointer		_inputImage;
			typename binaryImageType::Pointer		_binaryImage;
			
// 			itk::image<bool>::pointer binary;

			
		};
	};
};

#include "BinarizeMixPoisson.hxx"

#endif



