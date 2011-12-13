// ############################################################################################################################################################################
#ifndef _BinarizeMixPoisson_h_
#define _BinarizeMixPoisson_h_
// ############################################################################################################################################################################

#include <iostream>
#include <fstream>
#include <limits>
#include <stdio.h>
#include <math.h>

// ITK INCLUDES
#include <itkImage.h>
#include <itkImageFileWriter.h>


// GRAPH CUTS
//#include "../GraphCutsBoykov/maxflow-v3.01/graph.h"
#include "graph.h"

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
			double computePoissonProb( int intensity, double alpha);
			void graphCuts_2D();
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
			inputPixelType					_maxValueInputPixelType;			/*!< Maximum value of the given input pixel type */
			binaryPixelType					_maxValueBinaryPixelType;
			
			unsigned int 					_numberBins_mixPoisson;
			bool 						_getResultImg_mixPoisson;
			double						_sigmaNeighCost;				/*!< Neighborhood parameter cost. */
			
			unsigned int 					_numPoissonDist;				/*!< Number of poissson distribution */
			double						_alpha_B; // alpha_1 ???
			double						_alpha_1;
			double						_alpha_F; // alpha_2 ??
			double						_alpha_2; 
			double						_P_I;
			
			double						_alpha_C; // alpha_1 ???
			double						_alpha_3; 
			double						_P_I2;
			
			
			typename inputImageType::Pointer		_inputImage;
			typename inputImageType::PixelType 		*_inputImageArray; 
			typename binaryImageType::Pointer		_binaryImage;
			typename binaryImageType::PixelType		*_binaryImageArray;
			
// 			itk::image<bool>::pointer binary;

			
		};
	};
};

#include "BinarizeMixPoisson.hxx"

#endif



