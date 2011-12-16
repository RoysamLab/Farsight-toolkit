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
			void graphCuts_3D();
			typename itk::Image< binaryPixelType, 3 >::Pointer getBinarizedImage();
			ConnComp* getConnComponents();

		private:
			
			void runMinErrorThresholding();
			
			const Image::Info				*_info;						/*!< Information of the input image */
			long long					_numRows;					/*!< Number of rows in the input image */
			long long					_numColumns;					/*!< Number of colums in the input image */
			long long					_numStacks;					/*!< Number of stacks in the input image */
			long long					_totNumPixels;					/*!< Total number of pixels */
			double						_spacing_XY;
			double						_spacing_Z;
			unsigned int					_connComponentsConnectivity;			/*!< Componets connectivity */
			unsigned int 					_minObjectSize;					/*!< Minimum object size */
			unsigned int 					_numConnectedComponents;
			ConnComp					*_myConnComp;
				
			inputPixelType					_maxValueInputPixelType;			/*!< Maximum value of the given input pixel type */
			binaryPixelType					_maxValueBinaryPixelType;
			
			unsigned int 					_numberBins_mixPoisson;
			bool 						_getResultImg_mixPoisson;
			double						_sigmaNeighCost;				/*!< Neighborhood parameter cost. */
			double						_wNeigh;					/*!< W weight for neighborhood. */
			
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
			
			
		};
	};
};

#include "BinarizeMixPoisson.hxx"

#endif



