// ############################################################################################################################################################################
#ifndef _ftknuclearsegmentationnic_h_
#define _ftknuclearsegmentationnic_h_
// ############################################################################################################################################################################



// ############################################################################################################################################################################

// Rules to edit this class
// comment your changes after the line you edit // 
// folow doxygen standard to comment the file
// Instantiate any new function to make sure it does not contains errors (templates do not check for errors unles instantiated)
// All the new parameters should have a default value
// The order of the functions definition and declaration should be keept
// Also, the order of the parameters definitions should be keept
// ...
// ...

// ############################################################################################################################################################################


// STRUCT USED BY ALL THE FILES
#include "ftkGlobalIncludes.h"
#include "frkGlobalStructs.h"



// NUCLEAR SEGMENTATION NIC INCLUDES
#include "BinarizeMixPoisson/BinarizeMixPoisson.h"
//#include "BinarizeMixGaussian/BinarizeMixGaussian.h"
//#include "BinarizeOtsu/BinarizeOtsu.h"
//#include "OtherMethod/OtherMethod.h"
#include "SeedDetection/SeedDetectionLoG.h"
#include "MaxClustering/MaxClustering.h"
#include "AlphaExpansion/AlphaExpansion.h"
//#include "RadVoting/RadVoting.h"
//#include "CellShapeIni/CellShapeIni.h"
//#include "GVFforce/GVFforce.h"
//#include "AnotherForce/AnotherForce.h"
//#include "LevelSet/LevelSet.h"

namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
	namespace nucSecNic{
		//!  nuclear segmentation class
		/*!
		* this class implements the nuclear segmentation in basically 4 steps. binarization, seed detection, clustering, refine using alpha expansion.
		* details can be found in "improved automatic detection and segmentation of cell nuclei in histopathology images"
		* @see ftkbinarizemixpoisson 
		* @see ftkbinarizemixgaussian 
		* @see ftkbinarizeottsu 
		* @see ftkseeddetection 
		* @see ftkclustering
		* @see ftkfinalalphaexpansion
		* The number of been normally has to be a power of 2, otherwise "noise" can affect the performance
		*/
		template < typename inputPixelType, typename binaryPixelType = unsigned short, typename seedDetectPixelType = unsigned short, typename labelPixelType = unsigned short > // Label pixel type is especially designed in case images with more cells
		class ftkNuclearSegmentationNic
		{
		public:
			
			typedef double loGResponsePixelType;

			/**
			*	Image Input image type
			*/
			typedef itk::Image< inputPixelType, 2 > inputImageType_2;
			typedef itk::Image< inputPixelType, 3 > inputImageType_3;

			/**
			*	Binary image type
			*/
			typedef itk::Image< binaryPixelType, 3 > binaryImageType;
			
			/**
			 *	Seed Detection
			 */
			typedef itk::Image< seedDetectPixelType, 3 > seedDetectImageType;
			
			/**
			 *	LoG Response
			 */
			typedef itk::Image< loGResponsePixelType, 3 > loGResponseImageType;
			
			/**
			 * 	Label Image Type
			 */
			typedef itk::Image< seedDetectPixelType, 2 > labelImageType_2;	// THE TYPENAME LABEL PIXEL, WILL BE ELIMINATED, ONLY THE SEEDPIZEL TYPE WILL BE USED
			typedef itk::Image< seedDetectPixelType, 3 > labelImageType_3;
			
			/** 
			*	Struct to store the parameters name and their values
			*/
			typedef struct { std::string paramName; int value; } parameterInt;
			typedef struct { std::string paramName; double value; } parameterDoubl;

			/**
			*	constructor.
			*	details.
			*/
			ftkNuclearSegmentationNic();
			/**
			* a destructor.
			* a more elaborate description of the destructor.
			*/
			~ftkNuclearSegmentationNic();

			// ------------------------------------------------------------------ ERRASE MEMORY ------------------------------------------------------------------
			/**
			*	Reset all the pointers and deletes all the data
			*/
			void ResetRealeaseAll();

			/**
			*	This function errases the pointers to all the internal memory of the class. The class creates and destroys the pointers to the binary image, the seed Image, etc... when the user ask for a copy, the itkDuplicator function is used, to create a copy
			*	@see (all the pointers that allocate the image)
			*/
			void ReleaseMemory();

			/**
			*	Load the input image from a file
			*/
			//bool LoadInput(std::string fname, int chNumber = 0);							

			// ------------------------------------------------------------------ SET GLOBAL PARAMETERS ------------------------------------------------------------------
			/**
			*	Pass a pointer to the already loaded image	
			*	\param CH channel to be segmented
			* 	\param T time point to be segmented
			*/
			bool setInput(ftk::Image::Pointer dataImage, std::string dataFileName, unsigned int CH = 0, unsigned int T = 0); 

			/**
			 * 	Set Pipeline Number
			 */
			void setPipelineNumber( unsigned int pipelineNumber );
			
			/**
			*	Set the parameters filename	
			*/
			void setParameters(std::string paramFileName);

			/**
			*	Set the parameters according to the name
			*	@param value value of the given parameter (corresponding to name)
			*/
			void setParameter(std::string paramName, double value);
			
			/**
			 * 	Set up all the parameters to a default value
			 */
			void setParametersDefault();
			
			// ################################################################## SET PARTICULAR PARAMETERS ACCORDING TO THE MODULES OF THE PIPELINES ##################################################################
			// ------------------------------------------------------------------ PIPELINE 1 ------------------------------------------------------------------
			/**
			 *	1.1 Set BinarizePoisson Parameters
			 */
			void setBinarizeMixPoissonParameters_1(unsigned int numberBins_mixPoisson_1, bool getResultImg_mixPoisson_1 = false);
			void runBinarizeMixPoisson_1();					/*!< Run BinarizePoisson. */
			typename binaryImageType::Pointer getImageBinarizeMixPoisson_1(){};				/*!< Get . */
			
			/**
			 * 	1.2 Set BinarizeGaussian Parameters
			 */
			void setBinarizeMixGaussianParameters_1(unsigned int numberBins_mixGaussian_1, bool getResultImg_mixGaussian_1 = false);
			void runBinarizeMixGaussian_1(){};				/*!< Run BinarizeGaussian. */
			void getImageBinarizeMixGaussian_1(){};				/*!< Get BinarizeGaussian. */
			
			/**
			 * 	1.3 Set BinarizeOtsu Parameters
			 */
			void setBinarizeOtsuParameters_1(unsigned int numberBins_otsuBina_1, bool getResultImg_otsuBina_1 = false);
			void runBinarizeOtsu_1(){}; 					/*!< Run BinarizeOtsu. */
			void getImageBinarizeOtsu_1(){};				/*!< Get BinarizeOtsu. */
			
			/**
			 * 	1.4 Set Parameters of other binarization method (Adaptive, local, hard threshdold)
			 */
			void setOtherMethodParameters_1(unsigned int numberBins_otherMethod_1, bool getResultImg_otherMethod_1 = false);
			void runOtherMethod_1(){}; 					/*!< Run other binarization method (Adaptive, local, hard threshdold). */
			void getImageOtherMethod_1(){};					/*!< Get other binarization method (Adaptive, local, hard threshdold). */
					
			/**
			 * 	2.1 Set Parameters for Seed Detection
			 */
			void setSeedDetectionParameters_1(long long seedDetectMinScale, long long seedDetectMaxScale, bool getResultImg_seedDetect_1 = false);
			void runSeedDetection_1(); 					/*!< Run Seed Detection. */
			void getImageSeedDetection_1(){};				/*!< Get Seed Detection. */
			
			/**
			 * 	3.1 Set Max Clustering parameters
			 */
			void setMaxClusteringParameters_1(bool getResultImg_maxClust_1 = false);
			void runMaxClustering_1(); 					/*!< Run Max Clustering. */
			void getImageMaxClustering_1(){};				/*!< Get Max Clustering. */
			
			/**
			 * 	4.1 Run alpha expansion
			 */
			void setAlphaExpansionParameters_1(bool getResultImg_alphaExp_1 = false);
			void runAlphaExpansion_1(); 					/*!< Run alpha expansion. */
			void getImageAlphaExpansion_1(){};				/*!< Get alpha expansion. */
			
			// ------------------------------------------------------------------ PIPELINE 2 ------------------------------------------------------------------
			/**
			 * 	1.1 Set the parameter for seed detection by voting
			 */
			void setRadVotingParameters_2(unsigned int voteMinScale_2, unsigned int voteMaxScale_2, bool getResultImg_radVoting_2 = false);
			void runRadVoting_2(){}; 					/*!< Run seed detection by voting. */
			void getImageRadVoting_2(){};					/*!< Get seed detection by voting. */
			
			/**
			 * 	2.1 Set Initialization of the first cell shape Parameters
			 */
			void setCellShapeIniParameters_2( unsigned int numOutliers_cellShapeIni_2, bool getResultImg_cellShapeIni_2 = false);
			void runCellShapeIni_2(){}; 					/*!< Run Initialization of the first cell shape. */
			void getImageCellShapeIni_2(){};				/*!< Get Initialization of the first cell shape. */
			
			/**
			 * 	3.1 GVF External force !!! not sure is making this independent
			 */
			void setGVFforceParameters_2( double mu_GVFforce_2, double iter_GVFforce_2, bool getResultImg_GVFforce_2 = false);
			void runGVFforce_2(){}; 					/*!< Run GVF External force. */
			void getImageGVFforce_2(){};					/*!< Get GVF External force. */
			
			/**
			 * 	3.2 Another force !!! not sure is making this independent
			 */
			void setAnotherForceParameters_2( double alpha_1_AnotherForce_2, double alpha_2_AnotherForce_2, bool getResultImg_AnotherForce_2 = false);
			void runAnotherForce_2(){}; 					/*!< Run Another force. */
			void getImageAnotherForce_2(){};				/*!< Get Another force. */
			
			/**
			 * 	4.1 Set up the parameters of the Level Set equation
			 */
			void setLevelSetParameters_2( double betha_1_LevelSet_2, double betha_2_LevelSet_2);
			void runLevelSet_2(){}; 					/*!< Run Level Set equation. */
			void getImageLevelSet_2(){};					/*!< Get Level Set equation. */
			
			// ------------------------------------------------------------------ PIPELINE 3 ------------------------------------------------------------------
			/**
			 * 	1.1 Level Sets
			 */
			void setLevlSetParameters_3( double betha_1_LevlSet_3, double betha_2_LevlSet_3 );
			void runLevlSet_3(){};						/*!< Run Level Sets. */
			void getImageLevlSet_3(){};					/*!< Get Level Sets. */
			// ################################################################## DONE PIPELINES PARAMETERS ##################################################################
			// ###############################################################################################################################################################
			
			// ################################################################## MAIN FUNCTIONS ##################################################################
			/**
			 * 	Thiis is the main function, ALWAYS run this function inside a try and catch. This function actually Segments the Image 
			 *	@return This function does not 
			 */
			void runSegmentImage();
			void getSegmentImage(){};					/*!< Get Segmentation result. */
			
			
		protected:

			// ------------------------------------------------------------------ GLOBAL PARAMETERS ------------------------------------------------------------------
			std::string					_paramFileName;					/*!< Name of the file containing the parameters */
			std::vector<std::string>			_paramNames;					/*!< Vector of parameters names */
			std::vector<parameterDoubl>			_vectParameters;/*_myParameters;*/		/*!< !!! Not Sure */

			std::string					_dataFileName;					/*!< Name of the file that the data came from */
			ftk::Image::Pointer				_dataImage;					/*!< The data image */
			int						_CH;						/*!< Use this channel from the dataImage for segmentation */
			unsigned int					_T;						/*!< Time of the image been segmented. */
			const ftk::Image::Info 				*_info;						/*!< Information of the input image */
			unsigned int					_numRows;					/*!< Number of rows in the input image */
			unsigned int					_numColumns;					/*!< Number of colums in the input image */
			unsigned int					_numStacks;					/*!< Number of stacks in the input image */
			long long					_totNumPixels;					/*!< Total number of pixels */
			long long					_maxValueInputPixelType;			/*!< Maximum value of the given input pixel type */

			std::string					_labelFilename;					/*!< Name of the file that is the label image !!! Not Sure */ 
			ftk::Image::Pointer				_labelImage;					/*!< Labeled image, result of segmentation */
			typename itk::Image< inputPixelType, 3 >::Pointer	_itkPointerToInputImage_3;			/*!< Itk Pointer to the image to be segmented */

			std::map< int, ftk::Object::Box >		_bBoxMap;					/*!< Geometric Info, Edit: Bounding boxes */
			std::map< int, ftk::Object::Point >		_centerMap;					/*!< Centroids */

			bool						_EditsNotSaved;					/*!< Flag, indicate if the edits where saved !!! */

			//int						_currentTime;					/*!< Current Nucleus Editor time slider. */

			std::string					_errorMessage;					/*!< Sting that handles the error messages. */
  
			unsigned int					_numConnComp;					/*!< Number of connected componnents. */

			/**
			*	{0,1,2,3,4} for the stages in a nuclear segmentation. This is a internal variable, that keeps track of the process, and stop the algorithm to perform a step, if the previos one was not completed succesfully.
			*/
			int 						_lastRunStep;
			
			unsigned int 					_pipelineNumber;				/*!< Pipeline number of segmentation. */
			
			ConnComp					*_myConnComp;
			typename binaryImageType::Pointer		_binaryImage;
			typename seedDetectImageType::Pointer		_seedDetectImage;
			typename loGResponseImageType::Pointer		_loGResponseImage;
			typename labelImageType_3::Pointer		_maxClustImage;
			
// 			loGResponse
			
		private:
			
			// ################################################################## PARTICULAR PARAMETERS ACCORDING TO THE MODULES OF THE PIPELINES ##################################################################
			// ------------------------------------------------------------------ PIPELINE 1 ------------------------------------------------------------------
			BinarizeMixPoisson< inputPixelType, binaryPixelType >	*_objBinarizeMixPoisson_1;			/*!< 1.1 Object BinarizePoisson. */
			unsigned int 					_numberBins_mixPoisson_1;			/*!< 1.1 BinarizePoisson Parameters. */
			bool 						_getResultImg_mixPoisson_1;			/*!< 1.1 BinarizePoisson Parameters. */
			bool						_use_mixPoisson_1;				/*!< 1.1 Use up BinarizePoisson Parameters. See Note 1 */
			bool						_run_mixPoisson_1;				/*!< 1.1 Run BinarizePoisson Parameters. See Note 1 */
			
			//BinarizeMixPoisson_1				*_objBinarizeMixPoisson_1;			
			unsigned int 					_numberBins_mixGaussian_1;			/*!< 1.2 BinarizeGaussian Parameters. */
			bool 						_getResultImg_mixGaussian_1;			/*!< 1.2 BinarizeGaussian Parameters. */
			bool						_use_mixGaussian_1;				/*!< 1.2 Use up BinarizeGaussian Parameters. See Note 1 */
			bool						_run_mixGaussian_1;				/*!< 1.2 Run BinarizeGaussian Parameters. See Note 1 */
			
			//BinarizeOtsu_1
			unsigned int 					_numberBins_otsuBina_1;				/*!< 1.3 BinarizeOtsu Parameters */
			bool 						_getResultImg_otsuBina_1;			/*!< 1.3 BinarizeOtsu Parameters */
			bool						_use_otsuBina_1;				/*!< 1.3 Use up BinarizeOtsu Parameters. See Note 1 */
			bool						_run_otsuBina_1;				/*!< 1.3 Run BinarizeOtsu Parameters. See Note 1 */
			
			//OtherMethod_1
			unsigned int 					_numberBins_otherMethod_1;			/*!< 1.4 Parameters of other binarization method (Adaptive, local, hard threshdold) */
			bool 						_getResultImg_otherMethod_1;			/*!< 1.4 Parameters of other binarization method (Adaptive, local, hard threshdold) */
			bool						_use_otherMethod_1;				/*!< 1.4 Use up Binarize Other method Parameters. See Note 1 */
			bool						_run_otherMethod_1;				/*!< 1.4 Run Binarize Other method Parameters. See Note 1 */
			
			SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > *_objSeedDetectionLoG_1;	/*!< Seed Detection Object */
			long long					_seedDetectMinScale;				/*!< 2.1 Seed Detection Min Scale */
			long long					_seedDetectMaxScale;				/*!< 2.1 Seed Detection Max Scale */
			bool 						_getResultImg_seedDetect_1;			/*!< 2.1 Seed Detection */
			
			MaxClustering< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > *_objMaxClustering_1;	/*!< 3.1 Max Clustering Object */
			bool 						_getResultImg_maxClust_1;			/*!< 3.1 Max Clustering parameters */
			
			AlphaExpansion< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > *_objAlphaExpansion_1;
			bool 						_getResultImg_alphaExp_1;			/*!< 4.1 Run alpha expansion */
			
			// ------------------------------------------------------------------ PIPELINE 2 ------------------------------------------------------------------
			unsigned int 					_voteMinScale_2;				/*!< 1. parameter for seed detection by voting */
			unsigned int 					_voteMaxScale_2;				/*!< 1. parameter for seed detection by voting */
			bool 						_getResultImg_radVoting_2;			/*!< 1. parameter for seed detection by voting */
			
			unsigned int 					_numOutliers_cellShapeIni_2;			/*!< 2. Initialization of the first cell shape */
			bool 						_getResultImg_cellShapeIni_2;			/*!< 2. Initialization of the first cell shape */
			
			double						_mu_GVFforce_2;					/*!< 3.1 GVF External force !!! not sure is making this independent */
			double						_iter_GVFforce_2;				/*!< 3.1 GVF External force !!! not sure is making this independent */
			bool						_getResultImg_GVFforce_2;			/*!< 3.1 GVF External force !!! not sure is making this independent */
			
			double						_alpha_1_AnotherForce_2;			/*!< 3.2 Another force !!! not sure is making this independent */
			double						_alpha_2_AnotherForce_2;			/*!< 3.2 Another force !!! not sure is making this independent */
			bool 						_getResultImg_AnotherForce_2;			/*!< 3.2 Another force !!! not sure is making this independent */
			
			double						_betha_1_LevelSet_2;				/*!< 4. The parameters of the Level Set equation */
			double						_betha_2_LevelSet_2;				/*!< 4. The parameters of the Level Set equation */
			
			// ------------------------------------------------------------------ PIPELINE 3 ------------------------------------------------------------------
			double						_betha_1_LevlSet_3;				/*!< 1. Level Sets */
			double						_betha_2_LevlSet_3;				/*!< 1. Level Sets */
			
		};
	};	// End of nucSecNic namespace
};		// End of ftk namespace




#include "ftkNuclearSegmentationNic.hxx"

#endif
