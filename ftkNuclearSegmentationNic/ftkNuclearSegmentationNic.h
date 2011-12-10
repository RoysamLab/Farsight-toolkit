// ############################################################################################################################################################################
#ifndef _ftknuclearsegmentationnic_h_
#define _ftknuclearsegmentationnic_h_
// ############################################################################################################################################################################



// ############################################################################################################################################################################

// Rules to edit this class
// comment your changes after the line you edit // 
// folow doxygen standard to comment the file
// Instantiate any new function to make sure it does not contains errors (templates do not check for errors unles instantiated)



// ############################################################################################################################################################################
#include <vector>
#include <iostream>
#include <map>


#include <ftkObject.h>
#include <ftkImage/ftkImage.h>


#include <itkImage.h>


namespace ftk{
	/**
	* namespace corresponding to the segmentation algorithm 
	* implemented by nicolas
	*/
	namespace ftkNucSecNic{
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
		*/
		template < typename inputPixelType >
		class ftkNuclearSegmentationNic
		{
		public:

			/**
			*	Image Input image type
			*/
			typedef itk::Image< inputPixelType, 2 > inputImageType_2;
			typedef itk::Image< inputPixelType, 3 > inputImageType_3;

			/** 
			*	Struct to store the parameters name and their values
			*/
			typedef struct { std::string name; int value; } parameter;

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

			/**
			*	Pass a pointer to the already loaded image	
			*/
			bool SetInput(ftk::Image::Pointer inImg, std::string fname, int chNumber = 0); 

			/**
			*	Set the parameters filename	
			*/
			void SetParameters(std::string paramFileName);

			/**
			*	Set the parameters according to the name
			*	@param value value of the given parameter (corresponding to name)
			*/
			void SetParameter(std::string name, int value);

			//void setBinarizeParameters();
			//void setSeedDetectionParameters();



			//	-----------------------------------------------------------------------------------------------------------------------
			// 0. PreProcess
			//	-----------------------------------------------------------------------------------------------------------------------

			// Anisotropic diffusion
			void setTime( unsigned int time = 0 );

			//	-----------------------------------------------------------------------------------------------------------------------
			// 1. Binarization
			//	-----------------------------------------------------------------------------------------------------------------------
			/**
			*	binarization ussing a mixture of poisson distributions.
			*	@param binarizeparam parameters for the binarization
			*	1. number of poisson distrubution (2, or 3) is none of this values the number will be dettermined automatically
			*/
			void ftkBinarizeMixPoisson(unsigned int numberBins = 128, bool getResultImg = false);//, const std::vector<double> binarizeparam);

			/**
			*	binarization ussing a mixture of gaussian distributions
			*	@param binarizeparam parameters for the binarization
			*	1. number of gaussian distrubution (2, or 3) is none of this values the number will be dettermined automatically
			*/
			//void ftkBinarizeMixGaussian(unsigned int numberBins = 128, bool getResultImg = false, const std::vector<double> binarizeparam){};

			/**
			*	binarization using otsu method
			*/
			//void ftkBinarizeOtsu(unsigned int numberBins = 128,bool getResultImg = false, const std::vector<double> binarizeparam){};


			//	-----------------------------------------------------------------------------------------------------------------------
			// 2. Seed Detection
			//	-----------------------------------------------------------------------------------------------------------------------

			//	-----------------------------------------------------------------------------------------------------------------------
			// 3. Clusstering
			//	-----------------------------------------------------------------------------------------------------------------------

			//	-----------------------------------------------------------------------------------------------------------------------
			// 4. Finalize !!!
			//	-----------------------------------------------------------------------------------------------------------------------


		protected:

			std::string					_paramFileName;					/*!< Name of the file containing the parameters */
			std::vector<std::string>	_paramNames;					/*!< Vector of parameters names */
			std::vector<parameter>		_vectParameters;/*_myParameters;*/					/*!< !!! Not Sure */

			std::string					_dataFileName;					/*!< Name of the file that the data came from */
			ftk::Image::Pointer			_dataImage;						/*!< The data image */
			int							_CH;							/*!< Use this channel from the dataImage for segmentation */

			std::string					_labelFilename;					/*!< Name of the file that is the label image !!! Not Sure */ 
			ftk::Image::Pointer			_labelImage;					/*!< Labeled image, result of segmentation */

			std::map< int, ftk::Object::Box >	_bBoxMap;				/*!< Geometric Info, Edit: Bounding boxes */
			std::map< int, ftk::Object::Point >	_centerMap;				/*!< Centroids */

			bool						_EditsNotSaved;					/*!< Flag, indicate if the edits where saved !!! */

			int							_currentTime;					/*!< Current Nucleus Editor time slider. */

			std::string					_errorMessage;					/*!< Sting that handles the error messages. */

			unsigned int				_T;								/*!< Time of the image been segmented. */

			unsigned int				_numConnComp;					/*!< Number of connected componnents. */

			/**
			*	{0,1,2,3,4} for the stages in a nuclear segmentation. This is a internal variable, that keeps track of the process, and stop the algorithm to perform a step, if the previos one was not completed succesfully.
			*/
			int _lastRunStep;




		private:
			//	-----------------------------------------------------------------------------------------------------------------------
			// 1. Binarization
			//	-----------------------------------------------------------------------------------------------------------------------
			






		};

	};

};




#include "ftkNuclearSegmentationNic.hxx"

#endif