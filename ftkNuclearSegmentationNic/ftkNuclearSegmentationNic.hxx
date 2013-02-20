//#include "ftkNuclearSegmentationNic.h"

#include <itkImage.h>

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::ftkNuclearSegmentationNic()
{
	//NucleusSeg = NULL;
	this->ResetRealeaseAll();
	
	// Default parameters
	this->setParametersDefault();
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::~ftkNuclearSegmentationNic()
{

	ReleaseMemory();
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::ResetRealeaseAll()
{
	_paramFileName.clear();
	_vectParameters.clear();
  
	_dataFileName.clear();
	_dataImage = NULL;
	_CH = 0;
  
	_labelFilename.clear();
	_labelImage = NULL;
  
	_bBoxMap.clear();
	_centerMap.clear();
  
	_EditsNotSaved = false;
  
	_T = 0;
  
	_numConnComp = 0;
  
	_pipelineNumber = 0;
  
	_lastRunStep = 0;
  
	ReleaseMemory();
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::ReleaseMemory()
{
	// This should realease all the memory
	//itk::ImageIOBase::
	// !!! need to work on this

	/*if(NucleusSeg)
	{
		delete NucleusSeg;
		
		NucleusSeg = NULL;
	};*/
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//bool ftk::nucSecNic::ftkNuclearSegmentationNic::LoadInput( std::string fname, int chNumber )
//{
//	ftk::Image::Pointer tmpImg = ftk::Image::New();
//	if(!tmpImg->LoadFile(fname))	//Load for display
//	{
//		errorMessage = "Data Image failed to load";
//		return false;
//	};
//	return this->SetInput(tmpImg,fname,chNumber);
//};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
bool ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setInput(ftk::Image::Pointer dataImage, std::string dataFileName, unsigned int CH, unsigned int T)
{
	if( CH > dataImage->GetImageInfo()->numChannels)
	{
		_errorMessage = "Channel does not exist";
		return false;
	};
	//if(inImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
	//{
	//	errorMessage = "module only works for 8-bit unsigned char input data";
	//	return false;
	//};
	_dataFileName = dataFileName;
	_dataImage = dataImage;
	_CH = CH;
	_T = T;
	
	// There is a problem in the ftk::Image, the time should be a unsigned int (no negative time) so for that reason I create this temporal varialbe T_temp, to call the function, otherwise the function will be called using _T
	int T_temp = _T;
	_itkPointerToInputImage_3 = _dataImage->GetItkPtr< inputPixelType >( T_temp, _CH);
	
	_info = _dataImage->GetImageInfo();
	
	_numRows = _info->numRows;
	_numColumns = _info->numColumns;
	_numStacks = _info->numZSlices;
	
	_totNumPixels = (long long)_numRows*(long long)_numColumns*(long long)_numStacks;
	
	std::cout << std::endl << "Set up and input image of size, Rows: " << _numRows << ", Col: " << _numColumns << ", Slices: " << _numStacks << ", Voxels: " << _totNumPixels;
	
	_maxValueInputPixelType = std::numeric_limits<inputPixelType>::max();
	
	return true;
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setParameters(std::string paramFileName)
{
	_paramFileName = paramFileName;
	
	// This class should read the file, and then put those values in the respective function parameters

// 	char achBuffer[1024];
// 	char achBuffer2[1024];
// 	
// 	int iCounter = 0;
// 	
// 	std::ifstream inFile(_paramFileName.c_str());
// 
// 	if (! inFile)
// 	{
// 		cout << "Fatal Error: Could not open parameters file in the second time" << pFname
// 			<< " .... Terminating Program." << endl;
// 	}
// 	while (inFile)
// 	{
// 		inFile.getline(achBuffer, 1024, '\n');
// 		if (achBuffer[0] != '\0' && achBuffer[0] != '!')
// 		{
// 			strcpy(achBuffer2, achBuffer);
// 			pchStr = strtok(achBuffer, "\t ");
// 
// 			if (pchStr)
// 			{
// 				strcpy(m_pData[m_iNumOfElements].m_pName, pchStr);
// 			}
// 
// 			pchStr = strtok(NULL, "\t ");
// 
// 			if (pchStr)
// 			{
// 				if (*pchStr == ':')
// 					pchStr = strtok(NULL, "\t " );
// 
// 				m_pData[m_iNumOfElements].m_pValue= atoi(pchStr);
// 			}
// 			m_iNumOfElements++;
// 		}
// 	}
// 	inFile.close();	
// 
	// Reorder the parameters
// 	for(int i=0; i<=iCounter; i++)
// 	{
// 		if(!strcmp(m_pData[i].m_pName,"high_sensitivity"))
// 			params[0] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"LoG_size"))
// 			params[1] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"min_scale"))
// 			params[2] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"max_scale"))
// 			params[3] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"xy_clustering_res"))
// 			params[4] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"z_clustering_res"))
// 			params[5] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"finalize_segmentation"))
// 			params[6] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"sampling_ratio_XY_to_Z"))
// 			params[7] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"Use_Distance_Map"))
// 			params[8] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"refinement_range"))
// 			params[9] = m_pData[i].m_pValue;
// 		else if(!strcmp(m_pData[i].m_pName,"min_object_size"))
// 			params[10] = m_pData[i].m_pValue;
// 		else
// 			continue;
// 	}
// 
// 	setParams(params);
	
	
	
};

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setParameter(std::string paramName, double value)
{
	bool found = false;
	for(int i=0; i<(int)_vectParameters.size(); ++i)
	{
		if(_vectParameters.at(i).paramName == paramName)
		{
			_vectParameters.at(i).value = value;
			found = true;
			break;
		};
	};

	if(!found)		//Add the parameter
	{
		parameterDoubl newP;
		newP.paramName = paramName;
		newP.value = value;
		_vectParameters.push_back(newP);
	};
};

// // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
// void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setTime( unsigned int T )
// {
// 	_T = T;
// 	time = 1;
// };

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setParametersDefault( )
{
	_paramNames.push_back("high_sensitivity");		//0
	_paramNames.push_back("LoG_size");			//1
	_paramNames.push_back("min_scale");			//2
	_paramNames.push_back("max_scale");			//3
	_paramNames.push_back("xy_clustering_res");		//4
	_paramNames.push_back("z_clustering_res");		//5
	_paramNames.push_back("finalize_segmentation");		//6
	_paramNames.push_back("Sampling_ratio_XY_to_Z");	//7
	_paramNames.push_back("Use_Distance_Map");		//8
	_paramNames.push_back("refinement_range");		//9
	_paramNames.push_back("min_object_size");		//10
	
// ------------------------------------------------------------------ GLOBAL PARAMETER ------------------------------------------------------------------
	
	parameterDoubl paramDoubObj;
	
	paramDoubObj.paramName = _paramNames.at(0);
	paramDoubObj.value = 0;
	_vectParameters.push_back(paramDoubObj);
	std::cout << std::endl << _paramNames.at(0);
	paramDoubObj.paramName = _paramNames.at(1);
	paramDoubObj.value = 30;
	_vectParameters.push_back(paramDoubObj);
	std::cout << std::endl << _paramNames.at(1);
	paramDoubObj.paramName = _paramNames.at(2);
	paramDoubObj.value = 5;
	_vectParameters.push_back(paramDoubObj);
	std::cout << std::endl << _paramNames.at(2);
	paramDoubObj.paramName = _paramNames.at(3);
	paramDoubObj.value = 8;
	_vectParameters.push_back(paramDoubObj);
	std::cout << std::endl << _paramNames.at(3);
	paramDoubObj.paramName = _paramNames.at(4);
	paramDoubObj.value = 5;
	_vectParameters.push_back(paramDoubObj);
	std::cout << std::endl << _paramNames.at(4);
	paramDoubObj.paramName = _paramNames.at(5);
	paramDoubObj.value = 2;
	_vectParameters.push_back(paramDoubObj);
	
	paramDoubObj.paramName = _paramNames.at(6);
	paramDoubObj.value = 0;
	_vectParameters.push_back(paramDoubObj);
	
	paramDoubObj.paramName = _paramNames.at(7);
	paramDoubObj.value = 2;
	_vectParameters.push_back(paramDoubObj);
	
	paramDoubObj.paramName = _paramNames.at(8);
	paramDoubObj.value = 1;
	_vectParameters.push_back(paramDoubObj);
	
	paramDoubObj.paramName = _paramNames.at(9);
	paramDoubObj.value = 6;
	_vectParameters.push_back(paramDoubObj);
	
	paramDoubObj.paramName = _paramNames.at(10);
	paramDoubObj.value = 100;
	_vectParameters.push_back(paramDoubObj);
	
	
// ------------------------------------------------------------------ PIPELINE 1 ------------------------------------------------------------------
	_numberBins_mixPoisson_1 = 128;
	_getResultImg_mixPoisson_1 = false;
	_use_mixPoisson_1 = true;
	_run_mixPoisson_1 = false;


	_numberBins_mixGaussian_1 = 128;
	_getResultImg_mixGaussian_1 = false;
	_use_mixGaussian_1 = false;
	_run_mixGaussian_1 = false;

	_numberBins_otsuBina_1 = 128;
	_getResultImg_otsuBina_1 = false;
	_use_otsuBina_1 = false;
	_run_otsuBina_1 = false;

	_numberBins_otherMethod_1 = 128;
	_getResultImg_otherMethod_1 = false;
	_use_otherMethod_1 = false;
	_run_otherMethod_1 = false;

	_seedDetectMinScale = 4;
	_seedDetectMaxScale = 12;
	_getResultImg_seedDetect_1 = false;

	_getResultImg_maxClust_1 = false;

	_getResultImg_alphaExp_1 = false;

// ------------------------------------------------------------------ PIPELINE 2 ------------------------------------------------------------------
	_voteMinScale_2 = 1;
	_voteMaxScale_2 = 5;
	_getResultImg_radVoting_2 = false;

	_numOutliers_cellShapeIni_2 = 2;
	_getResultImg_cellShapeIni_2 = false;

	_mu_GVFforce_2 = 0.5;
	_iter_GVFforce_2 = 0.5;
	_getResultImg_GVFforce_2 = 0.7;

	_alpha_1_AnotherForce_2 = 0.8;
	_alpha_2_AnotherForce_2 = 0.9;
	_getResultImg_AnotherForce_2 = 1.0;

	_betha_1_LevelSet_2 = 1.1;
	_betha_2_LevelSet_2 = 1.2;

// ------------------------------------------------------------------ PIPELINE 3 ------------------------------------------------------------------
	_betha_1_LevlSet_3 = 1.3;
	_betha_2_LevlSet_3 = 1.4;

};

// ############################################################################################################################################################################
// ------------------------------------------------------------------ PIPELINE 1 ------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setBinarizeMixPoissonParameters_1(unsigned int numberBins_mixPoisson_1, bool getResultImg_mixPoisson_1)
{
	_numberBins_mixPoisson_1 = numberBins_mixPoisson_1;
	_getResultImg_mixPoisson_1 = getResultImg_mixPoisson_1;
	
	//_getResultImg_mixPoisson_1
	//_objBinarizeMixPoisson_1->setParameters( _numberBins_mixPoisson_1,_getResultImg_mixPoisson_1 )
	_use_mixPoisson_1 = true;
	_use_mixGaussian_1 = false;
	_use_otsuBina_1 = false;
	_use_otherMethod_1 = false;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::runBinarizeMixPoisson_1()
{
 	_objBinarizeMixPoisson_1 = new BinarizeMixPoisson< inputPixelType, binaryPixelType > ();
// 	BinarizeMixPoisson< inputPixelType, binaryPixelType >* objBinarizeMixPoisson_1 = new BinarizeMixPoisson< inputPixelType, binaryPixelType > ();
  	_objBinarizeMixPoisson_1->setParameters( _numberBins_mixPoisson_1, _getResultImg_mixPoisson_1 );
 	_objBinarizeMixPoisson_1->setInput( _info, _itkPointerToInputImage_3 );
 	_objBinarizeMixPoisson_1->runBinarization();
	_binaryImage = _objBinarizeMixPoisson_1->getBinarizedImage();
	_myConnComp = _objBinarizeMixPoisson_1->getConnComponents();
	
// 	int yy;
// 	cin >> yy;
	// delte the binarization object ??? not sure about this part
	delete _objBinarizeMixPoisson_1;
	
// 	cin>>yy;
	
};

// template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
// typename binaryImageType::Pointer ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::getImageBinarizeMixPoisson_1()
// {
// 	return _binaryImage;
// }

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setBinarizeMixGaussianParameters_1(unsigned int numberBins_mixGaussian_1, bool getResultImg_mixGaussian_1)
{
	_numberBins_mixGaussian_1 = numberBins_mixGaussian_1;
	_getResultImg_mixGaussian_1 = getResultImg_mixGaussian_1;
	
	_use_mixPoisson_1 = false;
	_use_mixGaussian_1 = true;
	_use_otsuBina_1 = false;
	_use_otherMethod_1 = false;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setBinarizeOtsuParameters_1(unsigned int numberBins_otsuBina_1, bool getResultImg_otsuBina_1){
	_numberBins_otsuBina_1 = numberBins_otsuBina_1;
	_getResultImg_otsuBina_1 = getResultImg_otsuBina_1;
	
	_use_mixPoisson_1 = false;
	_use_mixGaussian_1 = false;
	_use_otsuBina_1 = true;
	_use_otherMethod_1 = false;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setOtherMethodParameters_1(unsigned int numberBins_otherMethod_1, bool getResultImg_otherMethod_1){
	_numberBins_otherMethod_1 = numberBins_otherMethod_1;
	_getResultImg_otherMethod_1 = getResultImg_otherMethod_1;
	
	_use_mixPoisson_1 = false;
	_use_mixGaussian_1 = false;
	_use_otsuBina_1 = false;
	_use_otherMethod_1 = true;
};




template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setSeedDetectionParameters_1(long long seedDetectMinScale, long long seedDetectMaxScale, bool getResultImg_seedDetect_1){
	_seedDetectMinScale = seedDetectMinScale;
	_seedDetectMaxScale = seedDetectMaxScale;
	getResultImg_seedDetect_1 = getResultImg_seedDetect_1; 
	
};



template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::runSeedDetection_1(){
	_objSeedDetectionLoG_1 = new SeedDetectionLoG< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > ();
	_objSeedDetectionLoG_1->setParameters( _seedDetectMinScale, _seedDetectMaxScale, _getResultImg_seedDetect_1 );
	_objSeedDetectionLoG_1->setInput( _info, _itkPointerToInputImage_3, _binaryImage );
	
	_objSeedDetectionLoG_1->runSeedDetection();
	
	_seedDetectImage = _objSeedDetectionLoG_1->getSeedDetectImage();
	_loGResponseImage = _objSeedDetectionLoG_1->getLoGResponseImage();
	
	delete _objSeedDetectionLoG_1;
	
};




template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setMaxClusteringParameters_1(bool getResultImg_maxClust_1){
	_getResultImg_maxClust_1 = getResultImg_maxClust_1;
};


template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::runMaxClustering_1(  ){
	_objMaxClustering_1 = new MaxClustering< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > ();
	_objMaxClustering_1->setParameters( _getResultImg_seedDetect_1 );
	_objMaxClustering_1->setInput( _info, _itkPointerToInputImage_3, _binaryImage, _seedDetectImage, _loGResponseImage );
// 	
	_objMaxClustering_1->runMaxClustering();
	
	_maxClustImage = _objMaxClustering_1->getMaxClustImage();
	
	delete _objMaxClustering_1;
	
};




template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setAlphaExpansionParameters_1(bool getResultImg_alphaExp_1){
	_getResultImg_alphaExp_1 = getResultImg_alphaExp_1;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::runAlphaExpansion_1(){
	
	_objAlphaExpansion_1 = new AlphaExpansion< inputPixelType, binaryPixelType, seedDetectPixelType, loGResponsePixelType > ();
// 	_objMaxClustering_1->setParameters( _getResultImg_seedDetect_1 );
// 	_objMaxClustering_1->setInput( _info, _itkPointerToInputImage_3, _binaryImage, _seedDetectImage, _loGResponseImage );
// // 	
// 	_objMaxClustering_1->runMaxClustering();
// 	
// 	_maxClustImage = _objMaxClustering_1->getMaxClustImage();
// 	
// 	delete _objMaxClustering_1;
	
};

// ------------------------------------------------------------------ PIPELINE 2 ------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setRadVotingParameters_2(unsigned int voteMinScale_2, unsigned int voteMaxScale_2, bool getResultImg_radVoting_2){
	_voteMinScale_2 = voteMinScale_2;
	_voteMaxScale_2 = voteMaxScale_2;
	_getResultImg_radVoting_2 = getResultImg_radVoting_2;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setCellShapeIniParameters_2( unsigned int numOutliers_cellShapeIni_2, bool getResultImg_cellShapeIni_2){
	_numOutliers_cellShapeIni_2 = numOutliers_cellShapeIni_2;
	_getResultImg_cellShapeIni_2 = getResultImg_cellShapeIni_2;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setGVFforceParameters_2( double mu_GVFforce_2, double iter_GVFforce_2, bool getResultImg_GVFforce_2 ){
	_mu_GVFforce_2 = mu_GVFforce_2;
	_iter_GVFforce_2 = iter_GVFforce_2;
	_getResultImg_GVFforce_2 = getResultImg_GVFforce_2;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setAnotherForceParameters_2( double alpha_1_AnotherForce_2, double alpha_2_AnotherForce_2, bool getResultImg_AnotherForce_2 ){
	_alpha_1_AnotherForce_2 = alpha_1_AnotherForce_2;
	_alpha_2_AnotherForce_2 = alpha_2_AnotherForce_2;
	_getResultImg_AnotherForce_2 = getResultImg_AnotherForce_2;
};

template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setLevelSetParameters_2( double betha_1_LevelSet_2, double betha_2_LevelSet_2 ){
	_betha_1_LevelSet_2 = betha_1_LevelSet_2;
	_betha_2_LevelSet_2 = betha_2_LevelSet_2;
};

// ------------------------------------------------------------------ PIPELINE 3 ------------------------------------------------------------------
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::setLevlSetParameters_3( double betha_1_LevlSet_3, double betha_2_LevlSet_3 ){
	_betha_1_LevlSet_3 = betha_1_LevlSet_3;
	_betha_2_LevlSet_3 = betha_2_LevlSet_3;
};

// ###############################################################################################################################################################
template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::runSegmentImage(){
	
};
















// // ############################################################################################################################################################################
// template < typename inputPixelType, typename binaryPixelType, typename seedDetectPixelType, typename labelPixelType >
// void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::ftkBinarizeMixPoisson(unsigned int numberBins, bool getResultImg)
// //void ftk::nucSecNic::ftkNuclearSegmentationNic< inputPixelType, binaryPixelType, seedDetectPixelType, labelPixelType >::ftkBinarizeMixPoisson(unsigned int numberBins, bool getResultImg)//, const std::vector<double> 
// {
// 	// !!! Missing test for the correct data 
// 	typename inputImageType_3::Pointer dptr = _dataImage->GetItkPtr<inputPixelType>(_T,_CH,0);
// 	
// 	//itk::Image
// 	//itk::I
// 	//inputPixelType
// 
// 	//, bool getResultImg = false,const std::vector<double> binarizeparam)
// //{
// 	//ftk::Image::DataType dataType = dataImage->GetImageInfo()->dataType;
// 	//int t = 0;
// 
// 
// 	//_dataImage
// 
// 	//_dataImage->GetImageInfo()->dataType nic;
// 	//ftk::Image::itkPixelType
// 	//_dataImage->GetPixelT
// 	//	_dataImage->GetDataTypeITK
// 	//itk::Image<	_dataImage->GetImageInfo()->dataType,3>::Pointer dptr = _dataImage->GetItkPtr<_dataImage->GetImageInfo()->dataType>(_time,_channelNumber,0);	//Expects grayscale image 
// 	//if( 
// 	//	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(t,channelNumber,0);	//Expects grayscale image
// 
// };

// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################
// ############################################################################################################################################################################