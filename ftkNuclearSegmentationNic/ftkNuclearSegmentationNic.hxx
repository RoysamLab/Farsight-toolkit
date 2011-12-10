//#include "ftkNuclearSegmentationNic.h"

#include <itkImage.h>

// ############################################################################################################################################################################
template< typename inputPixelType >
ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::ftkNuclearSegmentationNic()
{
	//NucleusSeg = NULL;
	this->ResetRealeaseAll();

	_paramNames.push_back("High_sensitivity");
	_paramNames.push_back("LoG_size");
	_paramNames.push_back("min_scale");
	_paramNames.push_back("max_scale");
	_paramNames.push_back("xy_clustering_res");
	_paramNames.push_back("z_clustering_res");
	_paramNames.push_back("Finalize_segmentation");
	_paramNames.push_back("Sampling_ratio_XY_to_Z");
	_paramNames.push_back("Use_Distance_Map");
	_paramNames.push_back("Refinement_range");
	_paramNames.push_back("min_object_size");
	_currentTime = 0;
};

// ############################################################################################################################################################################
template< typename inputPixelType >
ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::~ftkNuclearSegmentationNic()
{

	ReleaseMemory();
};

// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::ResetRealeaseAll()
{
	_dataFileName.clear();
	_dataImage = NULL;
	_CH = 0;
	_labelImage = NULL;
	_bBoxMap.clear();
	_centerMap.clear();
	_paramFileName.clear();
	_vectParameters.clear();
	_EditsNotSaved = false;
	_lastRunStep = 0;
	_T = 0;

	ReleaseMemory();
};

// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::ReleaseMemory()
{
	//itk::ImageIOBase::
	// !!! need to work on this

	/*if(NucleusSeg)
	{
		delete NucleusSeg;
		
		NucleusSeg = NULL;
	}*/
};

//// ############################################################################################################################################################################
//bool ftk::ftkNucSecNic::ftkNuclearSegmentationNic::LoadInput( std::string fname, int chNumber )
//{
//	ftk::Image::Pointer tmpImg = ftk::Image::New();
//	if(!tmpImg->LoadFile(fname))	//Load for display
//	{
//		errorMessage = "Data Image failed to load";
//		return false;
//	}
//	return this->SetInput(tmpImg,fname,chNumber);
//}

// ############################################################################################################################################################################
template< typename inputPixelType >
bool ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::SetInput(ftk::Image::Pointer inImg, std::string fname, int chNumber)
{
	if(chNumber > inImg->GetImageInfo()->numChannels)
	{
		_errorMessage = "Channel does not exist";
		return false;
	}
	//if(inImg->GetImageInfo()->dataType != itk::ImageIOBase::UCHAR)
	//{
	//	errorMessage = "module only works for 8-bit unsigned char input data";
	//	return false;
	//}
	_dataFileName = fname;
	_dataImage = inImg;
	_CH = chNumber;


	return true;
};

// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::SetParameters(std::string paramFileName)
{
	_paramFileName = paramFileName;
}

// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::SetParameter(std::string name, int value)
{
	bool found = false;
	for(int i=0; i<(int)_vectParameters.size(); ++i)
	{
		if(_vectParameters.at(i).name == name)
		{
			_vectParameters.at(i).value = value;
			found = true;
			break;
		}
	}

	if(!found)		//Add the parameter
	{
		parameter newP;
		newP.name = name;
		newP.value = value;
		_vectParameters.push_back(newP);
	}
}

// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::setTime( unsigned int time )
{
	_T = time;
}


// ############################################################################################################################################################################
template< typename inputPixelType >
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::ftkBinarizeMixPoisson(unsigned int numberBins, bool getResultImg)
//void ftk::ftkNucSecNic::ftkNuclearSegmentationNic< inputPixelType >::ftkBinarizeMixPoisson(unsigned int numberBins, bool getResultImg)//, const std::vector<double> 
{
	// !!! Missing test for the correct data
	//inputImageType_3::Pointer dptr = _dataImage->GetItkPtr<inputPixelType>(_T,_CN,0);
	
	//itk::Image
	//itk::I
	//inputPixelType

	//, bool getResultImg = false,const std::vector<double> binarizeparam)
//{
	//ftk::Image::DataType dataType = dataImage->GetImageInfo()->dataType;
	//int t = 0;


	//_dataImage

	//_dataImage->GetImageInfo()->dataType nic;
	//ftk::Image::itkPixelType
	//_dataImage->GetPixelT
	//	_dataImage->GetDataTypeITK
	//itk::Image<	_dataImage->GetImageInfo()->dataType,3>::Pointer dptr = _dataImage->GetItkPtr<_dataImage->GetImageInfo()->dataType>(_time,_channelNumber,0);	//Expects grayscale image 
	//if( 
	//	unsigned char *dptr = dataImage->GetSlicePtr<unsigned char>(t,channelNumber,0);	//Expects grayscale image

}

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