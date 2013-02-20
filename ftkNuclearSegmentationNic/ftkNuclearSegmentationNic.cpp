#include "ftkNuclearSegmentationNic.h"

// ############################################################################################################################################################################
ftk::ftkNucSecNic::ftkNuclearSegmentationNic< typename inputPixelType >::ftkNuclearSegmentationNic()
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
}

// ############################################################################################################################################################################
ftk::ftkNucSecNic::ftkNuclearSegmentationNic::~ftkNuclearSegmentationNic()
{

	ReleaseMemory();
}

// ############################################################################################################################################################################
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::ResetRealeaseAll()
{
	_dataFileName.clear();
	_dataImage = NULL;
	_channelNumber = 0;
	_labelImage = NULL;
	_bBoxMap.clear();
	_centerMap.clear();
	_paramFileName.clear();
	_vectParameters.clear();
	_EditsNotSaved = false;
	_lastRunStep = 0;

	ReleaseMemory();
}

// ############################################################################################################################################################################
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::ReleaseMemory()
{
	//itk::ImageIOBase::
	// !!! need to work on this

	/*if(NucleusSeg)
	{
		delete NucleusSeg;
		
		NucleusSeg = NULL;
	}*/
}

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
bool ftk::ftkNucSecNic::ftkNuclearSegmentationNic::SetInput(ftk::Image::Pointer inImg, std::string fname, int chNumber)
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
	_channelNumber = chNumber;

	typedef 

	return true;
}

// ############################################################################################################################################################################
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::SetParameters(std::string paramFileName)
{
	_paramFileName = paramFileName;
}

// ############################################################################################################################################################################
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::SetParameter(std::string name, int value)
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
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::setTime( unsigned int time )
{
	_time = time;
}


// ############################################################################################################################################################################
void ftk::ftkNucSecNic::ftkNuclearSegmentationNic::ftkBinarizeMixPoisson(unsigned int numberBins, bool getResultImg)//, const std::vector<double> binarizeparam)
{

	//, bool getResultImg = false,const std::vector<double> binarizeparam)
//{
	//ftk::Image::DataType dataType = dataImage->GetImageInfo()->dataType;
	int t = 0;


	_dataImage

	_dataImage->GetImageInfo()->dataType nic;
	ftk::Image::itkPixelType
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