/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkNuclearSegmentation_h
#define __ftkNuclearSegmentation_h

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <ftkImage/ftkImage.h>
#include <SegmentationCommon/ftkLabelImageToFeatures.h>
#include <yousef_core/yousef_seg.h>
#include <tinyxml/tinyxml.h>
#include <ftkObject.h>
#include <map>
#include <set>

namespace ftk
{ 
/** \class NuclearSegmentation
 *  \brief For storage of a complete nuclear segmentation 
 *   
 *  Handles the execution, result, and editing of a nuclear segmentation
 *  
 */

class ftkPoint
{
public:
	ftkPoint(int x_in, int y_in, int z_in): x(x_in), y(y_in), z(z_in) {};
private:
	int x;
	int y;
	int z;
};
class NuclearSegmentation
{
public:
	NuclearSegmentation();

	//This is for beginning a completely new segmenation:
	bool SetInputs(std::string datafile, std::string paramfile);		
	bool LoadData();											//Will load the data image into memory
	bool Binarize();											//Will binarize the data image
	bool DetectSeeds();											//If binarization has been done it will detect seeds
	bool RunClustering();										//Will use binary image and seeds to do initial clustering
	bool Finalize();											//Will finilize the output using alpha expansion
	bool GetResultImage();										//Gets the result of last module and puts it in labelImage
	bool SaveOutput();											//Save the output of the last step executed (image format)
	//Segmentation is basically done at this point (hopefully), now move on to calculating the features and classification:
	bool LabelsToObjects(void);									//Will compute Intrinsic Features and create objects from the data and results images
	bool LoadAssociationsFromFile(std::string fName);			//Add the Associative Features to the objects
	bool LoadClassInfoFromFile( std::string fName );			//Add the Class info to the objects
	//Save in various formats:
	bool WriteToXML(std::string filename);
	bool WriteToMETA(std::string filename);
	bool WriteToLibSVM(std::string filename);
	bool SaveLabel();														//Save the changes made to the label image
	bool SaveLabelByClass();												//Will save a different label image for each class
	

	//We may also want to restore from previously found results:
	bool RestoreFromXML(std::string filename);					//Complete Restore from XML file
	bool LoadLabel();											//Load just the label image if the filename is already known
	bool LoadFromImages(std::string dfile, std::string rfile);	//Load from images -> then convert to objects
	void LoadFromDAT(std::string dfile, std::string rfile);		//Load from .dat -> then convert to objects
	bool LoadFromMETA(std::string META_file, std::string header_file, std::string data_file, std::string label_file);

	//Editing Functions  
	std::vector< int > Split(ftk::Object::Point P1, ftk::Object::Point P2);
	int Merge(vector<int> ids);
	bool Delete(vector<int> ids);
	bool Add( Object::Point p ){return 0;};
	bool editsNotSaved;				//Will be true if edits have been made and not saved to file.

	std::string GetErrorMessage() { return errorMessage; };
	std::string GetDataFilename() { return dataFilename; };
	std::string GetParamFilename() { return paramFilename; };

	std::vector<Object>* GetObjectsPtr(){ return &myObjects; };	//Use pointer
	Object* GetObjectPtr(int id);
	std::vector<std::string> GetFeatureNames(){ return featureNames; };

	ftk::Image::Pointer getDataImage(void){ return dataImage; };
	ftk::Image::Pointer getLabelImage(void){ return labelImage; };
	//*********************************************************************************************

	//Added by Yousef on 04-08-2009
	//This function will run graph coloring and will assign different colors for touching objects
	//For now, it will just write the list of labels into a text file, but this should be relaxed later
	bool RunGraphColoring(std::string labelname, std::string filename);	//Run Graph coloring on label image and save adjacency file as filename

private:
	std::string dataFilename;	//the filename of the data image		(full path)
	std::string labelFilename;	//the filename of the label image		(full path)
	std::string paramFilename;	//the filename of the parameter file	(full path)

	std::string errorMessage;

	ftk::Image::Pointer dataImage;
	ftk::Image::Pointer labelImage;
	yousef_nucleus_seg *NucleusSeg;
	int lastRunStep;

	typedef struct { string name; int value; } Parameter;
	std::vector<Parameter> myParameters;
	std::vector<Object> myObjects;
	std::vector<std::string> featureNames;
	std::vector<int> classes;

	//To help keep track of objects and their IDs:
	int GetObjectIndex(int objectID, std::string type);
	std::map<int,int> IdToIndexMap;
	int maxID;

	//Utilities for parsing XML:
	Object parseObject(TiXmlElement *object);
	Object::Point parseCenter(TiXmlElement *centerElement);
	Object::Box parseBound(TiXmlElement *boundElement);
	std::vector<float> parseFeatures(TiXmlElement *featureElement);

	//Utilities for writing XML:
	TiXmlElement *GetObjectElement(Object object);
	std::string NumToString(int i);
	std::string NumToString(double d);
	std::string NumToString(double d, int p);

	//General Utilites:
	bool FileExists(std::string filename);
	
	std::string TimeStamp();
	void ReassignLabels(std::vector<int> fromIds, int toId, ftk::Object::Box region);
	void ReassignLabel(int fromId, int toId);
	Object GetNewObject(int id, IntrinsicFeatures *features );
	ftk::Object::Box ExtremaBox(std::vector<int> ids);

//********************************************************************************************
//********************************************************************************************

}; // end NuclearSegmentation

}  // end namespace ftk

#endif	// end __ftkNuclearSegmentation_h

