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

#include <SegmentationCommon/ftkSegmentationResult.h>
#include <FTKImage/ftkImage.h>
#include <SegmentationCommon/ftkLabelImageToFeatures.h>
#include <yousef_core/yousef_seg.h>

namespace ftk
{ 
/** \class NuclearSegmentation
 *  \brief For storage of a complete nuclear segmentation 
 *   
 *  Handles the execution, result, and editing of a nuclear segmentation
 *  
 */
class NuclearSegmentation : public SegmentationResult
{
public:
	NuclearSegmentation(string projpath, string projname);	

	//********************************************************************************************
	//LEGACY FUNCTIONS TO INTERACT WITH MODULE WIDGET:
	string getPackageName(){ return string("nuclei"); };
	vector<string> getModuleNames() { return moduleNames; };
	ftk::Image* getDataImage(void){ return dataImage; };
	ftk::Image* getLabelImage(void){ return labelImage; };

	void setup(string imagefilename, string paramfilename);
	void executeModule(int);
	//*********************************************************************************************

	//For starting a new segmentation  -can set the parameters, then they will be in the xml for later 
	bool RunSegmentation(){return 0;};
	//Now for loading the segmentation results into Objects
	bool LoadFromResult(const char* dfile, const char* rfile);
	bool LoadFromMETA(std::string META_file, std::string header_file, std::string data_file, std::string label_file);
	void LoadClassInfoFromFile( std::string fName );
	bool LabelsToObjects(void);
	//bool CalculateFeatures(){return 0;};

	//Inherited, for loading a previously saved result:
	//bool RestoreFromXML();
	//Inherited, Writes segmentation info to XML
	//bool WriteToXML():
	//Inherited: Set/Get/Add for segmentation variables

	//Load up the data and result information into memory, Base class does not do this, it just gets the filename
	bool LoadData();				//Will load the data image into memory
	bool LoadLabel();				//Will load the label image into memory
	bool SaveLabel();				//Save the changes made to the label image
	bool SaveLabelByClass();		//Will save a different label image for each class
	bool SaveAll();					//Save all results to files (Label,XML,Seeds,META,etc);

	//Editing Functions  
	bool Split(int id){return 0;};
	int Merge(vector<int> ids);
	bool Delete(vector<int> ids);
	bool Add( Object::Point p ){return 0;};
	bool editsNotSaved;

private:

	bool FileExists(const char* fname);
	int GetObjectIndex(int objectID, string type);
	string TimeStamp();
	void ReassignLabels(vector<int> fromIds, int toId, ftk::Object::Box region);
	void ReassignLabel(int fromId, int toId);
	Object GetNewObject(int id, IntrinsicFeatures *features );
	void LoadAssociationsFromFile(std::string fName);
	ftk::Object::Box ExtremaBox(vector<int> ids);

	//Load up the data and result information into memory, Base class does not do this, it just gets the filename
	bool UnloadData(){return 0;};	//Release the data image from memory
	bool UnloadLabel(){return 0;}; //Release the label image from memory
//********************************************************************************************
//********************************************************************************************
	//LEGACY FUNCTIONS TO INTERACT WITH MODULE WIDGET:
	void initConstants(void);
	void createFTKLabelImg(int* data,int numColumns, int numRows, int numStacks);

	//LEGACY VARIABLES:
	vector<string> moduleNames;
	ftk::Image *dataImage;
	ftk::Image *labelImage;
	yousef_nucleus_seg *NucleusSeg;
//********************************************************************************************
//********************************************************************************************

}; // end NuclearSegmentation

}  // end namespace ftk

#endif	// end __ftkNuclearSegmentation_h

