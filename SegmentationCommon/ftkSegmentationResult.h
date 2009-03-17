/*=========================================================================

  Program:   Farsight Biological Image Segmentation and Visualization Toolkit
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: 0.00 $

=========================================================================*/
#ifndef __ftkSegmentationResult_h
#define __ftkSegmentationResult_h

#include <tinyxml/tinyxml.h>
#include <ftkObject.h>
#include <map>

namespace ftk
{ 
/** \class SegmentationResult
 *  \brief For storage of a complete segmentation result, XML, etc.
 *  
 *  Represents data stored in an XML FILE
 *  
 *  Inherit from this class to extend the functionality for a specific Segmentation type 
 * (to add edit functions, etc)
 */
class SegmentationResult
{
public:
	typedef struct { string name; int value; } Parameter;

	SegmentationResult(string projpath, string projname);

	bool RestoreFromXML();
	bool WriteToXML();
	bool WriteToMETA();

	string GetErrorMessage() { return errorMessage; };

	//Get and Set Functions : don't need to use this in inherited function
	string GetProgramName(){ return programName; };
	string GetProjectPath(){ return projectPath; };
	string GetProjectName(){ return projectName; };
	string PrependProjectPath(string file);
	vector<string> GetDataFilenames(){ return dataFilenames; };
	vector<string> GetResultFilenames(){ return resultFilenames; };
	vector<Parameter> GetParameters(){ return myParameters; };
	vector<Object>* GetObjectsPtr(){ return &myObjects; };	//Use pointer
	Object* GetObjectPtr(int id);
	vector<string> GetFeatureNames(){ return featureNames; };

	void SetProgramName(string name){ programName = name; };
	void SetProjectPath(string path){ projectPath = path; };
	void SetProjectName(string name){ projectName = name; };
	void SetAssociationFile(string f){ associationFile = f; };
	void SetClassFile(string f){ classFile = f; };
	void AddDataFile(string fname){ dataFilenames.push_back(fname); };
	void ClearDataFilenames(){ dataFilenames.clear(); };
	void AddResultFile(string fname){ resultFilenames.push_back(fname); };
	void ClearResultFilenames(){ resultFilenames.clear(); };
	void AddParameter(Parameter p){ myParameters.push_back(p); };
	void ClearParameters(){ myParameters.clear(); };
	//OBJECTS USE PTR
	void AddFeatureName(string name){ featureNames.push_back(name); };
	void ClearFeatureNames(){ featureNames.clear(); };

protected:
	string programName;  //The segmentation algorithm/program that was run
	string projectPath;  //Location of all files for this project
	string projectName;  //A string that is common to all filenames for this project (xml file is projectName.xml)
	vector<string> dataFilenames;	//the filenames of the data images  (no path, includes extension)
	vector<string> resultFilenames;	//the filenames of the result (no path, includes extension) -label image for nuclei, xml tree for ...
	string associationFile;
	string classFile;

	vector<Parameter> myParameters;
	vector<Object> myObjects;
	vector<string> featureNames;

	std::map<int,int> IdToIndexMap;
	int maxID;
	string errorMessage;

	string NumToString(int i);
	string NumToString(double d);
	string NumToString(double d, int p);

private:
	string XmlFilename(bool);

	Object parseObject(TiXmlElement *object);
	Object::Point parseCenter(TiXmlElement *centerElement);
	Object::Box parseBound(TiXmlElement *boundElement);
	vector<float> parseFeatures(TiXmlElement *featureElement);
	TiXmlElement *GetObjectElement(Object object);

}; // end SegmentationResult

}  // end namespace ftk

#endif	// end __ftkSegmentationResult_h
