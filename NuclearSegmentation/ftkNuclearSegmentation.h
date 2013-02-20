/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

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

#ifdef USE_TRACKING
#include <CellTrackerLib/ftkTrackFeatures.h>
#endif

#include <ftkImage/ftkImage.h>
#include <ftkFeatures/ftkLabelImageToFeatures.h>
#include <ftkCommon/ftkUtils.h>
#include <ftkFeatures/ftkObject.h>
#include <yousef_core/yousef_seg.h>

#include <vtkSmartPointer.h>
#include <vtkTable.h>

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
class NuclearSegmentation
{
public:
	typedef struct { std::string name; int value; } Parameter;

	NuclearSegmentation();
	~NuclearSegmentation();

	//This is for beginning a completely new segmenation:
	bool LoadInput(std::string fname, int chNumber = 0);							//Load the input image from a file
	bool SetInput(ftk::Image::Pointer inImg, std::string fname, int chNumber = 0); //Pass a pointer to the already loaded image
	void SetParameters(std::string paramfile);
	void SetParameter(std::string name, int value);
	bool Binarize(bool getResultImg = false);					//Will binarize the data image
	bool DetectSeeds(bool getResultImg = false);				//If binarization has been done it will detect seeds
	bool RunClustering(bool getResultImg = false);				//Will use binary image and seeds to do initial clustering
	bool Finalize();											//Will finilize the output using alpha expansion
	bool SegmentAllTimes(bool finalize = false);				//Will segment the nuclei of all time points.
	void ReleaseSegMemory();									//Delete the NucleusSeg object to release all of its memory.
	bool ComputeAllGeometries();								//Compute all geometries for the label image!!!
	bool ComputeAllGeometries(int numTimes);					//Compute all geometries and return the features for the label image across time!!!

	//For use when loading from pre-existing results:
	bool LoadLabelImage(std::string fname);						//Filename of the label image to load
	bool LoadFromDAT(std::string fname);						//Load label image from .dat
	bool SetLabelImage(ftk::Image::Pointer labImg, std::string fname);	//Pass pointer to already loaded label image

	//Save functions:
	bool SaveLabelImage(std::string fname = "");				//Save the output image of the last step executed (image format)
	//bool SaveLabelByClass();									//Will save a different label image for each class
	
	//Editing Functions 
	std::vector< int > Split(ftk::Object::Point P1, ftk::Object::Point P2, vtkSmartPointer<vtkTable> table = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	std::vector< std::vector<int> > BatchSplit(std::vector<int> ids, int numObjs, vtkSmartPointer<vtkTable> table = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	std::vector< int > SplitAlongZ(int objID, int cutSlice, vtkSmartPointer<vtkTable> table = NULL);
	std::vector< std::vector< int > > GroupMerge(std::vector<int> ids, vtkSmartPointer<vtkTable> table = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	int Merge(std::vector<int> ids, vtkSmartPointer<vtkTable> table = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	bool Delete(std::vector<int> ids, vtkSmartPointer<vtkTable> table = NULL);
	bool Exclude(int l, int r, int t, int b, int z1, int z2, vtkSmartPointer<vtkTable> table = NULL);
	int AddObject(int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table = NULL);
	bool FillObjects(std::vector<int> ids);
	bool FillAnObject(int id);
	bool EditsNotSaved;
	vtkSmartPointer<vtkTable> megaTable;

	//Edits applied on initial segmentation and updates LoG resp image
	ftk::Object::Point MergeInit(ftk::Object::Point P1, ftk::Object::Point P2, int* new_id); 
	std::vector< int > SplitInit(ftk::Object::Point P1, ftk::Object::Point P2);
	bool DeleteInit(ftk::Object::Point P1);

	//Misc string Gets
	std::string GetErrorMessage() { return errorMessage; };
	std::string GetDataFilename() { return dataFilename; };
	std::string GetParamFilename() { return paramFilename; };
	std::string GetLabelFilename() { return labelFilename; };
	int GetParameter(std::string name);
	std::vector<std::string> GetParameterNames(){ return paramNames; };
	int GetNumberOfObjects(void){ return (int)centerMap.size(); };

	//Get Data:
	ftk::Image::Pointer GetDataImage(void){ return dataImage; };	
	ftk::Image::Pointer GetLabelImage(void){ return labelImage; };
	std::map<int, ftk::Object::Point> * GetCenterMapPointer(){ return &centerMap; };
	std::map<int, ftk::Object::Box> * GetBoundingBoxMapPointer(){ return &bBoxMap; };	
	std::vector<std::map<int, ftk::Object::Point> >  centerMap4DImage;
	std::vector<std::map<int, ftk::Object::Box> > bBoxMap4DImage;	
	std::vector<std::vector<ftk::IntrinsicFeatures> > featureVector4DImage;
	std::vector< vtkSmartPointer<vtkTable> >  table4DImage;
	//*********************************************************************************************

	//ADDED BY YOUSEF/RAGHAV:
	std::vector<std::string> RunGraphColoring(std::string labelname, std::string filename);	//Run Graph coloring on label image and save adjacency file as filename
	std::vector<Seed> getSeeds();
	void updatetable4DImage(std::vector< vtkSmartPointer<vtkTable> > tableOfFeatures);
	void SetCurrentbBox(std::map<int, ftk::Object::Box> currentbBoxMap); // To set the current bBox to store edit information when time series is loaded
	void createMegaTable();	
	void AddTimeToMegaTable();

	//Set Data:/Amin
#ifdef USE_TRACKING
	void SetTrackFeatures(std::vector<std::vector<ftk::TrackPointFeatures> > trackfeatures);
	void SetTimeFeatures(std::vector<ftk::TrackFeatures>);
#endif
	void setCurrentTime(int t = 0){currentTime = t;};
	void ReassignLabels(std::vector<int> times, std::vector<int> ids, std::vector<int> new_ids);

protected:
	std::string errorMessage;

	std::string dataFilename;			//Name of the file that the data came from
	ftk::Image::Pointer dataImage;		//The data image
	int channelNumber;					//Use this channel from the dataImage for segmentation
	std::string labelFilename;			//Name of the file that is the label image
	ftk::Image::Pointer labelImage;		//My label image
	yousef_nucleus_seg *NucleusSeg;		//The Nuclear Segmentation module
	int lastRunStep;					//0,1,2,3,4 for the stages in a nuclear segmentation.
	int currentTime;					//current Nucleus Editor time slider.
			
	std::string paramFilename;
	std::vector<Parameter> myParameters;
	std::vector<std::string> paramNames;

	//Geometry information that is kept for editing purposes:
	std::map<int, ftk::Object::Box>		bBoxMap;			//Bounding boxes
	std::map<int, ftk::Object::Point>	centerMap;			//Centroids


	bool GetResultImage();									//Gets the result of last module and puts it in labelImage
	void GetParameters(void);								//Retrieve the Parameters from nuclear segmentation.
	void ResetAll(void);									//Clear all memory and variables
	void ConvertParameters(int params[12]);
	
	//Editing Utilities:
	long int maxID(void);										//Get the maximum ID in the table!
	bool addObjectToMaps(int ID, int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table = NULL);
	bool addObjectsToMaps(std::set<unsigned short> IDs, int x1, int y1, int z1, int x2, int y2, int z2, vtkSmartPointer<vtkTable> table = NULL, vtkSmartPointer<vtkTable> NucAdjTable = NULL);
	void removeObjectFromMaps(int ID, vtkSmartPointer<vtkTable> table);
	//FeatureCalcType::Pointer computeGeometries(int x1, int y1, int z1, int x2, int y2, int z2);
	void ReassignLabels(std::vector<int> fromIds, int toId);
	void ReassignLabel(int fromId, int toId);
	ftk::Object::Box ExtremaBox(std::vector<int> ids);
	ftk::Object::Box GrowBox(ftk::Object::Box b, int s);
	std::vector<int> GetNeighbors(int id);
	vtkSmartPointer<vtkTable> featureVectorTovtkTable(std::vector<ftk::IntrinsicFeatures> featurevector);
	
	
	//FOR PRINTING SEEDS IMAGE:
	void Cleandptr(unsigned short*x, std::vector<int> y);
	void Restoredptr(unsigned short* );
	std::list<int> negativeseeds;


#ifdef USE_TRACKING
	//Add Tracking Features:
	std::vector<std::vector<ftk::TrackPointFeatures> > nucsegTrackFeatures;
	std::vector<ftk::TrackFeatures> nucsegTimeFeatures;
	std::vector<ftk::TrackPointFeatures> currentTrackFeatures; // means current time
	void setCurrentTrackFeatures(int time);
#endif

//********************************************************************************************
//********************************************************************************************

}; // end NuclearSegmentation

}  // end namespace ftk

#endif	// end __ftkNuclearSegmentation_h

