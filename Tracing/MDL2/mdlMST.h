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
/**************************************************************************  
// BackboneExtract - Extract dendrite backbone from a list of 3D points
// Input format:  std::vector of 3D points
// Output format: vtk graph (.vtk 3D graph format)
// Author: Xiao liang, RPI  based on Xiaosong's MST 
// Date: 12/11/2009
 *  Adapted Jan. 2010 by Isaac Abbott 
 *        
 *************************************************************************/
#ifndef __mdlMST_h
#define __mdlMST_h

#include "mdlTypes.h"
#include "mdlUtils.h"
#include "WeightedMahalsnobisDistance.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "itkImage.h"


namespace mdl
{

class MST
{
public:
	MST(ImageType::Pointer inImage);
	
	~MST();
	//Setup:
	void SetDebug(bool inp = true){ debug = inp; };
	void SetUseVoxelRounding(bool inp = true){useVoxelRounding = inp;};
	void SetEdgeRange(int edge){ edgeRange = edge; };
	void SetAlpha(double alpha){Alpha = alpha;}	//For spines
	void SetPower(int p){ power = p; };
	void SetPruneThreshold(double p){PruneThreshold = p;}

	//Added by Xiao Liang:
	void SetVesselMap(ImageType::Pointer VesselMap);
	void SetFileofRealSpineFeature(char *FileofRealSpineFeature){RealSpineFeatureFilename = FileofRealSpineFeature;}
	void SetFileofNonSpineFeature(char *FileofNonSpineFeature){NonSpineFeatureFilename = FileofNonSpineFeature;}

	//void SetInputforSpineExtraction(ImageType::Pointer VesselMap,char *FileofRealSpineFeature,char *FileofNonSpineFeature);
	
	//Methods:
	void SetSkeletonPoints(std::vector<fPoint3D> * sp);
	bool CreateGraphAndMST(int type = 1 ); //Do first
	bool ErodeAndDialateNodeDegree(int morphStrength); //Do second (if desired)


	std::vector<pairE> SearchFirstandSecondLevelBranch(void);

	//The number is the edge are the node numbers (starting at 1)
	std::vector<pairE> BackboneExtract();
	std::vector<pairE> SpineExtract();		//MDL based Spine Extraction

	//Get Result:
	std::vector<fPoint3D> GetNodes(){ return nodes; };

private:
	//Parameters
	bool   debug;				//If debug is true, process in steps and print stuff
	bool   useVoxelRounding;  //Round Nodes to nearest integer (or voxel)
	int    edgeRange;
	double power;
	double Alpha;
	double PruneThreshold;

	//Images & size
	ImageType::Pointer m_inputImage; // for the intesity image
	ImageType::Pointer m_VesselMap; // for the VesselMap

	ImageType::RegionType region;
	int sizeX;
	int sizeY;
	int sizeZ;
	long numPix;

	//Skeleton points (input)
	std::vector<fPoint3D> * skeletonPoints;

	//Feature 

	char* RealSpineFeatureFilename;
    char* NonSpineFeatureFilename;

	//Intermediates:
	std::vector<pairE> edgeArray;	//for initial edges
	std::vector<float> edgeWeight;	//for initial edge weights

	Graph * nodeGraph;			//my full graph of all nodes
	Graph * mstGraph;			//min spanning tree graph

	bool skeletonPointsToNodes(bool roundToNearestVoxel=true);	//step 1
	bool nodesToEdges(int type = 1);//step 2
	float getXiaosongEdgeWeight(fPoint3D n1, fPoint3D n2, double pwr);
	float getEdgeWeight(fPoint3D n1, fPoint3D n2, ImageType::Pointer img);
	float getGeodesicEdgeWeight(fPoint3D n1, fPoint3D n2, ImageType::Pointer img);
	bool minimumSpanningTree();		//step 3
	Graph morphGraphPrune(Graph *graph, std::vector<fPoint3D> *nodes, float lengthThreshold);
	int roundToInt(double v);

	//output:
	//MST result
	std::vector<fPoint3D> nodes;
	std::vector< Edge > spanningTree;
	std::vector< int > nodeDegree;
};

}  // end namespace mdl

#endif
