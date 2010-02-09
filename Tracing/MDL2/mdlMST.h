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

//#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <stdio.h>
//#include <algorithm>

#include "itkImage.h"
//#include "itkImageRegionIterator.h"
//#include "itkImageRegionIteratorWithIndex.h"

namespace mdl
{

class MST
{
public:
	typedef std::pair<int, int>  E;

	MST(ImageType::Pointer inImage);
	~MST();
	//Setup:
	void SetDebug(bool inp = true){ debug = inp; };
	void SetEdgeRange(int edge){ edgeRange = edge; };
	void SetPower(int p){ power = p; };
	//Methods:
	void SetSkeletonPoints(std::vector<fPoint3D> * sp);
	bool CreateGraphAndMST();	//Do first
	bool ErodeAndDialateNodeDegree(int morphStrength); //Do second
	std::vector<E> BackboneExtract();
	


	//Get Result:

private:
	//Parameters
	bool debug;				//If debug is true, process in steps and print stuff
	int edgeRange;
	double power;

	//Images & size
	ImageType::Pointer m_inputImage;
	ImageType::RegionType region;
	int sizeX;
	int sizeY;
	int sizeZ;
	long numPix;

	//Skeleton points (input)
	std::vector<fPoint3D> * skeletonPoints;

	//Intermediates:
	E * edge_array;			//for initial edges
	float * edge_wght;		//for initial edge weights
	unsigned int num_edges; //number of initial edges
	Graph * g;				//my full graph at beginning

	bool skeletonPointsToNodes();	//step 1
	bool nodesToEdges();			//step 2
	bool minimumSpanningTree();		//step 3
	int roundToInt(double v);

	//output:
	//MST result
	std::vector<Point3D> nodes;
	std::vector< Edge > spanningTree;
	std::vector< int > nodeDegree;
};

}  // end namespace mdl

#endif
