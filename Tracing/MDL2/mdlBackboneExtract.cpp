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
#if defined(_MSC_VER)
#pragma warning(disable : 4996)
#endif

#include "mdlBackboneExtract.h"

namespace mdl
{

//Constructor
BackboneExtract::BackboneExtract(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	region = m_inputImage->GetBufferedRegion();
	sizeX = region.GetSize(0);
	sizeY = region.GetSize(1);
	sizeZ = region.GetSize(2);
	numPix = sizeX*sizeY*sizeZ;

	debug = false;
	edgeRange = 10;		//Some default values:
	timesErosion = 50;
	power = 1;

	nodes.clear();
	edgeArray.clear();
	edgeWeight.clear();

	skeletonPoints = NULL;
}

BackboneExtract::~BackboneExtract()
{
	m_inputImage=NULL;
	skeletonPoints=NULL;

	nodes.clear();
	edgeArray.clear();
	edgeWeight.clear();
}

void BackboneExtract::SetSkeletonPoints(std::vector<fPoint3D> * sp)
{
	skeletonPoints = sp;

	nodes.clear();
	edgeArray.clear();
	edgeWeight.clear();
}

bool BackboneExtract::Update()
{
	if(!m_inputImage || !skeletonPoints)
		return false;

	this->skeletonPointsToNodes();
	this->nodesToEdges();

	return true;
}

//Because the skeleton Points may be float values and may have the same
// voxel that they are closest to.
//So I make sure that each voxel is in my list only once!!!
bool BackboneExtract::skeletonPointsToNodes()
{
	if(!skeletonPoints)
		return false;

	bool * nodeAdded = new bool[numPix];
	if(!nodeAdded)	//Couldn't allocate memory
		return false;

	//- Initialize to zero
	for(long idx=0; idx<numPix; idx++)   
	{ 
		nodeAdded[idx] = false;
	}

	nodes.clear();

	for(int i=0; i<(int)skeletonPoints->size(); ++i)
	{
		fPoint3D fnode = skeletonPoints->at(i);
		Point3D node;
		node.x = (int)fnode.x;
		node.y = (int)fnode.y;
		node.z = (int)fnode.z;
		long idx = (node.z)*sizeX*sizeY + (node.y)*sizeX + (node.x);
		if(idx >= numPix)
			continue;

		if(!nodeAdded[idx])
		{
			nodeAdded[idx] = true;
			nodes.push_back(node);
		}
	}

	delete[] nodeAdded;

	return true;
}

//I'm going to convert the nodes into edges and edge weights.
//I only make edges between nodes that are within the edgeRange
bool BackboneExtract::nodesToEdges()
{
	if( (int)nodes.size() == 0)
		return false;

	edgeArray.clear();
	edgeWeight.clear();

	//Iterate through the nodes and create edges within the specified range.
	for(int i=0; i<nodes.size(); ++i)
	{
		for(int j=i+1; j<nodes.size(); ++j)
		{
			Point3D n1 = nodes.at(i);
			Point3D n2 = nodes.at(j);
			int dx = abs(n1.x - n2.x);
			int dy = abs(n1.y - n2.y);
			int dz = abs(n1.z - n2.z);

			if( dx > edgeRange )
				continue;
			if( dy > edgeRange )
				continue;
			if( dz > edgeRange )
				continue;

			//If I'm here, then I've found a close enough node
			edgeArray.push_back( E(i,j) ); //add an edge

			//Now compute edge weight:
			float densityFactor = 0;
			if(power >= 0.1)
			{
				//Get the two intensity values:
				ImageType::IndexType index1;
				index1[0] = n1.x;
				index1[1] = n1.y;
				index1[2] = n1.z;
				ImageType::IndexType index2;
				index2[0] = n2.x;
				index2[1] = n2.y;
				index2[2] = n2.z;
				PixelType pix1 = m_inputImage->GetPixel(index1);
				PixelType pix2 = m_inputImage->GetPixel(index2);

				//Adjust the intensity values
				float i1adj = (float)pow(double(pix1)+0.001, power);
				float i2adj = (float)pow(double(pix2)+0.001, power);

				//Compute the densityFactor:
				float term1 = (float)pow(double(i1adj+i2adj),double(1.05));
				// AAAHH THIS IS GROSS AND DOES NOT MAKE SENSE
				//float term2 = pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5));
				float term2 = 0;
				densityFactor = fabs(term1 + term2);
			}

			float num = (float)sqrt(float(dx*dx + dy*dy + dz*dz));
			float den = densityFactor*.02 + 1;
			float weight = num / den;
			edgeWeight.push_back(weight);
		}//end for j
	}//end for i
}


/*
//Search skeletonPoints for neighboring nodes
//Create 3 new vectors:
// 1. vertices
// 2. edge_array
// 3. edge_w
bool BackboneExtract::readNodes()
{
	if(!skeletonPoints)
		return false;

	if(debug)
		std::cerr << "Searching for neighboring nodes..." << std::endl;

	//Create voxelNodeIndex buffer for keeping track of indices.
	int * voxelNodeIndex = new int[numPix];
	if(!voxelNodeIndex)
	{
		//Could not allocate memory
		return false;
	}
	else
	{
		//- Initialize to zero
		for(int idx=0; idx<numPix; idx++)   
			voxelNodeIndex[idx]=0;
	}

	int slsz = sizeX*sizeY;		//Pixels in a slice

	int num_nodes = 0;  // initial
	float densityFactor = 0;

	std::vector<Point3D> vertexPoints;
	std::vector<float> edge_w;	//edge weights
	std::vector<E> edge_array;	//edge pairs

	for(int itr=0; itr<skeletonPoints->size(); ++itr)
	{
		Point3D nPos = skeletonPoints->at(itr);	//nodePosition
		long idx = int(nPos.z)*slsz + int(nPos.y)*sizeX + int(nPos.x);
	
		if(voxelNodeIndex[idx] == 0)
		{
			num_nodes++;
			voxelNodeIndex[idx] = num_nodes;

			vertexPoints.push_back( nPos );

			// Find all neighbor nodes within edgeRange
			for (int kk = -edgeRange; kk <= edgeRange; kk++)
			{
				for (int jj = -edgeRange; jj <= edgeRange; jj++)
				{
					for (int ii = -edgeRange; ii <= edgeRange; ii++)
					{
						//don't consider current point
						if (ii==0 && jj==0 && kk==0)
							continue;

						//check if the block is inside the image
						if (int(nPos.x)+ii < 0 || int(nPos.x)+ii >= sizeX)
							continue;
						if (int(nPos.y)+jj < 0 || int(nPos.y)+jj >= sizeY)
							continue;
						if (int(nPos.z)+kk < 0 || int(nPos.z)+kk >= sizeZ)
							continue;

						long iidx = kk*slsz + jj*sizeX + ii;
						int iidxMid1 = int(kk/3.0) * slsz + int(jj/3.0)*sizeX + int(ii/3.0);
						int iidxMid2 = int(kk*2.0/3.0)*slsz + int(jj*2.0/3.0)*sizeX + int(ii*2.0/3.0);
						iidx = idx + iidx;
						iidxMid1 = idx + iidxMid1;
						iidxMid2 = idx + iidxMid2;
						if (iidx<0 || iidx >= numPix)
							continue;

						if(voxelNodeIndex[iidx] != 0)
						{
							edge_array.push_back( E(num_nodes, voxelNodeIndex[iidx]) );

							densityFactor = 0;
							if(power >= 0.1)
							{
								// test by square densityFactor
								PixelType pix = m_inputImage->GetPixel((int)nPos.x,(int)nPos.y,(int)nPos.z);
								PixelType pix2 = m_inputImage->GetPixel(int(nPos.x)+ii,int(nPos.y)+jj,int(nPos.z)+kk);
								float temp1 = (float)pow(double(pix+0.001), power);
								temp1 += (float)pow(double(pix2+0.001),power);
								double temp2 = pow(double(temp1), double(1.05));
								double temp3 = pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5);
								densityFactor =(float)fabs(temp2 + temp3);
							}
							
							float div = 1.0 + (densityFactor*0.02);
							float dist = sqrt(float(kk*kk + jj*jj + ii*ii));
							float weight = (float)(dist / div);
							edge_w.push_back(weight);
						}//end if voxelNodeIndex
					}//end for ii
				}//end for jj
			}//end for kk

		}//end if(voxelNodeIndex[idx])
	}//end for itr

	return true;
}
*/

}