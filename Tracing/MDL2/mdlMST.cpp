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

#include "mdlMST.h"

namespace mdl
{

//Constructor
MST::MST(ImageType::Pointer inImage)
{
	m_inputImage = inImage;

	region = m_inputImage->GetBufferedRegion();
	sizeX = region.GetSize(0);
	sizeY = region.GetSize(1);
	sizeZ = region.GetSize(2);
	numPix = sizeX*sizeY*sizeZ;

	debug = false;
	edgeRange = 10;		//Some default values:
	power = 1;

	//input
	skeletonPoints = NULL;

	edge_array = NULL;
	edge_wght = NULL;
	num_edges = 0;
	g = NULL;

	//output
	nodes.clear();
	spanningTree.clear();
	nodeDegree.clear();
}

MST::~MST()
{
	m_inputImage=NULL;
	skeletonPoints=NULL;

	if(edge_array)
	{
		delete[] edge_array;
		edge_array = NULL;
	}
	if(edge_wght)
	{
		delete[] edge_wght;
		edge_wght = NULL;
	}
	if(g)
	{
		delete g;
		g = NULL;
	}
}

void MST::SetSkeletonPoints(std::vector<fPoint3D> * sp)
{
	skeletonPoints = sp;
}

bool MST::CreateGraphAndMST()
{
	if(!m_inputImage || !skeletonPoints)
		return false;

	this->skeletonPointsToNodes();
	this->nodesToEdges();
	this->minimumSpanningTree();

	return true;
}

int MST::roundToInt(double v)
{
	double intpart;

	if( modf(v, &intpart) < 0.5 )
		return (int)floor(v);
	else
		return (int)ceil(v);
}

//Because the skeleton Points may be float values and may have the same
// voxel that they are closest to.
//So I make sure that each voxel is in my list only once!!!
bool MST::skeletonPointsToNodes()
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
		node.x = roundToInt(fnode.x);
		node.y = roundToInt(fnode.y);
		node.z = roundToInt(fnode.z);
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

	int num_nodes = (int)nodes.size();

	if(debug)
		std::cerr << "Number of Nodes = " << num_nodes << std::endl;

	return true;
}

//I'm going to convert the nodes into edges and edge weights.
//I only make edges between nodes that are within the edgeRange
bool MST::nodesToEdges()
{
	if( (int)nodes.size() == 0)
		return false;

	std::vector<E> edgeArray;
	std::vector<float> edgeWeight;

	//int num_nodes = (int)nodes.size()-1;
	int num_nodes = (int)nodes.size();

	//Iterate through the nodes and create edges within the specified range.
	for(int i=0; i<num_nodes; ++i)
	{
		//for(int j=1; j<i; ++j)
		for(int j=0; j<i; ++j)
		{
			Point3D n1 = nodes.at(i);
			Point3D n2 = nodes.at(j);
			int dx = n1.x - n2.x;
			int dy = n1.y - n2.y;
			int dz = n1.z - n2.z;

			if( abs(dx) > edgeRange )
				continue;
			if( abs(dy) > edgeRange )
				continue;
			if( abs(dz) > edgeRange )
				continue;

			//If I'm here, then I've found a close enough node
			edgeArray.push_back( E(i+1,j+1) ); //add an edge (count starts at 1 for nodes)

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

	//Convert std c++ vector to c array:
	num_edges = (unsigned int)edgeArray.size();
	
	if(edge_array)
		delete[] edge_array;
	if(edge_wght)
		delete[] edge_wght;
	edge_array = new E[num_edges];
	edge_wght = new float[num_edges];

	for(unsigned int i = 0; i<num_edges; ++i)
	{
		edge_array[i] = edgeArray.at(i);
		edge_wght[i] = edgeWeight.at(i);
	}

	edgeWeight.clear();
	edgeArray.clear();

	if(debug)
		std::cerr << "Finished making edges = " << num_edges << std::endl;

	return true;
}

bool MST::minimumSpanningTree()
{
	if( num_edges <= 0 )
		return false;

	int num_nodes = (int)nodes.size();

	spanningTree.clear();

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	Graph * g = new Graph(num_nodes);
	for (size_t j = 0; j < num_edges; ++j) 
	{
		Edge e; 
		bool inserted;
		tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, *g);
		weightmap[e] = edge_w[j];
	}
#else
	if(g)
		delete g;
	g = new Graph(edge_array, edge_array + num_edges, edge_wght, num_nodes);
#endif

	if(!g)
		return false;

	boost::property_map<Graph, boost::edge_weight_t>::type weight = get(boost::edge_weight, *g);
	
	// MST algorithm
	boost::kruskal_minimum_spanning_tree(*g, back_inserter(spanningTree));

	if(debug)
		std::cerr << "kruskal_minimum_spanning_tree(MST) is finished!" << std::endl;

	// Create a graph for the initial MST
	Graph msTree(num_nodes+1);
	int num_edge_MST = 0;
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		add_edge(source(*ei, *g), target(*ei, *g), msTree);
		num_edge_MST++;
	}

	if(debug)
		std::cerr << "MST edges = " << num_edge_MST << std::endl;

	return true;
}

bool MST::ErodeAndDialateNodeDegree(int mophStrength)
{
	if((int)spanningTree.size() == 0 || !g)
	{
		std::cerr << "missing tree or g\n";
		return false;
	}

	int num_nodes = (int)nodes.size();

	int * degree_nodes = new int[num_nodes+1];
	int * degree_nodes_buffer = new int[num_nodes+1];
	if(!degree_nodes || !degree_nodes_buffer)
	{
		std::cerr << "Could not allocate degree buffers\n";
		return false;
	}

	//Inititialize for all the vertices.
	//Actual vertex index starts from 1
	for(int i=0; i<num_nodes+1; ++i)
	{
		degree_nodes[i] = 0;
		degree_nodes_buffer[i] = 0;
	}

	// create initial degree_nodes array
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
	{
		degree_nodes[source(*ei, *g)] ++;
		degree_nodes[target(*ei, *g)] ++;
	}

	// -- Erosion and Dilation of MST
	int times_erosion = mophStrength;

	int * edge_eroded = new int[num_nodes*2];
	if(!edge_eroded)
	{
		std::cerr << "Could not allocate edge_eroded\n";
		return false;
	}

	int num_edge_eroded = 0;
	while (times_erosion != 0) 
	{
		times_erosion--;
		for (int i=1; i<=num_nodes; i++)   
			degree_nodes_buffer[i] = degree_nodes[i];

		for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei) 
		{
			if (degree_nodes_buffer[source(*ei, *g)]>0 && degree_nodes_buffer[target(*ei, *g)]>0)  
			{
				if (degree_nodes_buffer[source(*ei, *g)]==1 || degree_nodes_buffer[target(*ei, *g)]==1)  
				{
					degree_nodes[source(*ei, *g)] --;
					degree_nodes[target(*ei, *g)] --;
          
					// Save the edges eroded in a stack-like array. Each edge takes two elements of the array
					edge_eroded[num_edge_eroded*2]  = (int)source(*ei, *g); // Saving of eroded edges
					edge_eroded[num_edge_eroded*2+1] = (int)target(*ei, *g);
					num_edge_eroded ++;
				}//end if	
			}// end if
		}// end for (vector...
	}// end while (times_erosion != 0)

	if(debug)
		std::cerr << "Erosion of MST is finished!" << std::endl;

	// Dilation the MST by counting up the degree of nodes
	while (num_edge_eroded !=0) 
	{
		num_edge_eroded --;
		int edge_source = edge_eroded[num_edge_eroded*2];  // Read the stored eroded edges
		int edge_target = edge_eroded[num_edge_eroded*2+1];
		if ((degree_nodes[edge_source]+ degree_nodes[edge_target]) == 1)  
		{  // if a branch tip edge
			degree_nodes[edge_source] ++;
			degree_nodes[edge_target] ++;
		}
	}

	if(debug)
		std::cerr << "Dilation of MST is finished!" << std::endl;

	//copy int output vectors:
	nodeDegree.clear();
	for(int i=0; i<num_nodes; ++i)
	{
		nodeDegree.push_back( degree_nodes[i+1] );
	}

	delete[] edge_eroded;
	delete[] degree_nodes_buffer;
	delete[] degree_nodes;
		
	if(debug)
		std::cerr << " DIALATED NODES = " << nodeDegree.size() << std::endl;

	if(debug)
	{
		FILE * fout = fopen("degrees.txt", "w");
		if (fout != NULL)
		{
			for(int i=0; i<(int)nodeDegree.size(); ++i)
			{
				fprintf(fout,"%d %f\n", i, (float)nodeDegree.at(i));
			}
			fclose(fout); 
		}
	}
	return true;
}

std::vector<MST::E> MST::BackboneExtract()
{
	std::vector<E> retLines;

	if((int)spanningTree.size() == 0 || !g || (int)nodeDegree.size() == 0)
	{
		std::cerr << "missing tree or g or node degrees\n";
		return retLines;
	}

	//Create a msTree graph for backbone from the MST generated above
	//BackBone graph created or Graph for the parts of total skeletons 
	Graph msTreeBB(1);
	for (std::vector < Edge >::iterator ei = spanningTree.begin(); ei != spanningTree.end(); ++ei)
	{
		if (nodeDegree.at(source(*ei, *g)-1)!=0 && nodeDegree.at(target(*ei, *g)-1)!=0)
		{  
			add_edge(source(*ei, *g), target(*ei, *g), msTreeBB);
		}
    }

	//Get the lines from the backbone mst:
	typedef boost::graph_traits < Graph >::edge_iterator Edge_iter;
	Edge_iter ei, ei_end; 
	for (tie(ei, ei_end) = edges(msTreeBB); ei != ei_end; ++ei)
	{
		E ne( (int)source(*ei, msTreeBB)-1, (int)target(*ei, msTreeBB)-1 );
		retLines.push_back( ne );
   }

	if(debug)
		std::cerr << "Number of backbone lines = " << retLines.size() << std::endl;

	if(debug)
	{
		FILE * fout = fopen("BackboneCandidate.vtk", "w");
		if (fout != NULL)
		{

			int num_nodes = (int)nodes.size();
			int num_lines = (int)retLines.size();

			fprintf(fout, "# vtk DataFile Version 3.0\n");
			fprintf(fout,"MST of skel\n");
			fprintf(fout,"ASCII\n");
			fprintf(fout,"DATASET POLYDATA\n");
			fprintf(fout,"POINTS %d float\n",num_nodes);

			for(int i=0; i<num_nodes; ++i)
			{
				mdl::Point3D nd = nodes.at(i);
				fprintf(fout,"%f %f %f\n", (float)nd.x, (float)nd.y, (float)nd.z);
			}

			fprintf(fout,"LINES %d %d\n", num_lines, num_lines*3);

			for(int i=0; i<num_lines; ++i)
			{
				E e = retLines.at(i);
				fprintf(fout, "2 %d %d\n", e.first, e.second);
			}

			fclose(fout);
		}
	}

	return retLines;
}

}