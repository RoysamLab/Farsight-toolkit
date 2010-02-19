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

// BackboneExtract - Extract dendrite bacbone from a list of 3D points
// Input format:  .skel  (3D points)
// Output format: .vtk  (3D graph format)
// Author: Xiao liang
// Date: 12/11/2009

#include "MST.h"

int main(int argc, char *argv[])
{
  ifstream fin;
  FILE *volfile;
  FILE *fout;
  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int i;
  int ii, jj, kk;
  long idx, iidx;
  long slsz, sz;
  int p;
  int edgeRange;
  int line_count;
  int *voxelNodeIndex;  // dimension of sizeX * sizeY * sizeZ, for nodes index
  DATATYPEIN *volin;
  VoxelPosition nodePosition;
  VoxelPosition *vertexPos;
  int num_nodes;
  int idx_edge;
  

  E *edge_array;
  float *edge_w;
  int *degree_nodes;  // record degree for each node
  int *degree_nodes_buffer;
  int times_erosion;
  int times_dilation;
  float densityFactor;
  int *edge_eroded;
  int num_edge_eroded;
  int edge_source, edge_target;
  if (argc < 10)
    {
    cout << argv[0] << " <data dir> <skel pt file> <vol file> <xs> <ys> <zs> "
         << "<edgeRange> <Erosion seps>"
         << "<out backbone VTK file >" << " Power " <<endl;
    return 1;
    }
  
  //- open skeleton points file 
  std::string infilename = argv[1];
  infilename += argv[2];
  fin.open(infilename.c_str());
  if (!fin)
    {
    cerr << "couldn't open skel file " << infilename << " for input" << endl;
    return -1;
    }

  //- open intensty image file 
  std::string volfilename = argv[1];
  volfilename += argv[3];
  if( (volfile=fopen(volfilename.c_str(), "rb") ) == NULL)  // open vol file
    {
    cerr << "couldn't open volfile " << volfilename << " for input" << endl;
    return -1;
    }
  
  //- get the size of intensity images 
  sizeX = atoi(argv[4]);
  sizeY = atoi(argv[5]);
  sizeZ = atoi(argv[6]);

  //- Search Range
  edgeRange = atoi(argv[7]);
  times_erosion = atoi(argv[8]);

  if ((fout = fopen(argv[9], "w")) == NULL)
    {
    cerr << "Cannot open " << argv[9] << " for writing" << endl;
    return 1;
    }
  double pr=atof(argv[10]); 

  voxelNodeIndex = new int[sizeX*sizeY*sizeZ];
  int *nodeIndicesInitialized = new int[sizeX*sizeY*sizeZ];
  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));

  // -slice size
  slsz = sizeX*sizeY;   
  sz = slsz*sizeZ;
  edge_array = new E[MAX_NUM_EDGE];
  edge_w = new float[MAX_NUM_EDGE];

  // -read in vol file
  if ( fread(volin, sizeof(DATATYPEIN), sz, volfile) < (unsigned long)sz)  
  {
    cerr << "File size is not the same as volume size" << endl;
    return 1;
  }
  fclose(volfile);


  //- Initialize to zero
  for(idx=0; idx<sz; idx++)   { 
    voxelNodeIndex[idx]=0;
    nodeIndicesInitialized[idx]=0;

  }
   
  num_nodes = 0;  // initial
 //int Allnumnodes =0; 
 int itr = 0;
 fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
  while (!fin.eof() ) 
    {
    itr++;
    idx = int(nodePosition.z)*slsz + int(nodePosition.y)*sizeX + int(nodePosition.x);
    if (nodeIndicesInitialized[idx] == 0)
      {
      nodeIndicesInitialized[idx] = 1;
      num_nodes++;
      }
    fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
    }
 
  vertexPos = new VoxelPosition[num_nodes+1];
  //print header information in the main output file

  fprintf(fout, "# vtk DataFile Version 3.0\n");
  fprintf(fout,"MST of skel\n");
  fprintf(fout,"ASCII\n");
  fprintf(fout,"DATASET POLYDATA\n");
  fprintf(fout,"POINTS %d float\n",num_nodes);


  //reinitialize the file and variables used to loop through it

  fin.clear();
  fin.seekg(0, ios::beg);
  num_nodes = 0;  // initial
  idx_edge = 0;   // initial
  itr = 0;
  fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
  cout << "searching input file for neighboring nodes..." << endl;
  while (!fin.eof() ) 
    {
    itr++;
    idx = int(nodePosition.z)*slsz + int(nodePosition.y)*sizeX + int(nodePosition.x);
    if (voxelNodeIndex[idx] == 0)
      {
      num_nodes++;
      // Save the index of node (vertex) at each voxel position
      voxelNodeIndex[idx] = num_nodes;
      // output the node positions to vtk file
      fprintf(fout,"%f %f %f\n", nodePosition.x, nodePosition.y, nodePosition.z);  
	  vertexPos[voxelNodeIndex[idx]].x =nodePosition.x;
	  vertexPos[voxelNodeIndex[idx]].y =nodePosition.y;
	  vertexPos[voxelNodeIndex[idx]].z =nodePosition.z;

      // Find all neighbor nodes within edgeRange
      for (kk = -edgeRange; kk <= edgeRange; kk++)
        {
        for (jj = -edgeRange; jj <= edgeRange; jj++)
          {
          for (ii = -edgeRange; ii <= edgeRange; ii++)
            {
            if (ii==0 && jj==0 && kk==0)
              {
              //don't consider current point
              continue;
              }
            if (int(nodePosition.x)+ii < 0 || int(nodePosition.x)+ii >= sizeX)
              {
              //check if the block is inside the image
              continue;
              }
            if (int(nodePosition.y)+jj < 0 || int(nodePosition.y)+jj >= sizeY)
              {
              continue;
              }
            if (int(nodePosition.z)+kk < 0 || int(nodePosition.z)+kk >= sizeZ)
              {
              continue;
              }
            iidx = kk * slsz + jj*sizeX + ii;
            int iidxMid1 = int(kk/3.0) * slsz + int(jj/3.0)*sizeX + int(ii/3.0);
            int iidxMid2 = int(kk*2.0/3.0) * slsz + int(jj*2.0/3.0)*sizeX + int(ii*2.0/3.0);
            iidx = idx + iidx;
            iidxMid1 = idx + iidxMid1;
            iidxMid2 = idx + iidxMid2;
            if (iidx<0 || iidx >= sz)
              {
              continue;
              }
            if (voxelNodeIndex[iidx] != 0)
              {
              // put into edge array of MST
              edge_array[idx_edge] = E(num_nodes, voxelNodeIndex[iidx]);
			  densityFactor = 0;
              if(pr<0.1) // when the input is very small 
			  {
			  densityFactor =0; 
			  }
			  else 
			  {
              densityFactor =(float)  pow(double(volin[idx]+0.001), pr); 
              densityFactor += (float)pow(double(volin[iidx]+0.001),pr);      
			  densityFactor =(float) fabs(pow(double(densityFactor),double(1.05)) +pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5));   // test by square densityFactor
			  }
       
              edge_w[idx_edge] = (float)(sqrt(float(kk*kk + jj*jj + ii*ii)) / (densityFactor*0.02+1));
              idx_edge++;
              }
            }
          }
        }
      }
    fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
    }

   cout << "finished reading nodes" << endl; 
   degree_nodes = new int[num_nodes+1];
   degree_nodes_buffer = new int[num_nodes+1];


  // Initialize for all the vertices. Actual vertex index starts from 1.
  for (i=0; i<=num_nodes; i++)
    {
    degree_nodes[i] = 0;  //initialize to zeros
    }

  size_t num_edges = idx_edge;   // sizeof(edge_array) / sizeof(E);
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
  Graph g(num_nodes);
  property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
  for (size_t j = 0; j < num_edges; ++j) {
    Edge e; bool inserted;
    tie(e, inserted) = add_edge(edge_array[j].first, edge_array[j].second, g);
    weightmap[e] = edge_w[j];
  }
#else
  Graph g(edge_array, edge_array + num_edges, edge_w, num_nodes);
#endif

  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
  vector < Edge > spanning_tree;

  // MST algorithm
  kruskal_minimum_spanning_tree(g, back_inserter(spanning_tree));

  cout << "kruskal_minimum_spanning_tree(MST) is finished!" << endl;
  
  // create initial degree_nodes array
  for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    degree_nodes[source(*ei, g)] ++;
    degree_nodes[target(*ei, g)] ++;
  }

  // Create a graph for the initial MST
  Graph msTree(num_nodes+1);
  int num_edge_MST = 0;
  for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
      //if (degree_nodes_initialMST[source(*ei, g)]!=0 && degree_nodes_initialMST[target(*ei, g)]!=0) 
    { 
          add_edge(source(*ei, g), target(*ei, g), msTree);
          num_edge_MST++;
    }
  }

  // -- Erosion and Dilation of MST
  times_dilation = times_erosion;
  edge_eroded = new int[num_nodes*2];
  num_edge_eroded = 0;
  while (times_erosion !=0) {
    times_erosion--;
    for (i=1; i<=num_nodes; i++)   degree_nodes_buffer[i] = degree_nodes[i];
    for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
      if (degree_nodes_buffer[source(*ei, g)]>0 && degree_nodes_buffer[target(*ei, g)]>0)  {
              if (degree_nodes_buffer[source(*ei, g)]==1 || degree_nodes_buffer[target(*ei, g)]==1)  {
                  degree_nodes[source(*ei, g)] --;
                  degree_nodes[target(*ei, g)] --;
          // Save the edges eroded in a stack-like array. Each edge takes two elements of the array
          edge_eroded[num_edge_eroded*2]  = source(*ei, g); // Saving of eroded edges
          edge_eroded[num_edge_eroded*2+1] = target(*ei, g);
          num_edge_eroded ++;
              }
      }
    }
  }
  cout << "Erosion of MST is finished!" << endl;
  // Dilation the MST by counting up the degree of nodes
  while (num_edge_eroded !=0) {
    //times_dilation--;
    num_edge_eroded --;
      //for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    edge_source = edge_eroded[num_edge_eroded*2];  // Read the stored eroded edges
    edge_target = edge_eroded[num_edge_eroded*2+1];
      if ((degree_nodes[edge_source]+ degree_nodes[edge_target]) == 1)  {  // if a branch tip edge
           degree_nodes[edge_source] ++;
           degree_nodes[edge_target] ++;
    }
  }

  cout << "Dilation of MST is finished!" << endl;

  num_nodes                 =0;
  // Create a msTree graph for backbone from the MST generated above
  Graph msTreeBB(num_nodes+1);        // BackBone graph created or Graph for the parts of total skeletons 
  for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei)
    {
    if (degree_nodes[source(*ei, g)]!=0 && degree_nodes[target(*ei, g)]!=0)
      {  
      add_edge(source(*ei, g), target(*ei, g), msTreeBB);
	  num_nodes++;
      }
    }
 
 cout << "Lines generating is finished!" << endl;

 Edge_iter   ei, ei_end; 
 line_count = 0;

 for (tie(ei, ei_end) = edges(msTreeBB); ei != ei_end; ++ei)
   {
   line_count++; // count the number of lines output in vtk file
   }
 fprintf(fout,"LINES %d %d\n", line_count, line_count*3);

 for (tie(ei, ei_end) = edges(msTreeBB); ei != ei_end; ++ei)
   {
   // output lines into vtk file
   //fprintf(fout, "2 %zu %zu\n", source(*ei, msTreeBB)-1, target(*ei, msTreeBB)-1);
	   fprintf(fout, "2 %ld %ld\n", source(*ei, msTreeBB)-1, target(*ei, msTreeBB)-1);
   }
  
 cout << "Sucess!" << endl;

  fclose(fout); 

 
  //cleanup the memory we've allocated...
  free(volin); 
  volin = NULL;
  delete [] nodeIndicesInitialized;
  delete [] voxelNodeIndex;
  delete [] edge_array;
  delete [] edge_w;
  delete [] vertexPos;
  delete [] degree_nodes;
  delete [] degree_nodes_buffer;
  delete [] edge_eroded;

  cout << "-- Number of backbone nodes is : " << num_nodes << endl;
  cout << "-- Number of backbone Lines is : " << line_count << endl;

  return EXIT_SUCCESS;

}

