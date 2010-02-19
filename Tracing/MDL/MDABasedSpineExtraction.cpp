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
// MinSpanTree - Generate graph structure from a list of 3D points
// Input format: .skel  (3D points)
// Output format: .vtk  (3D graph format)
// Author: Xiaosong Yuan and Xiao liang, RPI
// Date: 19/12/2009
// Status: under modification of MDL

#include "MST.h"
#include "WeightedMahalsnobisDistance.h"

#define InterMedial 0
#define FeatureNumber  3;

double mahalanobisDist(double meanDensityBranch, double length_leaf, double meanVesselBranch, int spineOne);

int main(int argc, char *argv[])
{
  ifstream fin;
  FILE *volfile, *vesselfile;
  DATATYPEIN *volin, *volvessel;
  FILE *fout_MDL;  // add by xiao
  FILE *fout_Spine; // add by xiao
  int sizeX,sizeY,sizeZ;         // Sizes in x,y,z dimensions
  int i, j, k;
  int ii, jj, kk;
  long idx, iidx;
  long slsz, sz;
  int p;
  int edgeRange;
  int line_count;
  int branchChosen;
  int indVert, indVert_last;

  int *voxelNodeIndex;  // dimension of sizeX * sizeY * sizeZ, for nodes index
  VoxelPosition nodePosition;
  VoxelPosition *vertexPos;
  int num_nodes;
  int idx_edge;
  int num_leaves;

  E *edge_array;
  float *edge_w;
  int *degree_nodes;  // record degree for each node
  int *degree_nodes_tree; //record degree for final node
  int *degree_nodes_buffer;
  int *degree_nodes_initialMST;
  int times_erosion = 0;
  int times_dilation;
  float densityFactor;
  double meanDensityBranch[MAXNumBranch];   // Suppose at most MAXNumBranch branches at the 2nd level branch from BB
  double meanVesselBranch[MAXNumBranch];
  double length_leaf[MAXNumBranch];

  double aveDensityBranch[MAXNumBranch];   // Suppose at most MAXNumBranch branches at the 2nd level branch from BB
  double aveVesselBranch[MAXNumBranch];
  double length_2leaf[MAXNumBranch];  // length of two level branches

  float length_edge, leaf_length;
  int *edge_eroded;
  int num_edge_eroded;
  int edge_source, edge_target;
  int *vertBackbone;
  int numBranch_on_Backbone;
  int index_vert;
  double mahalanobis_dist[MAXNumBranch];
  double mahalanobis_dist_nonSpine[MAXNumBranch];
  double mahalanobis_dist_min;
  int mahalanobis_dist_minIndex;
  double MDL;
  double MDL_min;
  int MDL_minIndex;
  double sum_mahalanobis_nonSpine;
  double alpha;

  FILE *fclass_identify;  // this file is used to record the spine candidate' possible feature; 

  FILE *foutSpineCandidate; // this file is used to record the spine candidate;

  if (argc < 14)
    {
    cerr << argv[0] << " <data dir> <skel pt file> <vol file> <xs> <ys> <zs> "
       << "<edgeRange> <graph prune size> <morph strength> <weight factor> "
       << "<full path to vessel file> <out Feature txt file> "
       << "<out spine VTK graph>" << endl;
    return 1;
    }

  std::string infilename = argv[1];
  infilename += argv[2];
  fin.open(infilename.c_str());
  if (!fin)
    {
    cerr << "couldn't open skel file " << infilename << " for input" << endl;
    return -1;
    }

  std::string volfilename = argv[1];
  volfilename += argv[3];

  if( (volfile=fopen(volfilename.c_str(), "rb") ) == NULL)  // open vol file
    {
    cerr << "couldn't open volfile " << volfilename << " for input" << endl;
    return -1;
    }
  
  sizeX = atoi(argv[4]);
  sizeY = atoi(argv[5]);
  sizeZ = atoi(argv[6]);
  edgeRange = atoi(argv[7]);
  leaf_length =(float) atof(argv[8]);
  times_erosion = atoi(argv[9]);
  alpha = atof(argv[10]);

  if((vesselfile=fopen(argv[11], "rb"))==NULL)  // open vol file
    {
    cerr << "couldn't open vessel file for input" << endl;
    return -1;
    }

  if ((fout_MDL = fopen(argv[12], "w")) == NULL) // 
    {
    cerr << "Cannot open " << argv[12] << " for writing" << endl;
    return 1;
    }
  
  if ((fout_Spine = fopen(argv[13], "w")) == NULL)
    {
    cerr << "Cannot open " << argv[13] << " for writing" << endl;
    return 1;
    }
  
  if ((fclass_identify = fopen("CLASSIFIER_TRAINING.txt", "w")) == NULL)  
    {
    cerr << "Cannot open CLASSIFIER_TRAINING.txt for writing" << endl;
    return 1;
    }

  if ((foutSpineCandidate = fopen("Spine_Candiate_MedialResult.vtk ", "w")) == NULL)  
    {
    cerr << "Cannot open Spine_Candiate_MedialResult.vtk for writing" << endl;
    return 1;
    }

  #if InterMedial
    FILE *tempfile1;
    if ((tempfile1 = fopen("MST.vtk ", "w")) == NULL)  
     {
     cerr << "Cannot open MST.vtk for writing" << endl;
     return 1;
     }
  #endif

  voxelNodeIndex = new int[sizeX*sizeY*sizeZ];
  int *nodeIndicesInitialized = new int[sizeX*sizeY*sizeZ];
  volin = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  volvessel = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  slsz = sizeX*sizeY;   // slice size
  sz = slsz*sizeZ;
  edge_array = new E[MAX_NUM_EDGE];
  edge_w = new float[MAX_NUM_EDGE];
  //edge_array = (E*)malloc(MAX_NUM_EDGE*sizeof(E));

  if ( fread(volin, sizeof(DATATYPEIN), sz, volfile) < (unsigned long)sz)  // read in vol file
    {
    cerr << "File size is not the same as volume size" << endl;
    return 1;
    }

  if ( fread(volvessel, sizeof(DATATYPEIN), sz, vesselfile) < (unsigned long)sz)  // read in vessel file
    {
    cerr << "File size is not the same as vessel size" << endl;
    return 1;
    }

  fclose(volfile);
  fclose(vesselfile);


  for(idx=0; idx<sz; idx++)   {  //Initialize to zero
    voxelNodeIndex[idx]=0;
    nodeIndicesInitialized[idx]=0;

  }

  num_nodes = 0;  // initial
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

  //print header information in the main output file
 
  fprintf(fout_Spine, "# vtk DataFile Version 3.0\n");
  fprintf(fout_Spine,"MST of skel\n");
  fprintf(fout_Spine,"ASCII\n");
  fprintf(fout_Spine,"DATASET POLYDATA\n");
  fprintf(fout_Spine,"POINTS %d float\n",num_nodes);

  fprintf(foutSpineCandidate,"# vtk DataFile Version 3.0\n");
  fprintf(foutSpineCandidate,"MST of skel\n");
  fprintf(foutSpineCandidate,"ASCII\n");
  fprintf(foutSpineCandidate,"DATASET POLYDATA\n");
  fprintf(foutSpineCandidate,"POINTS %d float\n",num_nodes);

  #if InterMedial
    if(times_erosion > 1)
    {
    fprintf(tempfile1,"# vtk DataFile Version 3.0\n");
    fprintf(tempfile1,"MST of skel\n");
    fprintf(tempfile1,"ASCII\n");
    fprintf(tempfile1,"DATASET POLYDATA\n");
    fprintf(tempfile1,"POINTS %d float\n",num_nodes);
    }
  #endif
   
  //reinitialize the file and variables used to loop through it
  fin.clear();
  fin.seekg(0, ios::beg);
  num_nodes = 0;  // initial
  idx_edge = 0;   // initial
  itr = 0;
  fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
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
      fprintf(fout_Spine,"%f %f %f\n", nodePosition.x, nodePosition.y, nodePosition.z);  
      fprintf(foutSpineCandidate,"%f %f %f\n", nodePosition.x, nodePosition.y, nodePosition.z); 
      #if InterMedial
        if(times_erosion > 1)
        {
        fprintf(tempfile1,"%f %f %f\n", nodePosition.x, nodePosition.y, nodePosition.z);
        }
      #endif
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
              densityFactor =(float)  pow(double(volin[idx]+0.001), 1); // by xiao (float)
              densityFactor += (float)pow(double(volin[iidx]+0.001),1);
              densityFactor =(float) fabs(pow(double(densityFactor),double(1.05)) +pow(double(voxelNodeIndex[iidxMid1]+voxelNodeIndex[iidxMid2]+1), 0.5));   // test by square densityFactor
              edge_w[idx_edge] = (float)(sqrt(float(kk*kk + jj*jj + ii*ii)) / (densityFactor*0.02+1));
              idx_edge++;
              }
            }
          }
        }
      }
    fin >> nodePosition.x >> nodePosition.y >> nodePosition.z >> p;
    }
  
  // Save vertex position in an array
  vertexPos = new VoxelPosition[num_nodes+1];
  for (k = 0; k < sizeZ; k++)
     for (j = 0; j < sizeY; j++)
        for (i = 0; i < sizeX; i++) {
      idx = k*slsz + j*sizeX + i; 
      if (voxelNodeIndex[idx] != 0) {
        vertexPos[voxelNodeIndex[idx]].x = (float)i;
        vertexPos[voxelNodeIndex[idx]].y = (float)j;
        vertexPos[voxelNodeIndex[idx]].z = (float)k;
      }
  }

  degree_nodes = new int[num_nodes+1];
  degree_nodes_tree = new int[num_nodes+1];
  degree_nodes_buffer = new int[num_nodes+1];
  degree_nodes_initialMST = new int[num_nodes+1];
  vertBackbone = new int[num_nodes+1];
  // Initialize for all the vertices. Actual vertex index starts from 1.
  
  // -Create Graph 
  for (i=0; i<=num_nodes; i++)
    {
    degree_nodes[i] = 0;  //initialize to zeros
    degree_nodes_tree[i] = 0;  //initialize to zeros
    degree_nodes_initialMST[i] = 0;  //initialize to zeros
    vertBackbone[i] = 0;  //initialize to zeros
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

  //- MST algorithm
  kruskal_minimum_spanning_tree(g, back_inserter(spanning_tree));
 
 
  // create initial degree_nodes array
  for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    degree_nodes[source(*ei, g)] ++;
    degree_nodes[target(*ei, g)] ++;
    degree_nodes_initialMST[source(*ei, g)] ++;
    degree_nodes_initialMST[target(*ei, g)] ++;
  }

  // Create a graph for the initial MST
  Graph msTree(num_nodes+1);
  int num_edge_MST = 0;
  for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei)
  {
          add_edge(source(*ei, g), target(*ei, g), msTree);
          num_edge_MST++;
  }
  
   // - initial MST Writer
 
   Edge_iter   ei, ei_end;

#if InterMedial
   {
   if(times_erosion > 1)
   {
    line_count = 0;
    for (tie(ei, ei_end) = edges(msTree); ei != ei_end; ++ei)
     {
     line_count++; // count the number of lines output in vtk file
     }
    fprintf(tempfile1,"LINES %d %d\n", line_count, line_count*3);
    for (tie(ei, ei_end) = edges(msTree); ei != ei_end; ++ei)
     {
     fprintf(tempfile1, "2 %ld %ld\n", source(*ei, msTree) - 1,
             target(*ei, msTree) - 1);
     }
   fclose (tempfile1);
   } // end if
   }
 #endif 

  times_dilation = times_erosion;
  edge_eroded = new int[num_nodes*2];
  num_edge_eroded = 0;
  while (times_erosion !=0) {
    times_erosion--;
    for (i=1; i<=num_nodes; i++)   degree_nodes_buffer[i] = degree_nodes[i];
    for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei)
  {
      if (degree_nodes_buffer[source(*ei, g)]>0 && degree_nodes_buffer[target(*ei, g)]>0) 
    {
              if (degree_nodes_buffer[source(*ei, g)]==1 || degree_nodes_buffer[target(*ei, g)]==1) 
        {
                  degree_nodes[source(*ei, g)] --;
                  degree_nodes[target(*ei, g)] --;
          // Save the edges eroded in a stack-like array. Each edge takes two elements of the array
                  edge_eroded[num_edge_eroded*2]  = source(*ei, g); // Saving of eroded edges
                  edge_eroded[num_edge_eroded*2+1] = target(*ei, g);
                  num_edge_eroded ++;
              } //end if
      }// end if 
    }// end for
  }// end while

  // Dilation the MST by counting up the degree of nodes
  while (num_edge_eroded !=0) 
  {
    //times_dilation--;
    num_edge_eroded --;
      //for (vector < Edge >::iterator ei = spanning_tree.begin(); ei != spanning_tree.end(); ++ei) {
    edge_source = edge_eroded[num_edge_eroded*2];  // Read the stored eroded edges
    edge_target = edge_eroded[num_edge_eroded*2+1];
      if ((degree_nodes[edge_source]+ degree_nodes[edge_target]) == 1)  
    {  // if a branch tip edge
           degree_nodes[edge_source] ++;
           degree_nodes[edge_target] ++;
      } // end if
   } // end while 

  // Create a Backbone vertice flag array
  for (i=0; i<=num_nodes; i++)   {
    if (degree_nodes[i] >= 1)  
        vertBackbone[i] = 1;
  }

  //-- Label branches on Backbone of MST generated     
  Graph msTreeSpineCandidate(num_nodes+1);     // Spine Candidate graph created, is to save the possible spine
  Graph DetectedSpine(num_nodes+1); 
  
  int NumberNodesofSpineCandidate=0;
  int NumberNodesofRealSpine=0;
  //Edge_iter   ei, ei_end;
  Vertex_iter vi, vend;
  graph_traits<Graph>::out_edge_iterator  outei, outedge_end, outei2, outedge_end2, outei3, outedge_end3;

  int vertsCurBranch2[MAXNumBranch][2000];  // For the 2nd level branch at BB (at most MAXNumBranch 2nd level branches, at most 2000 vertices)
  int vertsCurBr_Index2[MAXNumBranch];      // when array at [0], it is the 1st level branch at BB

  for (i=0; i<MAXNumBranch; i++)  vertsCurBr_Index2[i]=0;  // Initialize to zeros

  // CONSIDER ALL branches on the BackBone
  num_leaves = 0;
  fprintf(fclass_identify, "ID  meanDensityBranch, length  meanVesselBranch\n");
  // PRUNING short branches on the initial MST under certain threshold (e.g. 5)

  msTree = morphGraphPrune(msTree, num_nodes, vertexPos, leaf_length);

  MDLClassifier LDA_RealSpine(3);
  MDLClassifier LDA_NonSpine(3);

  int LDA_t1= LDA_RealSpine.MeanVectorandVarianceMatrix((char *)"RealSpinePrior.txt");
  int LDA_t2= LDA_NonSpine.MeanVectorandVarianceMatrix((char *)"NonSpinePrior.txt");
  if (LDA_t1 >0)
  {
	  std::cout << " There is Real-Spine Feature Sample,We do machine learning based classification" << std::endl;
  }
  else 
  {

	  std::cout << " There is not Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
  }
  if (LDA_t2 >0)
  {
	  std::cout << " There is Non-Spine Feature Sample,We do machine learning based classification" << std::endl;
  }
  else 
  {
	  std::cout << " There is Non-Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
  }

  double sample[3];
  //ONLY run this first!
  typedef property_map<Graph, vertex_index_t>::type IndexMap;
  IndexMap index = get(vertex_index, msTree);  // get index map of vertices
  
  for(boost::tie(vi, vend) = vertices(msTree); vi != vend; ++vi)
    { // for all vertex in the graph
    vertsCurBr_Index2[0] = 0;
    index_vert = int(index[*vi]);   // Get index from the graph IndexMap
    if (vertBackbone[index_vert]==0) continue;  // if the vertex is not on Backbone, continue
    int outdegree = out_degree(*vi, msTree);
    numBranch_on_Backbone = outdegree - degree_nodes[index_vert];
      
    if (numBranch_on_Backbone <= 0)  continue;  // if it has at least one branch out of BackBone
    // For each out branch (edge) on the BackBone
    for (boost::tie(outei, outedge_end) = out_edges(*vi, msTree); outei != outedge_end; outei++)
      {
      int targ = target(*outei, msTree);
      if (vertBackbone[targ] > 0)  continue;  // the edge is on BackBone, continue
      meanDensityBranch[0] = 0;
      meanVesselBranch[0] = 0;
      vertsCurBranch2[0][vertsCurBr_Index2[0]] = source(*outei, msTree);
      vertsCurBr_Index2[0]++;
      vertsCurBranch2[0][vertsCurBr_Index2[0]] = target(*outei, msTree);
      num_leaves++;

      while (out_degree(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], msTree), msTree) == 2)
        {
        for (boost::tie(outei2, outedge_end2) =
             out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], msTree),
                       msTree);
             outei2 != outedge_end2; ++outei2)
          {
          if (target(*outei2, msTree) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
            continue;
          vertsCurBranch2[0][vertsCurBr_Index2[0]+1] = target(*outei2, msTree);
          }
        vertsCurBr_Index2[0]++;
        }
      branchChosen = 1;
      // Evaluate with MDL if the branch is chosen
      //if (i == prunetimes-1) num_leaves++;
      length_leaf[0] = 0;
      for (j = 0; j <= vertsCurBr_Index2[0]; j++)
        {
        indVert = vertsCurBranch2[0][j];
    // output the locations 
    //fprintf(fclass_identify, "%d  %6.2f %6.2f %6.2f\n", num_leaves, vertexPos[indVert].x, vertexPos[indVert].y, vertexPos[indVert].z);
        
       if (j==0)
          {
          indVert_last = indVert;
          }
        length_edge = (vertexPos[indVert].x-vertexPos[indVert_last].x)*(vertexPos[indVert].x-vertexPos[indVert_last].x);
        length_edge+= (vertexPos[indVert].y-vertexPos[indVert_last].y)*(vertexPos[indVert].y-vertexPos[indVert_last].y);
        length_edge+= (vertexPos[indVert].z-vertexPos[indVert_last].z)*(vertexPos[indVert].z-vertexPos[indVert_last].z);
        length_edge = sqrt(length_edge);
        length_leaf[0] += length_edge;
        indVert_last = indVert;

        idx =(long) (vertexPos[indVert].z *slsz + vertexPos[indVert].y *sizeX + vertexPos[indVert].x); // add long by xiao
        meanDensityBranch[0] += volin[idx];
        meanVesselBranch[0] += volvessel[idx];
        }
        meanDensityBranch[0] = meanDensityBranch[0]/ (vertsCurBr_Index2[0]+1);
        meanVesselBranch[0] = meanVesselBranch[0] / (vertsCurBr_Index2[0]+1);

        //mahalanobis_dist[0]    =      mahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 1);  
        //mahalanobis_dist_nonSpine[0] = mahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 0); 
  
       sample[0] =  meanDensityBranch[0];
       sample[1] =  length_leaf[0];
       sample[2] =  meanVesselBranch[0];
       if (LDA_t1>0)
        {
          mahalanobis_dist[0] = LDA_RealSpine.MahalanobisDist(sample); 
        }
       else 
        {  //std::cout << " There is no Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
          mahalanobis_dist[0] = LDA_RealSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 1);
        }
       if (LDA_t2>0)
        {
    
          mahalanobis_dist_nonSpine[0] = LDA_NonSpine.MahalanobisDist(sample);
         }
       else 
        {
          //std::cout << " There is no Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
          mahalanobis_dist_nonSpine[0] = LDA_NonSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 0);
        }
   
        mahalanobis_dist_min = mahalanobis_dist[0];
        mahalanobis_dist_minIndex = 0;
        // output the spine candidate feature sample;  
        //fprintf(fclass_identify, "%d  %f %f %f\n", -num_leaves, meanDensityBranch[0], length_leaf[0], meanVesselBranch[0]);
        fprintf(fclass_identify, "%d  %f %f %f\n", num_leaves, meanDensityBranch[0], length_leaf[0], meanVesselBranch[0]);
        
        if (branchChosen == 1)  
        {
         for (j = 1; j <= vertsCurBr_Index2[0]; j++) 
         {
           add_edge(vertsCurBranch2[0][j-1], vertsCurBranch2[0][j], msTreeSpineCandidate);   // add branch for the 1nd level
            NumberNodesofSpineCandidate++;
         }
        }

        // ## Begin to check the 2nd level of branches located at BB
        int ind2Brch = 0;
 
       // For each 2nd level branch starting from the end of 1st level branch
        for (boost::tie(outei2, outedge_end2) =out_edges(vertex(vertsCurBranch2[0][vertsCurBr_Index2[0]], msTree),msTree);
           outei2 != outedge_end2; ++outei2)
        {
        if (target(*outei2, msTree) == (unsigned int)vertsCurBranch2[0][vertsCurBr_Index2[0]-1])
          continue;  // continue if the out edge belongs to the old branch
        ind2Brch++;
        vertsCurBranch2[ind2Brch][0] = vertsCurBranch2[0][vertsCurBr_Index2[0]];
        vertsCurBr_Index2[ind2Brch] = 1;
        vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]] = target(*outei2, msTree);
        // Search for the end of 2nd level branch
        while (out_degree(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], msTree), msTree) == 2)
          {
          //curVert = vertex(vertsCurBranch[vertsCurBr_Index], msTree);
          for (boost::tie(outei3, outedge_end3) = out_edges(vertex(vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]], msTree), msTree);
                                                outei3 != outedge_end3; ++outei3) {
            if (target(*outei3, msTree) == (unsigned int)vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]-1])
              continue;
            vertsCurBranch2[ind2Brch][vertsCurBr_Index2[ind2Brch]+1] = target(*outei3, msTree);
            }
          vertsCurBr_Index2[ind2Brch]++;
          }
        // Compute the feature-based description length for the 2nd level branch
        length_leaf[ind2Brch] = 0;
        meanDensityBranch[ind2Brch] = 0;
        meanVesselBranch[ind2Brch] = 0;
        num_leaves++;  // for training purpose
        for (j = 0; j <= vertsCurBr_Index2[ind2Brch]; j++)
         {
          indVert = vertsCurBranch2[ind2Brch][j];
          // second level feature 
          //fprintf(fclass_identify, "%d  %6.2f %6.2f %6.2f\n", num_leaves, vertexPos[indVert].x, vertexPos[indVert].y, vertexPos[indVert].z);
         
          if (j==0)
            {
            indVert_last = indVert;
            }
          length_edge = (vertexPos[indVert].x-vertexPos[indVert_last].x)*(vertexPos[indVert].x-vertexPos[indVert_last].x);
          length_edge+= (vertexPos[indVert].y-vertexPos[indVert_last].y)*(vertexPos[indVert].y-vertexPos[indVert_last].y);
          length_edge+= (vertexPos[indVert].z-vertexPos[indVert_last].z)*(vertexPos[indVert].z-vertexPos[indVert_last].z);
          length_edge = sqrt(length_edge);
          length_leaf[ind2Brch] += length_edge;
          indVert_last = indVert;

          idx = (long)(vertexPos[indVert].z *slsz + vertexPos[indVert].y *sizeX + vertexPos[indVert].x);
           meanDensityBranch[ind2Brch] += volin[idx];
          meanVesselBranch[ind2Brch] += volvessel[idx];
          }
          meanDensityBranch[ind2Brch] = meanDensityBranch[ind2Brch]/ (vertsCurBr_Index2[ind2Brch]+1);
          meanVesselBranch[ind2Brch] = meanVesselBranch[ind2Brch] / (vertsCurBr_Index2[ind2Brch]+1);

         // Compute the average features of two level branches
         length_2leaf[ind2Brch] = length_leaf[ind2Brch] + length_leaf[0];
         aveDensityBranch[ind2Brch] = (meanDensityBranch[ind2Brch]*length_leaf[ind2Brch] +  meanDensityBranch[0]*length_leaf[0]);
         aveDensityBranch[ind2Brch] = aveDensityBranch[ind2Brch] / length_2leaf[ind2Brch];
         aveVesselBranch[ind2Brch] = (meanVesselBranch[ind2Brch]*length_leaf[ind2Brch] + meanVesselBranch[0]*length_leaf[0]);
         aveVesselBranch[ind2Brch] = aveVesselBranch[ind2Brch] / length_2leaf[ind2Brch];

         // mahalanobis distance is based on features of two-level branches
         //mahalanobis_dist[ind2Brch]    = mahalanobisDist(aveDensityBranch[ind2Brch], length_2leaf[ind2Brch], aveVesselBranch[ind2Brch], 1);
         //mahalanobis_dist_nonSpine[ind2Brch]=mahalanobisDist(aveDensityBranch[ind2Brch], length_2leaf[ind2Brch], aveVesselBranch[ind2Brch], 0);
        
	     sample[0] =  aveDensityBranch[ind2Brch];
         sample[1] =  length_2leaf[ind2Brch];
         sample[2] =  aveVesselBranch[ind2Brch];
         if (LDA_t1>0)
         {
          mahalanobis_dist[ind2Brch] = LDA_RealSpine.MahalanobisDist(sample); 
          }
         else 
         {     //std::cout << " There is no Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
          mahalanobis_dist[ind2Brch] = LDA_RealSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 1);
          }
         if (LDA_t2>0)
          {
    
          mahalanobis_dist_nonSpine[ind2Brch] = LDA_NonSpine.MahalanobisDist(sample);
          }
         else 
          {
          //std::cout << " There is no Spine Feature file for Machine Learning, thus we do the default classification" << std::endl;
          mahalanobis_dist_nonSpine[ind2Brch] = LDA_NonSpine.MahalanobisDist(meanDensityBranch[0], length_leaf[0], meanVesselBranch[0], 0);
          }

		// out put 
        //fprintf(fclass_identify, "%d  %f %f %f\n", -num_leaves, aveDensityBranch[ind2Brch], length_2leaf[ind2Brch], aveVesselBranch[ind2Brch]);
        fprintf(fclass_identify, "%d  %f %f %f\n", num_leaves, aveDensityBranch[ind2Brch], length_2leaf[ind2Brch], aveVesselBranch[ind2Brch]);      
     
        for (j = 1; j <= vertsCurBr_Index2[ind2Brch]; j++) 
        {
         // test: add any branches to the Spine Candidate
         add_edge(vertsCurBranch2[ind2Brch][j-1], vertsCurBranch2[ind2Brch][j], msTreeSpineCandidate);   // add branch for the 2nd level
         NumberNodesofSpineCandidate++;
        }

       } // End of 2nd level branch

      // ** Minimal MDL is solved by looking at each choice
      // Empty model set is chosen (no branch) 
      //----------------------------MDL fitness for spine --------------------------------------//
 
      MDL_minIndex = -1; // Indicate that empty model set is chosen
      sum_mahalanobis_nonSpine = 0;
      for (i = 0; i<= ind2Brch; i++)
        {
        sum_mahalanobis_nonSpine += mahalanobis_dist_nonSpine[i];
        }
      MDL_min = sum_mahalanobis_nonSpine;

      // 2. Only 1st level branch model set
      MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[0] + mahalanobis_dist[0];
      MDL = MDL + (1-alpha)*(1/alpha)*(-14.0/3.0);      // alpha represents the model description length of one branch
      if (MDL < MDL_min)
        {
        MDL_min = MDL;
        MDL_minIndex = 0;  
        }
      // 3. Two level branch model set, including 1st level and 2nd level branches
      for (i = 1; i <= ind2Brch; i++)
        {
        MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[i] + mahalanobis_dist[i];
        MDL = MDL + (1-alpha)*(1/alpha)*(-14.0/3.0);
        if (MDL <  MDL_min)
          {
          MDL_min = MDL;
          MDL_minIndex =  i;
          }
        } //end for
    
      // Adding Spine model - Only add the branches with the minimal MDL
       if (MDL_minIndex >= 0)
        {
        for (j = 1; j <= vertsCurBr_Index2[0]; j++)
          {
          //add_edge(vertsCurBranch2[0][j-1], vertsCurBranch2[0][j], msTreeSpineCandidate); 
      add_edge(vertsCurBranch2[0][j-1], vertsCurBranch2[0][j], DetectedSpine); 
          NumberNodesofRealSpine++;
          }
        if (MDL_minIndex >= 1)
          {
          for (j = 1; j <= vertsCurBr_Index2[MDL_minIndex]; j++)
            {
             add_edge(vertsCurBranch2[MDL_minIndex][j-1],
                     vertsCurBranch2[MDL_minIndex][j],  DetectedSpine);
       //add_edge(vertsCurBranch2[MDL_minIndex][j-1],
             //        vertsCurBranch2[MDL_minIndex][j],  msTreeSpineCandidate);
             NumberNodesofRealSpine++;
            } // end for
          }// end if
        }// end if(MDL_minIndex >= 0)

      // 1. Empty model set is chosen (no branch)
      sum_mahalanobis_nonSpine = 0;
      for (i = 0; i<= ind2Brch; i++) 
      {
        sum_mahalanobis_nonSpine += mahalanobis_dist_nonSpine[i];
      }

      // 2. Only 1st level branch model set
      MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[0] + mahalanobis_dist[0];
      MDL_minIndex = 0; // Indicate that first level branch model set is chosen
      MDL_min = MDL;

      // 3. Two level branch model set, including 1st level and 2nd level branches
      for (i = 1; i <= ind2Brch; i++) 
      {
        MDL = sum_mahalanobis_nonSpine - mahalanobis_dist_nonSpine[i] + mahalanobis_dist[i];
        if (MDL <  MDL_min) 
        {
          MDL_min = MDL;
          MDL_minIndex =  i;
        }
      }// end for 
      fprintf(fout_MDL, "%10d %20f %20f\n", ind2Brch, sum_mahalanobis_nonSpine, MDL_min);  
  
    
   }    // End of each out edge 
  } // End of all vertice

  //----------------------------End MDL fitness for spine --------------------------------------//

  //msTreeBB = msTreeBB_buffer;
 
   //fprintf(fclass_identify, "0 0 0 0 end");
   
  
   // - Spine-Candidate Writer
   line_count = 0;

   for (tie(ei, ei_end) = edges(msTreeSpineCandidate); ei != ei_end; ++ei)
     {
     line_count++; // count the number of lines output in vtk file
     }
   fprintf(foutSpineCandidate,"LINES %d %d\n", line_count, line_count*3);

   for (tie(ei, ei_end) = edges(msTreeSpineCandidate); ei != ei_end; ++ei)
     { 
       fprintf(foutSpineCandidate, "2 %ld %ld\n", source(*ei, msTreeSpineCandidate) - 1,
             target(*ei, msTreeSpineCandidate) - 1);  
     }


   // --------  prune the spine andidate to get the finall spine
   //msTreeSpineCandidate = morphGraphPrune(msTreeSpineCandidate, num_nodes, vertexPos, leaf_length);


   // - Spine Writer
   line_count = 0;
   for (tie(ei, ei_end) = edges(DetectedSpine); ei != ei_end; ++ei)
     {
     line_count++; // count the number of lines output in vtk file
     }
   fprintf(fout_Spine,"LINES %d %d\n", line_count, line_count*3);

   for (tie(ei, ei_end) = edges(DetectedSpine); ei != ei_end; ++ei)
     {
     // output lines into vtk file
     //fprintf(fout_Spine, "2 %zu %zu\n", source(*ei, msTreeSpineCandidate) - 1,
     fprintf(fout_Spine, "2 %ld %ld\n", source(*ei, DetectedSpine) - 1,
             target(*ei, DetectedSpine) - 1);
     }


  //cleanup the memory we've allocated...
  delete [] nodeIndicesInitialized;
  delete [] voxelNodeIndex;
  delete [] edge_array;
  delete [] edge_w;
  delete [] vertexPos;
  delete [] degree_nodes;
  delete [] degree_nodes_tree;
  delete [] degree_nodes_buffer;
  delete [] degree_nodes_initialMST;
  delete [] vertBackbone;
  delete [] edge_eroded;

  // release memory for malloc   by xiao liang
  free(volin);// = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  free(volvessel);// = (DATATYPEIN*)malloc(sizeX*sizeY*sizeZ*sizeof(DATATYPEIN));
  fclose(foutSpineCandidate);
  fclose(fout_Spine);
  fout_Spine = NULL;
  foutSpineCandidate = NULL;
  volin = NULL;
  volvessel = NULL;
  fclose(fout_MDL);
  cout << "Spine Extraction is done" << endl; 
  return EXIT_SUCCESS;
}


double mahalanobisDist(double meanDensityBranch, double length_leaf, double meanVesselBranch, int spineOne) {
    double x1, x2, x3;
    double mahalanobis_dist = 0;

  /* For dataset time330
  if (spineOne == 1)  {
        x1 = meanDensityBranch - 58.4550;  // minus mean_feature from matlab 
        x2 = length_leaf -    14.0834;              //if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch -  205.4605;
    mahalanobis_dist = x1*x1*0.0032+ 2*x1*x2*(0.0018)+ 2*x1*x3*(-0.001)+ x2*x2*0.0562+ 2*x2*x3*(-0.0017) +x3*x3*0.0006;
  }
  else  {
        x1 = meanDensityBranch - 13.1488;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf -     6.9848;                             
        x3 = meanVesselBranch -   59.0966;
    mahalanobis_dist = x1*x1*0.0101+ 2*x1*x2*(0.0011)+ 2*x1*x3*(-0.0025)+ x2*x2*0.0202+ 2*x2*x3*(-0.0011) +x3*x3*0.0009;
  }
  */

  // For dataset Trach6A

 // For dataset Trach6A
  /*
  if (spineOne == 1)  {
        x1 = meanDensityBranch - 44.54;  // minus mean_feature from matlab 
        x2 = length_leaf - 14.43;               //if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch - 198.62;
    mahalanobis_dist = x1*x1*0.0019+ 2*x1*x2*(0.0031)+ 2*x1*x3*(-0.0002)+ x2*x2*0.0252+ 2*x2*x3*0.0002  +x3*x3*0.0004;
  }
  else  {
        x1 = meanDensityBranch - 37.52;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf - 9.42;                             
        x3 = meanVesselBranch - 178.26;
    mahalanobis_dist = x1*x1*0.0018+ 2*x1*x2*(-0.0002)+ 2*x1*x3*(-0.0003)+ x2*x2*0.0182+ 2*x2*x3*0.0010 +x3*x3*0.0004;
  }
  
*/
  

  // For dataset MBFsp

  if (spineOne == 1)  {
        x1 = meanDensityBranch - 58.4550;  // minus mean_feature from matlab 
        x2 = length_leaf -    14.0834;              //if (x2>10)  x2=10; //keep long dendrite 
        x3 = meanVesselBranch -  58.4550;
    mahalanobis_dist = x1*x1*0.0032+ 2*x1*x2*(0.0018)+ 2*x1*x3*(-0.001)+ x2*x2*0.0562+ 2*x2*x3*(-0.0017) +x3*x3*0.0006;
  }
  else  {
        x1 = meanDensityBranch - 13.1488;  // assume non-spine has close-to-zero distribution, but the same variance as spines 
        x2 = length_leaf -     6.9848;                             
        x3 = meanVesselBranch -   13.1488;
    mahalanobis_dist = x1*x1*0.0101+ 2*x1*x2*(0.0011)+ 2*x1*x3*(-0.0025)+ x2*x2*0.0202+ 2*x2*x3*(-0.0011) +x3*x3*0.0009;
  }
  
  return mahalanobis_dist;
}

