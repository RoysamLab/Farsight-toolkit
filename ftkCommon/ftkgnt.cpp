#include "ftkgnt.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////
// BUILD THE REGION ADJACENCY graph1 IN AN ADJACENCY_LIST 
// LOGIC : ALL NODES TRAVERSED AND IF NO CHANGE IN NUMBER OF VERTICES AND EDGES, RETURN THE graph1
///////////////////////////////////////////////////////////////////////////////////////////////////////


ftkgnt::ftkgnt()
{
	hypotheses.clear();
	MAX_VOL = 1e5;
	MAX_DEPTH = 6;
 }


void ftkgnt::runLabFilter(InputImageType::Pointer input, OutputImageType::Pointer output, bool CytoImage)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// NEED TO CALCULATE THE NEIGHBORHOOD INFORMATION FOR ALL IDs
	///////////////////////////////////////////////////////////////////////////////////////////////////////	
	cyto_image = CytoImage;
	labelFilter = FeatureCalcType::New();
	labelFilter->SetCompleteImageInputs( input, output, cyto_image );
	labelFilter->SetLevel(3);
	labelFilter->ComputeHistogramOn();
	labelFilter->Update();

}

void ftkgnt::setmaxVol(double vol)
{
	this->MAX_VOL = vol;

}

void ftkgnt::setmaxDepth(int depth)
{
	this->MAX_DEPTH = depth;
}


ftkgnt::RAGraph ftkgnt::BuildRAG(unsigned short id)
{
	Initialize(id);
	unsigned int counter = 0; 
	boost::property_map<RAGraph, boost::vertex_name_t>::type node_name = get(boost::vertex_name, this->RAG);
	
	while (counter != num_vertices(this->RAG))
	{
		id = static_cast<unsigned long>(atoi(node_name[counter].c_str()));
		std::vector<unsigned short> RAG_cells = labelFilter->GetContactNeighbors(id);
		RAG_cells.erase (RAG_cells.begin());
		
		//Get the Source Vertex for the iteration
		int root = GetNodeIndex(id,this->RAG);
		if(RAG_cells.size()>0)
		{
			for(unsigned int i =0 ; i<RAG_cells.size() ; i++)
			{	
				int tail = GetNodeIndex(RAG_cells[i],this->RAG);
				if(tail ==-1)
				{
					node V;
					V = add_vertex(this->RAG);
					node_name[V] = convert2string(RAG_cells[i]);
					tail = num_vertices(this->RAG)-1;
				}
				bool bRet;
				Edge e;
				boost::tie(e,bRet) = edge(root,tail,this->RAG);
				if(!bRet)
				{
					add_edge(root,tail,this->RAG);
				}
			}
		}
		counter = counter + 1 ;
	} 
	return this->RAG;
}



void ftkgnt::Initialize(unsigned short id)
{		
		//Clear all the edges and vertices of the graphs if present. 
		this->RAG.clear();	
		// Property accessors
		node_name nodeName  = get(boost::vertex_name,this->RAG);
		std::stringstream out;
		out << id;
		std::string s = out.str();

		//Vertex V
		node V;
		V = add_vertex(this->RAG);
		nodeName[V] = s;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
// CONVERT THE INPUT DATA TYPE TO STRING 
// NEED TO TEMPLATE
///////////////////////////////////////////////////////////////////////////////////////////////////////
std::string ftkgnt::convert2string(unsigned short id)
{
	std::stringstream out;
	out << id;
	std::string s = out.str();
	return s;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
// CHECK IF THE NODE ALREADY EXISTS IN THE graph1
// WORKS ON THE ASSUMPTION THAT THERE ARE NO DUPLICATE NODES IN THE graph1 (TRUE FOR VALID RAGs)
///////////////////////////////////////////////////////////////////////////////////////////////////////

int ftkgnt::GetNodeIndex(unsigned short id,ftkgnt::RAGraph graph1)
{
	
	std::string id1 = convert2string(id);

	int index;
	int flag = -1;
	boost::property_map<ftkgnt::RAGraph, boost::vertex_name_t>::type nodes = get(boost::vertex_name, graph1);
	for ( unsigned int i = 0 ; i < num_vertices(graph1) ; ++i)
	{
		if(nodes[i]==id1)
		{
			flag = 0;
			index = i;
			//NEED TO CHECK IF A BREAK CAN BE ADDED ?
		}
	}
	
	if (flag==0) 
	{
		return index;
	}
	else
	{
		return flag;
	}
}


void ftkgnt::setFeats(std::vector<FeaturesType> allFeat,std::vector< unsigned short > labelIndex)
{
	this->allFeat = allFeat;
	this->labelIndex = labelIndex;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// BUILD THE MERGE TREE FOR THE REGION ADJACENY graph1 WITH "id" AS THE ROOT
///////////////////////////////////////////////////////////////////////////////////////////////////////

ftkgnt::MTreeType ftkgnt::BuildMergeTreeDcon(ftkgnt::RAGraph R1, unsigned short id,std::vector< std::set<int> > hypothesis)
{
	//Create the Merge Tree
	ftkgnt::MTreeType mTree; 
	std::set<int> currRPS;
	std::set<int> nextRPS;	
	
	//To store the RPSs @ the current depth
	std::vector<std::set<int> > curr_depth_RPS;
	
	//To store the RPSs @ all depths
	std::vector< std::vector< std::set<int> > > depth_map;  
	
	//This vector stores the current vector of labels (in int ) 
	std::set< int > curr_members;
	
	unsigned int depth;
	
	//Add the root node to the Merge Tree
	std::string s = convert2string(id);
	ftkgnt::node_mt V;
	V = add_vertex(mTree);
	
	//Current Root Path Set 	
	currRPS.insert(id);
	
	mTree[V].label =  s;
	mTree[V].RPS = currRPS;
	curr_depth_RPS.push_back(currRPS);
	depth_map.push_back(curr_depth_RPS);
	depth = 0;
	curr_depth_RPS = depth_map[depth];
	
	//Add the root as a member and get the current volume
	curr_members.insert(static_cast<int>(id));
	
	// Adjacency Iterators will iterate through the adjacent vertex of the 
	// Region Adjacency graph1 a.k.a the Adjacency_list
	ftkgnt::AdjVertIt avi, avinext, av_end;
	std::set<int>::iterator it;
	std::set<int>::iterator RPSIterator;
	std::set<int>::iterator volIterator;
	std::set<int>::iterator nRPSIterator;
	std::set<int>::iterator RPSIterator2;
	boost::property_map<ftkgnt::RAGraph, boost::vertex_name_t>::type nodes = get(boost::vertex_name, R1);
	
	//Start the loop here
	
	// Logic: For each node in the tree go through the Root Path Set 
	// For every element in the root path set get the neighbors  
	// Add the neighbor to the tree if valid.
	// Stop when All nodes traversed and if no change in number of vertices and 
	// number of edges,return the graph1
	
	unsigned int counter = 0; 
	while (counter != num_vertices(mTree))
	{
		currRPS =  mTree[counter].RPS;
		
		ftkgnt::node_mt V2 = vertex(counter,mTree);
		for(RPSIterator = currRPS.begin(); RPSIterator != currRPS.end(); RPSIterator++)
		{	
			int vertex_index = GetNodeIndex(static_cast<unsigned short>(*RPSIterator),R1);
			ftkgnt::node v = vertex(vertex_index,R1);
			boost::tie(avi, av_end)=adjacent_vertices(v,R1); 
			
			for (avi=avi; avi < av_end ; ++avi)
			{
				nextRPS = currRPS;	   
				ftkgnt::node X = *avi;
				int neighbor = atoi(nodes[X].c_str());
				
				//if "it" points to currRPS.end(), then this node is not present in 
				// the current RPS. RPS condition in Gang's paper   
				nextRPS.insert(neighbor);	
				
				it=currRPS.find(neighbor);
				depth = nextRPS.size() - 1 ;		   	    
				bool depth_cond = true;
				
				// Check if nextRPS is present in the depthmap for the current depth
				//This is the depth condition in Gang's paper.
				if(depth <= depth_map.size()-1)    
				{
					curr_depth_RPS= depth_map[depth]; 
					depth_cond = (curr_depth_RPS.end() == find(curr_depth_RPS.begin(), curr_depth_RPS.end(), nextRPS));	
				}  
				
				if(it==currRPS.end() && depth_cond)
				{
					//This condition checks if the current node is not @  
					// a new level/depth in the tree in the tree
					if(depth <= depth_map.size()-1) 
					{
							curr_depth_RPS= depth_map[depth];
							curr_depth_RPS.push_back(nextRPS);   
							depth_map[depth] =  curr_depth_RPS;     
					}	
					
					// If it is at the new depth.. first check the minimum volume @ the max depth
					// If this value is > than the limit of the cells... return the tree	   
					// If not update the 	   
					else
					{
						bool dcon = (depth<MAX_DEPTH);
						if(dcon)
						{
								curr_depth_RPS.clear();
								curr_depth_RPS.push_back(nextRPS);   
								depth_map.push_back(curr_depth_RPS);
							}
													
						else
						{
							return mTree;
						}
						
					}
					
					//Check if this hypothesis has been checked previously i.e. if this combination of nodes occured 
					// in a previous merge tree 
					bool hypo_cond;
					hypo_cond = (hypothesis.end() == find(hypothesis.begin(), hypothesis.end(), nextRPS));
						
					double vol = 0;	
					for(volIterator = nextRPS.begin(); volIterator != nextRPS.end(); volIterator++)
					{	
						std::vector<unsigned short>::iterator posn1 = find(labelIndex.begin(), labelIndex.end(), *volIterator);
						ftk::IntrinsicFeatures  features  =	allFeat[posn1 - labelIndex.begin()];
						vol+= features.ScalarFeatures[ftk::IntrinsicFeatures::VOLUME];
					}

					
					//Check for the volume condition 
					//Prevents unnecessary extension of tree branches in clusters   
					bool vol_cond = (vol<MAX_VOL);

					if(hypo_cond && vol_cond)
					{
						ftkgnt::node_mt V;
						V = add_vertex(mTree);
						mTree[V].label = nodes[X];
						mTree[V].RPS = nextRPS;
						int tail = num_vertices(mTree)-1;
						add_edge(counter,tail,mTree);
					}
					
				}
				
				//Delete nextRPS
				nextRPS.clear();
				
			}
			
		}		
		counter = counter +1;
		
	}
	
	return mTree;
	
}





