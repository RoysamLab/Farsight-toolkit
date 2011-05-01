# include "kNearestObjects.h"


kNearestObjects::kNearestObjects(std::map< unsigned int, std::vector<float> > centroidMap)
{
	allIds = 0;
	featureTable = NULL;
	this->centerMap = centroidMap;
	sample = SampleType::New();
	treeGenerator = TreeGeneratorType::New();

	sample->SetMeasurementVectorSize( 3 );

	MeasurementVectorType mv;
	unsigned int id;
	for (It = centerMap.begin(); It != centerMap.end(); ++It )
	{
		id = (*It).first;
		mv[0] = centerMap[id].at(0);
		mv[1] = centerMap[id].at(1);
		mv[2] = centerMap[id].at(2);
		idToCentroidMap[id] = mv;
		sample->PushBack( mv );
    }

	treeGenerator->SetSample( sample );
	treeGenerator->SetBucketSize( 16 );
	treeGenerator->Update();

	tree = treeGenerator->GetOutput();
}


//***********************************************************************************************
// K Nearest Neighbors
//***********************************************************************************************

std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects::k_nearest_neighbors_All(unsigned int k, unsigned short Class_dest, unsigned short Class_src)
{	
	std::vector<std::vector< std::pair<unsigned int, double> >> kNearestIDs;
	for(IdIt = idToCentroidMap.begin(); IdIt != idToCentroidMap.end(); ++IdIt )
	{
		std::vector< std::pair<unsigned int, double> > kNearestIds;
		unsigned int ID = (*IdIt).first;
		if(Class_src == 0)
			kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
		else
		{
			for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
			{
				if( (ID == featureTable->GetValueByName(row, "ID").ToUnsignedInt())
					&& (featureTable->GetValueByName(row, "prediction_default1").ToUnsignedShort() == Class_src) )
				{
					kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
					kNearestIDs.push_back(kNearestIds);
					break;
				}
			}
		}		
	}

	return kNearestIDs;
}

std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects::k_nearest_neighbors_IDs(std::vector<unsigned int> IDs, unsigned int k, unsigned short Class_dest)
{
	std::vector<std::vector< std::pair<unsigned int, double> >> kNearestIDs;
	for(int n=0; n<(int)IDs.size(); ++n)
	{
		std::vector< std::pair<unsigned int, double> > kNearestIds;
		unsigned int ID = IDs.at(n);
		kNearestIds = k_nearest_neighbors_ID(ID, k, Class_dest);
		kNearestIDs.push_back(kNearestIds);				
	}

	return kNearestIDs;
}

std::vector< std::pair<unsigned int, double> > kNearestObjects::k_nearest_neighbors_ID(unsigned int id, unsigned int k, unsigned short Class_dest)
{
	std::vector< std::pair<unsigned int, double> > kNearestIds;
	
	if(Class_dest == 0)
	{
		DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
		DistanceMetricType::OriginType origin( 3 );
		for (unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i)
		{
			origin[i] = idToCentroidMap[id][i];
		}
		distanceMetric->SetOrigin( origin );

		TreeType::InstanceIdentifierVectorType kNeighbors;
		tree->Search( idToCentroidMap[id], k+1, kNeighbors ); 
		kNearestIds.push_back( std::make_pair(id,0) );
		for ( unsigned int i = 0 ; i < kNeighbors.size() ; ++i )
		{
		    std::map<int, MeasurementVectorType>::iterator iter;
			MeasurementVectorType mvt = tree->GetMeasurementVector(kNeighbors[i]);
			for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
			{
				if(((*iter).second == mvt) && ((*iter).first!=id))
				{
					double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( kNeighbors[i] ));
					std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
					kNearestIds.push_back(Neighbor);
					break;				
				}
			}
		}
	}

	else
	{
		kNearestIds = neighborsWithinRadius_ID(id, 1e5, Class_dest);
		unsigned int vector_size = (unsigned int)kNearestIds.size();
		double max_dist;
		unsigned int max_pos;
		while((unsigned int)kNearestIds.size() != k+1)
		{
			max_dist = 0;
			for(unsigned int x=1; x<k+1; ++x)
			{
				if(kNearestIds.at(x).second > max_dist)
				{
					max_dist = kNearestIds.at(x).second;
					max_pos = x;
				}
			}
			if(kNearestIds.at(k+1).second < max_dist )
				kNearestIds.erase(kNearestIds.begin()+max_pos);
			else
				kNearestIds.erase(kNearestIds.begin()+k+1);
		}
	}

	return kNearestIds;
}

//***********************************************************************************************
// Neighbors with certain radius
//***********************************************************************************************

std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects::neighborsWithinRadius_All(double radius, unsigned short Class_dest, unsigned short Class_src)
{	
	std::vector<std::vector< std::pair<unsigned int, double> >> inRadiusIDs;
	for(IdIt = idToCentroidMap.begin(); IdIt != idToCentroidMap.end(); ++IdIt )
	{
		std::vector< std::pair<unsigned int, double> > inRadiusIds;
		unsigned int ID = (*IdIt).first;
		if(Class_src == 0)
			inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
		else
		{
			for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
			{
				if( (ID == featureTable->GetValueByName(row, "ID").ToUnsignedInt())
					&& (featureTable->GetValueByName(row, "prediction_default1").ToUnsignedShort() == Class_src) )
				{
					inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
					inRadiusIDs.push_back(inRadiusIds);
					break;
				}
			}
		}		
	}

	return inRadiusIDs;
}

std::vector< std::vector< std::pair<unsigned int, double> > > kNearestObjects::neighborsWithinRadius_IDs(std::vector<unsigned int> IDs, double radius, unsigned short Class_dest)
{
	std::vector<std::vector< std::pair<unsigned int, double> >> inRadiusIDs;
	for(int n=0; n<(int)IDs.size(); ++n)
	{
		std::vector< std::pair<unsigned int, double> > inRadiusIds;
		unsigned int ID = IDs.at(n);
		inRadiusIds = neighborsWithinRadius_ID(ID, radius, Class_dest);
		inRadiusIDs.push_back(inRadiusIds);				
	}

	return inRadiusIDs;
}

std::vector< std::pair<unsigned int, double> > kNearestObjects::neighborsWithinRadius_ID(unsigned int id, double radius, unsigned short Class_dest)
{
	DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
	DistanceMetricType::OriginType origin( 3 );
	for ( unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i )
	{
		origin[i] = idToCentroidMap[id][i];
    }
	distanceMetric->SetOrigin( origin );

	TreeType::InstanceIdentifierVectorType radNeighbors;
	tree->Search( idToCentroidMap[id], radius, radNeighbors ); 
	std::vector< std::pair<unsigned int, double> > inRadiusIds;
	inRadiusIds.push_back( std::make_pair(id,0) );
	for ( unsigned int i = 0 ; i < radNeighbors.size() ; ++i )
	{
	    std::map<int, MeasurementVectorType>::iterator iter;
		MeasurementVectorType mvt = tree->GetMeasurementVector(radNeighbors[i]);
		for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
		{
			if(((*iter).second == mvt) && ((*iter).first!=id))
			{
				if(Class_dest == 0)
				{
					double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( radNeighbors[i] ));
					std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
					inRadiusIds.push_back(Neighbor);
					break;
				}
				else
				{
					for(int row = 0; row<(int)featureTable->GetNumberOfRows(); ++row)
					{
						if( ((*iter).first == featureTable->GetValueByName(row, "ID").ToUnsignedInt())
							&& (featureTable->GetValueByName(row, "prediction_default1").ToUnsignedShort() == Class_dest) )
						{
							double dist = distanceMetric->Evaluate( tree->GetMeasurementVector( radNeighbors[i] ));
							std::pair<unsigned int, double> Neighbor = std::make_pair((*iter).first, dist);
							inRadiusIds.push_back(Neighbor);
							break;
						}
					}
					break;
				}
			}
		}
	}

	return inRadiusIds;
}


//***********************************************************************************************
// Build Graph from Vectors
//***********************************************************************************************
vtkSmartPointer<vtkTable> kNearestObjects::vectorsToGraphTable(std::vector< std::vector< std::pair<unsigned int, double> > > NeighborIDs)
{
	graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);

	allIds = 1;
	for(unsigned int j=0 ; j<NeighborIDs.size() ; j++)
	{
		this->vectorsToGraphTable(NeighborIDs[j]);
	}
	allIds = 0;

	return graphtable;
}

vtkSmartPointer<vtkTable> kNearestObjects::vectorsToGraphTable(std::vector< std::pair<unsigned int, double> > NeighborIds)
{
	if(allIds == 0)
	{
		graphtable = vtkSmartPointer<vtkTable>::New();
		graphtable->Initialize();
		vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
		column1->SetName( "Source" );
		graphtable->AddColumn(column1);
		vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
		column2->SetName( "Target" );
		graphtable->AddColumn(column2);
	}
	
	for(unsigned int i=1 ; i<(int)NeighborIds.size() ; i++)
	{
		vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	model_data1->InsertNextValue(static_cast<vtkVariant>(NeighborIds.at(0).first));
		model_data1->InsertNextValue(static_cast<vtkVariant>(NeighborIds.at(i).first));
		graphtable->InsertNextRow(model_data1);	
	}
	return graphtable;
}

//***********************************************************************************************
// Build Graph from Vectors
//***********************************************************************************************

kNearestObjects::NeighborGraph kNearestObjects::getNeighborGraph(std::vector< std::vector< std::pair<unsigned int, double> > > NeighborIDs)
{
	nodeName  = get(vertex_name, NG);	
	allIds = 1;
	for(unsigned int j=0 ; j<NeighborIDs.size() ; j++)
	{
		this->getNeighborGraph(NeighborIDs[j]);
	}
	allIds = 0;

	return NG;
}

kNearestObjects::NeighborGraph kNearestObjects::getNeighborGraph(std::vector< std::pair<unsigned int, double> > NeighborIds)
{
	if(allIds == 0)
	{
		NG.clear();
		nodeName  = get(vertex_name, NG);
	}

	int src = GetNodeIndex(NeighborIds.at(0).first);
	if(src == -1)
	{
		node S;
		S = add_vertex(NG);
		nodeName[S] = convert2string(NeighborIds.at(0).first);
		src = num_vertices(NG)-1;
	}

	for(unsigned int i=1 ; i<(int)NeighborIds.size() ; i++)
	{
		int trg = GetNodeIndex(NeighborIds[i].first);
		if(trg ==-1)
		{
			node T;
			T = add_vertex(NG);
			nodeName[T] = convert2string(NeighborIds[i].first);
			trg = num_vertices(NG)-1;
		}
		bool bRet;
		Edge e;
		tie(e,bRet) = edge(src,trg,NG);
		if(!bRet)
		{
			add_edge(src,trg,NG);
		}
	}

	return NG;
}

std::string kNearestObjects::convert2string(unsigned int id)
{
	std::stringstream out;
	out << id;
	std::string s = out.str();
	return s;
}

int kNearestObjects::GetNodeIndex(unsigned int id)
{	
	std::string id1 = convert2string(id);
	int index;
	int flag = -1;
	for ( unsigned int i = 0 ; i < num_vertices(NG) ; ++i)
	{
		if(nodeName[i]==id1)
		{
			flag = 0;
			index = i;
			break;
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

//***********************************************************************************************
// Build Table from graph
//***********************************************************************************************

vtkSmartPointer<vtkTable> kNearestObjects::graphToTable(NeighborGraph g)
{
	graphtable = vtkSmartPointer<vtkTable>::New();
	graphtable->Initialize();
	vtkSmartPointer<vtkStringArray> column1 = vtkSmartPointer<vtkStringArray>::New();
    column1->SetName( "Source" );
    graphtable->AddColumn(column1);
    vtkSmartPointer<vtkStringArray> column2 = vtkSmartPointer<vtkStringArray>::New();
    column2->SetName( "Target" );
	graphtable->AddColumn(column2);

	for (tie(vi,vi_end) = vertices(g) ; vi != vi_end ; ++vi)
    {
	    for (tie(ai,ai_end) = adjacent_vertices(*vi, g) ; ai != ai_end ; ++ai)
	    {
         //   check = 1;
	        //for(int i=0; i<(int)graphtable->GetNumberOfRows(); ++i)
	        //{
		       // if(((static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,1))) || 
		 	     //  ((static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,1))))
		       // {
		 	     //   check = 0;
		 	     //   break;
		       // }
 	       // }
	        //if(check == 1)
	        //{
	            vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*vi]));
		        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*ai]));
		        graphtable->InsertNextRow(model_data1);		  
	        //}	   
        }
    }
    return graphtable;
}
	





