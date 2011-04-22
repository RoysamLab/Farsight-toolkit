# include "kNearestObjects.h"


kNearestObjects::kNearestObjects(std::map< unsigned int, std::vector<float> > centroidMap)
{
	allIds = 0;
	this->centerMap = centroidMap;
	sample = SampleType::New();
	treeGenerator = TreeGeneratorType::New();
}

std::vector<std::vector<unsigned int>> kNearestObjects::k_nearest_neighbors_All(unsigned int k)
{
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
	
	std::vector<std::vector<unsigned int>> kNearestIDs;
	for(IdIt = idToCentroidMap.begin(); IdIt != idToCentroidMap.end(); ++IdIt )
	{
		std::vector<unsigned int> kNearestIds;
		kNearestIds = k_nearest_neighbors_ID((*IdIt).first, k);
		kNearestIDs.push_back(kNearestIds);
	}

	return kNearestIDs;
}

std::vector<unsigned int> kNearestObjects::k_nearest_neighbors_ID(unsigned int id, unsigned int k)
{
	DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();
	DistanceMetricType::OriginType origin( 3 );
	for ( unsigned int i = 0 ; i < (unsigned int)sample->GetMeasurementVectorSize() ; ++i )
	{
		origin[i] = idToCentroidMap[id][i];
    }
	distanceMetric->SetOrigin( origin );

	TreeType::InstanceIdentifierVectorType neighbors;
	tree->Search( idToCentroidMap[id], k+1, neighbors ); 
	std::vector<unsigned int> kNearestIds;
	kNearestIds.push_back(id);
	for ( unsigned int i = 0 ; i < k+1 ; ++i )
	{
	    std::map<int, MeasurementVectorType>::iterator iter;
		MeasurementVectorType mvt = tree->GetMeasurementVector(neighbors[i]);
		for(iter = idToCentroidMap.begin(); iter != idToCentroidMap.end(); ++iter)
		{
			if((*iter).second == mvt)
			{
				if((*iter).first!=id)
				{
					kNearestIds.push_back((*iter).first);
					break;
				}
			}
		}
	}

	return kNearestIds;
}

kNearestObjects::kNeighborGraph kNearestObjects::kNearestGraph(std::vector<unsigned int> kNearestIds)
{
	if(allIds == 0)
	{
		KNG.clear();
		nodeName  = get(vertex_name, KNG);
	}

	int src = GetNodeIndex(kNearestIds.at(0));
	if(src == -1)
	{
		node S;
		S = add_vertex(KNG);
		nodeName[S] = convert2string(kNearestIds.at(0));
		src = num_vertices(KNG)-1;
	}

	for(unsigned int i=1 ; i<kNearestIds.size() ; i++)
	{
		int trg = GetNodeIndex(kNearestIds[i]);
		if(trg ==-1)
		{
			node T;
			T = add_vertex(KNG);
			nodeName[T] = convert2string(kNearestIds[i]);
			trg = num_vertices(KNG)-1;
		}
		bool bRet;
		Edge e;
		tie(e,bRet) = edge(src,trg,KNG);
		if(!bRet)
		{
			add_edge(src,trg,KNG);
		}
	}

	return KNG;
}

kNearestObjects::kNeighborGraph kNearestObjects::kNearestGraph(std::vector<std::vector<unsigned int>> kNearestIDs)
{
	nodeName  = get(vertex_name, KNG);	
	allIds = 1;
	for(unsigned int j=0 ; j<kNearestIDs.size() ; j++)
	{
		this->kNearestGraph(kNearestIDs[j]);
	}
	allIds = 0;

	return KNG;
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
	for ( unsigned int i = 0 ; i < num_vertices(KNG) ; ++i)
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

vtkSmartPointer<vtkTable> kNearestObjects::kNeighborTable(kNeighborGraph g)
{
	vtkSmartPointer<vtkTable> graphtable = vtkSmartPointer<vtkTable>::New();
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
            check = 1;
	        for(int i=0; i<(int)graphtable->GetNumberOfRows(); ++i)
	        {
		        //std::cout<<i<<std::endl;
		        if(((static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,1))) || 
		 	       ((static_cast<vtkVariant>(nodeName[*ai])==graphtable->GetValue(i,0))&&(static_cast<vtkVariant>(nodeName[*vi])==graphtable->GetValue(i,1))))
		        {
		 	        check = 0;
		 	        break;
		        }
 	        }
	        if(check == 1)
	        {
	            vtkSmartPointer<vtkVariantArray> model_data1 = vtkSmartPointer<vtkVariantArray>::New();
     	        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*vi]));
		        model_data1->InsertNextValue(static_cast<vtkVariant>(nodeName[*ai]));
		        graphtable->InsertNextRow(model_data1);		  
	        }	   
        }
    }
    return graphtable;
}
	





