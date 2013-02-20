#include "ftkGraphs/kNearestObjects.h"
#include "ftkCommon/ftkUtils.h"

double average(std::vector< std::pair<unsigned int, double> >);
int main(int argc, char* argv[])
{	

	if(argc < 3)
	{
		std::cout << "SpatialMeasurements.exe <InputTableName> <NumberOfCellTypes>" << std::endl;
		return -1;
	}

	std::string MyName = argv[0];					
	std::string InputTableName = argv[1];					
	int NumberOfCellTypes = atoi(argv[2]);

	//std::string InputTableName = "G:/DATA/0131/microglia_types.txt";					
	//int NumberOfCellTypes = 4;



	vtkSmartPointer<vtkTable> table = ftk::LoadTable(InputTableName);

	for(int row=(int)table->GetNumberOfRows()-1; row>0; --row)
	{
		int electrode_number = table->GetValueByName(row,"prediction_active_mg_neu").ToInt() - NumberOfCellTypes;
		if(electrode_number <= 0) break;
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		std::stringstream ss;
		ss << electrode_number;
		std::string column_name = "Dist_Electrode_" + ss.str();
		column->SetName(column_name.c_str());
		column->SetNumberOfValues(table->GetNumberOfRows());
		table->AddColumn(column);

		double x1 = 0.267 * table->GetValue(row,1).ToDouble();
		double y1 = 0.267 * table->GetValue(row,2).ToDouble();
		double z1 = 0.3 * table->GetValue(row,3).ToDouble();
		for(int r=0; r<(int)table->GetNumberOfRows(); ++r)
		{
			double x2 = 0.267 * table->GetValue(r,1).ToDouble();
			double y2 = 0.267 * table->GetValue(r,2).ToDouble();
			double z2 = 0.3 * table->GetValue(r,3).ToDouble();
			double dist = std::sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
			table->SetValue(r, table->GetNumberOfColumns()-1, dist);
			//std::cout << "ID: " << table->GetValue(r,0).ToInt() << "\r";
		}
	}

	std::vector< std::map< unsigned int, std::vector<double> > > centroidMaps;
	centroidMaps.resize(NumberOfCellTypes);
	std::vector< std::vector<unsigned int> > ID_vectors;
	ID_vectors.resize(NumberOfCellTypes);
	std::vector< kNearestObjects<3>* > kd_trees;
	kd_trees.resize(NumberOfCellTypes);
	std::vector< std::vector<std::vector< std::pair<unsigned int, double> > > > kNeighborID_vectors;
	kNeighborID_vectors.resize(NumberOfCellTypes);

	for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
	{
		int cell_class = table->GetValueByName(row,"prediction_active_mg_neu").ToInt();
		if(cell_class > NumberOfCellTypes) break;
		unsigned int id = table->GetValue(row,0).ToUnsignedInt();
		ID_vectors[cell_class-1].push_back(id);
		std::vector<double> c;
		c.push_back(0.267 * table->GetValue(row,1).ToDouble());
		c.push_back(0.267 * table->GetValue(row,2).ToDouble());
		c.push_back(0.3 * table->GetValue(row,3).ToDouble());		
		centroidMaps[cell_class-1][id] = c;
		//std::cout << "ID: " << table->GetValue(row,0).ToInt() << "\r";
	}

	for(int i=0; i<NumberOfCellTypes; ++i)
	{
		kd_trees[i] = new kNearestObjects<3>(centroidMaps[i]);
		kNeighborID_vectors[i] = kd_trees[i]->k_nearest_neighbors_IDs(ID_vectors[i], 5, 0);
		std::stringstream ss2;
		ss2 << i+1;
		std::string full_string = "Class_" + ss2.str() + "_proximity";
		vtkSmartPointer<vtkDoubleArray> column = vtkSmartPointer<vtkDoubleArray>::New();
		column->SetName(full_string.c_str());
		column->SetNumberOfValues((int)table->GetNumberOfRows());
		table->AddColumn(column);
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			table->SetValueByName(row, full_string.c_str(), 0);
		}
		for(int j=0; j < (int)kNeighborID_vectors[i].size(); ++j)
		{
			int Id = kNeighborID_vectors[i].at(j).at(0).first;
			double avg_dist = average(kNeighborID_vectors[i].at(j));
			for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
			{
				if(table->GetValue(row,0).ToInt() == Id)
				{
					table->SetValueByName(row, full_string.c_str(), vtkVariant(avg_dist));
					break;
				}
			}
		}
	}

	std::string::iterator it;
	it = InputTableName.end() - 4;
	InputTableName.erase(it, it+4);
	std::string OutputTableName = InputTableName + "_spat_stats.txt";
	ftk::SaveTable(OutputTableName, table);

	return 0;

}


double average(std::vector< std::pair<unsigned int, double> > ID)
{
	double dist = 0;
	for(int k=1; k<(int)ID.size(); ++k)
	{
		dist += ID.at(k).second;
	}
	double average = dist/(int)(ID.size()-1);
	return average;
}
