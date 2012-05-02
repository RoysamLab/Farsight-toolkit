#include "Heatmap.h"

BiHeatmap::BiHeatmap(QWidget *parent)
: QMainWindow(parent)
{
	this->treedata1 = false;
	this->treedata2 = false;
	this->dragLineFlag1 = false;
	this->dragLineFlag2 = false;
	this->removeActorflag = false;	
	this->localselection = false;
	this->num_selection_area = 0;
}

BiHeatmap::~BiHeatmap()
{
}

void BiHeatmap::setDataForHeatmap(std::vector<int > & order1,std::vector<int > & order2)
{
	this->rowMapFromOriginalToReorder.clear();
	this->columnMapFromOriginalToReorder.clear();

	this->order1 = order1;
	this->order2 = order2;

	for(int i=0; i<num_rows; i++)
		this->rowMapFromOriginalToReorder.insert( std::pair< int, int>(order1[i], i));

	for(int i=0; i<num_cols; i++)
		this->columnMapFromOriginalToReorder.insert( std::pair< int, int>(order2[i], i));
}

void BiHeatmap::creatDataForHeatmap()
{
	this->data.resize(this->table->GetNumberOfRows());
	for(int i = 0; i < this->table->GetNumberOfRows(); i++)
		for(int j = 1; j < this->table->GetNumberOfColumns(); j++)
			this->data[i].push_back(this->table->GetValue(i,j).ToDouble());
	this->normalize();

	std::vector<double > temp;
	temp.resize(num_cols);

	std::vector<std::vector<double > > tempdata;

	for(int i = 0; i < this->num_rows; i++)
	{
		for(int j = 0; j < this->num_cols; j++)
		{
			temp[j] = this->data[i][order2[j]];
		}
	tempdata.push_back(temp);
	}

	for(int i = 0; i < this->num_rows; i++)
		this->data[this->num_rows - i - 1] = tempdata[order1[i]]; 
	if(this->treedata1 == true)
		this->creatDataForTree1();
	if(this->treedata2 == true)
		this->creatDataForTree2();
}
void BiHeatmap::normalize()
{
	for(int i = 0; i < this->num_rows; i++)
	{
		for(int j = 0; j < this->num_cols; j++)
		{
			if(this->data[i][j] > 1)
				this->data[i][j] = 1 + log(this->data[i][j]);
			else if(this->data[i][j] < -1)
				this->data[i][j] = -1 - log(-(this->data[i][j]));
			else
				this->data[i][j] = this->data[i][j];	
		}
	}

	//int numr = this->data.size();
	//int numc = this->data[0].size();

	//for(int j = 0; j < numc; j++)
	//{
	//	double mean = 0.0;
	//	for(int i = 0; i < numr; i++)
	//	{
	//		mean += this->data[i][j];
	//	}
	//	mean /= numr;

	//	double sum = 0.0;
	//	for(int i = 0; i < numr; i++)
	//	{
	//		sum += ( this->data[i][j] - mean ) * ( this->data[i][j] - mean );
	//	}
	//	double std = sqrt(sum / numr);

	//	for(int i = 0; i < numr; i++)
	//	{
	//		this->data[i][j]= ((this->data[i][j] - mean) / std);
	//	}
	//}
}

void BiHeatmap::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels, ObjectSelection * sels2)
{
	this->table = table;
	for( int i = 0; i < this->table->GetNumberOfRows(); i++)
	{
		for( int j = 1; j < this->table->GetNumberOfColumns(); j++)
		{
			 double var = this->table->GetValue(i, j).ToDouble();
			 if( !boost::math::isnan(var))
			 {}
			 else
			 {
				 if(i>0)
					 vtkVariant value = this->table->GetValue(i-1, j);
				 else 
					 vtkVariant value = this->table->GetValue(i+1, j);
				 this->table->SetValue(i, j, 0);
				 break;
			 }
		}
	}

	this->num_rows = this->table->GetNumberOfRows();
	this->num_cols = this->table->GetNumberOfColumns() - 1;

	this->indMapFromVertexToInd.clear();
	this->indMapFromIndToVertex.clear();
	if( this->table)
	{
		for( int i = 0; i < this->table->GetNumberOfRows(); i++)
		{
			int var = this->table->GetValue( i, 0).ToInt();
			this->indMapFromVertexToInd.insert( std::pair< int, int>(var, i));
			this->indMapFromIndToVertex.push_back( var);
		}
	}

	if(!sels)
		this->Selection = new ObjectSelection();
	else
		this->Selection = sels;

	if(!sels2)
		this->Selection2 = new ObjectSelection();
	else
		this->Selection2 = sels2;

	connect(Selection, SIGNAL(changed()), this, SLOT(GetSelecectedIDs()));
}

void BiHeatmap::showHeatmap()
{	
	this->view = vtkSmartPointer<vtkGraphLayoutView>::New();
	this->creatDataForHeatmap();
	if(this->treedata1 == true && this->treedata2 == true)
		this->drawPoints();
	this->SetInteractStyle();

	this->plane = vtkSmartPointer<vtkPlaneSource>::New();
    this->plane->SetXResolution(this->num_cols);
    this->plane->SetYResolution(this->num_rows);

	this->cellData = vtkSmartPointer<vtkFloatArray>::New();
	int index = 0;
	for (int i = 0; i < this->num_rows; i++)
		for(int j = 0; j < this->num_cols; j++)
			cellData->InsertNextValue(index++);

	this->celllut = vtkSmartPointer<vtkLookupTable>::New();
	this->celllut->SetNumberOfTableValues(this->num_rows*this->num_cols);
	this->celllut->SetTableRange(0, this->num_rows*this->num_cols - 1);   
	this->celllut->Build();
	int k = 0;
	for(int i = 0; i < this->num_rows; i++)
	{
		for(int j = 0; j < this->num_cols; j++)
		{
			rgb rgb = GetRGBValue( this->data[num_rows - i - 1][j]);
			celllut->SetTableValue(k++, rgb.r, rgb.g, rgb.b);
		}
	}

	this->WriteFile("data_for_maping.txt");
	this->plane->Update();
	this->plane->GetOutput()->GetCellData()->SetScalars(cellData);

	this->cellmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->cellmapper->SetInputConnection(plane->GetOutputPort());
	this->cellmapper->SetScalarRange(0, this->num_rows*this->num_cols - 1);
	this->cellmapper->SetLookupTable(celllut);

	this->cellactor = vtkSmartPointer<vtkActor>::New();
	this->cellactor->SetMapper(cellmapper);

	//show scalar bar
	vtkSmartPointer<vtkLookupTable> scalarbarLut = vtkSmartPointer<vtkLookupTable>::New();
	scalarbarLut->SetTableRange (-6, 13);
	scalarbarLut->SetNumberOfTableValues(COLOR_MAP_SIZE);
	for(int index = 0; index<COLOR_MAP_SIZE;index++)
	{
		rgb rgbscalar = COLORMAP[index];
		scalarbarLut->SetTableValue(index, rgbscalar.r, rgbscalar.g, rgbscalar.b);
	}
	scalarbarLut->Build();

	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(scalarbarLut);
	scalarBar->SetTitle("Color Map");
	scalarBar->SetNumberOfLabels(10);
	scalarBar->GetTitleTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->GetLabelTextProperty()->SetColor(0,0,0);
	scalarBar->GetTitleTextProperty()->SetFontSize (10);
	scalarBar->SetMaximumHeightInPixels(1000);
	scalarBar->SetMaximumWidthInPixels(100);

	this->view->GetRenderer()->AddActor(cellactor);
	this->view->GetRenderer()->AddActor2D(scalarBar);
	this->showFeatureNames();
	this->view->Render();
}

void BiHeatmap::SetInteractStyle()
{
	this->myCellPicker = vtkSmartPointer<vtkCellPicker>::New();
	this->view->GetInteractor()->SetPicker(this->myCellPicker);
	this->myCellPicker->SetTolerance(0.004);

	this->selectionCallbackrightdown =vtkSmartPointer<vtkCallbackCommand>::New();
	this->selectionCallbackrightdown->SetClientData(this);
	this->selectionCallbackrightdown->SetCallback(SelectionCallbackFunctionRightButtonDown);

	this->selectionCallbackrightup =vtkSmartPointer<vtkCallbackCommand>::New();
	this->selectionCallbackrightup->SetClientData(this);
	this->selectionCallbackrightup->SetCallback(SelectionCallbackFunctionRightButtonUp);

	this->keyPress = vtkSmartPointer<vtkCallbackCommand>::New();
	this->keyPress->SetCallback(HandleKeyPress);
	this->keyPress->SetClientData(this);

	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::RightButtonReleaseEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyPressEvent);
	this->view->GetInteractor()->RemoveObservers(vtkCommand::KeyReleaseEvent);
	this->view->GetInteractor()->AddObserver(vtkCommand::RightButtonPressEvent, selectionCallbackrightdown);
	this->view->GetInteractor()->AddObserver(vtkCommand::RightButtonReleaseEvent, selectionCallbackrightup);
	this->view->GetInteractor()->AddObserver(vtkCommand::KeyPressEvent, this->keyPress);

	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();

	this->view->GetRenderer()->GradientBackgroundOff();
	this->view->GetRenderer()->SetBackground(1,1,1);
	this->view->GetInteractor()->Start();
}

rgb BiHeatmap::GetRGBValue(double val)
{
	int index = 6 * (val+6.5) - 1;   
	if( index >= COLOR_MAP_SIZE)
	{
		index = COLOR_MAP_SIZE - 1;
	}
	else if( index < 0)
	{
		index = 0;
	}
	return COLORMAP[index];
}

void BiHeatmap::showTree1()
{
	int num_nodes = this->tree1.coordinates.size();

	this->treecolors1 = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->treecolors1->SetNumberOfComponents(3);
	this->treecolors1->SetName("treecolors1");
	this->treepoints1 = vtkSmartPointer<vtkPoints>::New();	
	this->treelines1 = vtkSmartPointer<vtkCellArray>::New();
	unsigned char color[3] = {0, 0, 0};
//
	for(int i = 0; i < num_nodes - 1; i++)
		this->treecolors1->InsertNextTupleValue(color);
	int counter = this->num_rows;
	int num_pp_id = 0;
	for(int i = 0; i < this->tree1.treeid.size(); i++)
	{
		double pp[3];
		pp[0] = this->tree1.coordinates[counter][0];
		pp[1] = this->tree1.coordinates[counter][1];
		pp[2] = this->tree1.coordinates[counter++][2]; 
		this->treepoints1->InsertNextPoint(pp);
		for(int j = 1; j < this->tree1.treeid[i].size(); j++)
		{
			double ps[3];			
			ps[0] = this->tree1.coordinates[this->tree1.treeid[i][j]][0];
			ps[1] = this->tree1.coordinates[this->tree1.treeid[i][j]][1];
			ps[2] = this->tree1.coordinates[this->tree1.treeid[i][j]][2];
			this->treepoints1->InsertNextPoint(ps);

			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			line->GetPointIds()->SetId(0, num_pp_id);
			line->GetPointIds()->SetId(1, num_pp_id + j);
			this->treelines1->InsertNextCell(line);
		}
		num_pp_id += this->tree1.treeid[i].size();
	}
//

/////////////////////////////////////////////////////////////////////
	//int num_line = num_nodes + this->tree1.treeid.size() - 1;
	//for(int i = 0; i < num_line; i++)
	//	this->treecolors1->InsertNextTupleValue(color);
	//int num_point_id = 0;
	//int counter = this->num_rows;
	//for(int i = 0; i < this->tree1.treeid.size(); i++)
	//{
	//	double xoffset;
	//	xoffset = this->tree1.coordinates[counter++][0];
	//	for(int j = 1; j < this->tree1.treeid[i].size(); j++)
	//	{
	//		double p1[3];	
	//		double p2[3];
	//		p1[0] = this->tree1.coordinates[this->tree1.treeid[i][j]][0];
	//		p1[1] = this->tree1.coordinates[this->tree1.treeid[i][j]][1];
	//		p1[2] = this->tree1.coordinates[this->tree1.treeid[i][j]][2];				
	//		p2[0] = xoffset;
	//		p2[1] = p1[1];
	//		p2[2] = p1[2];
	//		this->treepoints1->InsertNextPoint(p1);
	//		this->treepoints1->InsertNextPoint(p2);

	//		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	//		line->GetPointIds()->SetId(0, 2*(num_point_id + j - 1));
	//		line->GetPointIds()->SetId(1, 2*(num_point_id + j - 1) + 1);
	//		this->treelines1->InsertNextCell(line);
	//	}
	//	num_point_id += (this->tree1.treeid[i].size() - 1);
	//}

	//counter = this->num_rows;
	//for(int i = 0; i < this->tree1.treeid.size(); i++)
	//{
	//	double xoffset = this->tree1.coordinates[counter++][0];
	//	double p1[3];
	//	double p2[3];
	//	p1[0] = xoffset;
	//	p1[1] = this->tree1.coordinates[this->tree1.treeid[i][1]][1];
	//	p1[2] = 0;
	//	p2[0] = xoffset;
	//	p2[1] = this->tree1.coordinates[this->tree1.treeid[i][this->tree1.treeid[i].size() - 1]][1];
	//	p2[2] = 0;

	//	this->treepoints1->InsertNextPoint(p1);
	//	this->treepoints1->InsertNextPoint(p2);
	//	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
	//	line->GetPointIds()->SetId(0, 2*(num_nodes - 1 + i));
	//	line->GetPointIds()->SetId(1, 2*(num_nodes - 1 + i) + 1);
	//	this->treelines1->InsertNextCell(line);
	//}
/////////////////////////////////////////////////////////////////////
	this->treelinesPolyData1 =vtkSmartPointer<vtkPolyData>::New();
	this->treelinesPolyData1->SetPoints(this->treepoints1);
	this->treelinesPolyData1->SetLines(this->treelines1);
	this->treelinesPolyData1->GetCellData()->SetScalars(this->treecolors1);
	
	this->treemapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->treemapper1->SetInput(this->treelinesPolyData1);
	this->treemapper1->SetScalarRange(0, 255);

	this->treeactor1 = vtkSmartPointer<vtkActor>::New();
	this->treeactor1->SetMapper(this->treemapper1);

	this->view->GetRenderer()->AddActor(this->treeactor1);
	this->view->Render();
}

void BiHeatmap::showTree2()
{
	int num_nodes = this->tree2.coordinates.size();

	this->treecolors2 = vtkSmartPointer<vtkUnsignedCharArray>::New();
	this->treecolors2->SetNumberOfComponents(3);
	this->treecolors2->SetName("treecolors2");
	unsigned char color[3] = {0, 0, 0};
	for(int i = 0; i < num_nodes - 1; i++)
		this->treecolors2->InsertNextTupleValue(color);

	int counter = this->num_cols;
	int num_pp_id = 0;
	this->treepoints2 = vtkSmartPointer<vtkPoints>::New();	
	this->treelines2 = vtkSmartPointer<vtkCellArray>::New();
	for(int i = 0; i < this->tree2.treeid.size(); i++)
	{
		double pp[3];
		pp[0] = this->tree2.coordinates[counter][0];
		pp[1] = this->tree2.coordinates[counter][1];
		pp[2] = this->tree2.coordinates[counter++][2];
		this->treepoints2->InsertNextPoint(pp);
		for(int j = 1; j < this->tree2.treeid[i].size(); j++)
		{
			double ps[3];			
			ps[0] = this->tree2.coordinates[this->tree2.treeid[i][j]][0];
			ps[1] = this->tree2.coordinates[this->tree2.treeid[i][j]][1];
			ps[2] = this->tree2.coordinates[this->tree2.treeid[i][j]][2];
			this->treepoints2->InsertNextPoint(ps);

			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			line->GetPointIds()->SetId(0, num_pp_id);
			line->GetPointIds()->SetId(1, num_pp_id + j);
			this->treelines2->InsertNextCell(line);
		}
		num_pp_id += this->tree2.treeid[i].size();
	}

	this->treelinesPolyData2 =vtkSmartPointer<vtkPolyData>::New();
	this->treelinesPolyData2->SetPoints(this->treepoints2);
	this->treelinesPolyData2->SetLines(this->treelines2);
	this->treelinesPolyData2->GetCellData()->SetScalars(this->treecolors2);
	
	this->treemapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->treemapper2->SetInput(this->treelinesPolyData2);
	this->treemapper2->SetScalarRange(0, 255);

	this->treeactor2 = vtkSmartPointer<vtkActor>::New();
	this->treeactor2->SetMapper(this->treemapper2);

	this->view->GetRenderer()->AddActor(this->treeactor2);
	this->view->Render();
}

void BiHeatmap::setDataForTree1(Level_id levels1)
{
	this->levels1 = levels1;
	this->treedata1 = true;
}

void BiHeatmap::setDataForTree2(Level_id levels2)
{
	this->levels2 = levels2;
	this->treedata2 = true;
}

void BiHeatmap::creatDataForTree1()
{
	double xoffset = 0.15;
	std::vector<double > tempco;
	for(int i = 0; i < this->num_rows; i++)
	{
		tempco.clear();
		double var = (double)rowMapFromOriginalToReorder.find(i)->second;
		tempco.push_back(-0.5);
		tempco.push_back(-0.5 + var*(1.0/this->num_rows) + 0.5/this->num_rows);
		tempco.push_back(0);
		this->tree1.coordinates.push_back(tempco);
	}

	int counter = this->num_rows;
	int num_before_two_rows = this->num_rows;
	for(int i = 1; i <= this->levels1.numm_levels; i++)
	{
		for(int j = 0; j < this->levels1.level_id[i].size(); j++)
		{
			std::vector<int > tempid;
			tempid.push_back(counter++);
			for(int k = 0; k < this->levels1.level_id[i][j].size(); k++)
			{
				if(i == 1)
					tempid.push_back(this->levels1.level_id[i][j][k]);	
				else
					tempid.push_back(this->levels1.level_id[i][j][k] + num_before_two_rows);
			}
			this->tree1.treeid.push_back(tempid);

			tempco.clear();
			int id1 = tempid[1];
			int id2 = tempid[tempid.size() - 1];
			double mid = (this->tree1.coordinates[id1][1] + this->tree1.coordinates[id2][1])/2;
			tempco.push_back(-0.5 - xoffset * i);
			tempco.push_back(mid);
			tempco.push_back(0);
			this->tree1.coordinates.push_back(tempco);
		}
		if(i >= 2)
			num_before_two_rows += this->levels1.level_id[i - 1].size();
	}
}

void BiHeatmap::creatDataForTree2()
{
	double yoffset = 0.15;
	std::vector<double > tempco;
	for(int i = 0; i < this->num_cols; i++)
	{
		tempco.clear();
		double var = (double)columnMapFromOriginalToReorder.find(i)->second;
		tempco.push_back(-0.5 + var*(1.0/this->num_cols) + 0.5/this->num_cols);
		tempco.push_back(0.5);	
		tempco.push_back(0);
		this->tree2.coordinates.push_back(tempco);
	}

	int counter = this->num_cols;
	int num_before_two_rows = this->num_cols;
	for(int i = 1; i <= this->levels2.numm_levels; i++)
	{
		for(int j = 0; j < this->levels2.level_id[i].size(); j++)
		{
			std::vector<int > tempid;
			tempid.push_back(counter++);
			for(int k = 0; k < this->levels2.level_id[i][j].size(); k++)
			{
				if(i == 1)
					tempid.push_back(this->levels2.level_id[i][j][k]);	
				else
					tempid.push_back(this->levels2.level_id[i][j][k] + num_before_two_rows);
			}
			this->tree2.treeid.push_back(tempid);

			tempco.clear();
			int id1 = tempid[1];
			int id2 = tempid[tempid.size() - 1];
			double mid = (this->tree2.coordinates[id1][0] + this->tree2.coordinates[id2][0])/2;			
			tempco.push_back(mid);
			tempco.push_back(0.5 + yoffset * i);
			tempco.push_back(0);
			this->tree2.coordinates.push_back(tempco);
		}
		if(i >= 2)
			num_before_two_rows += this->levels2.level_id[i - 1].size();
	}
}
void BiHeatmap::GetSelecectedIDs()
{
	std::cout<<"get selected"<<endl;
	std::set<long int> selectedIDs2 = this->Selection2->getSelections();
	std::set<long int> selectedIDs1 = this->Selection->getSelections();	

	if(selectedIDs2.empty ()&& selectedIDs1.empty())
	{
		if(this->removeActorflag == true)
		{
			if(this->dragLineFlag1)
				this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor1);
			if(this->dragLineFlag2)
				this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor2);

			for(int i = 0; i < this->num_selection_area; i++)
				this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
			this->removeActorflag = false;
			this->num_selection_area = 0;

			if(this->dragLineFlag1)
				this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor1);
			if(this->dragLineFlag2)
				this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor2);
		}
		return;
	}
	std::set<long int>::iterator iter1 = selectedIDs1.begin();
	std::set<long int>::iterator iter2 = selectedIDs2.begin();
	vtkSmartPointer<vtkIdTypeArray> getcellids = vtkSmartPointer<vtkIdTypeArray>::New();
	getcellids->SetNumberOfComponents(1);

	int num1 = selectedIDs1.size();
	int num2 = selectedIDs2.size();

	std::vector<int > IDs1;
	std::vector<int > IDs2;
	IDs1.resize(num1);
	IDs2.resize(num2);

	int count1 = 0;
	int count2 = 0;
	//////////////////////////////////
	int min1 = this->num_rows;
	int max1 = 0;
	//////////////////////////////
	while(iter1 != selectedIDs1.end())
	{
		int index1 = *iter1;
		int var = indMapFromVertexToInd.find(index1)->second;
		int id1 = rowMapFromOriginalToReorder.find(var)->second;
		IDs1[count1++] = id1;
		iter1++;
		///////////////////////////////////////
		if (id1 < min1)
		{
			min1 = id1;
		}
		if (id1 > max1)
		{
			max1 = id1;
		}
		////////////////////////////////
	}

	//////////////////////////////////
	int min2 = this->num_cols;
	int max2 = 0;
	//////////////////////////////
	while(iter2 != selectedIDs2.end())
	{
		int index2 = *iter2;
		int id2 = columnMapFromOriginalToReorder.find(index2)->second;
		IDs2[count2++] = id2;
		iter2++;
			///////////////////////////////////////
		if (id2 < min2)
		{
			min2 = id2;
		}
		if (id2 > max2)
		{
			max2 = id2;
		}
		////////////////////////////////
	}

	if( num1 == 0 && num2 != 0)
	{
		num1 = this->num_rows;
		IDs1.resize(num1);
		for( int i = 0; i < this->num_rows; i++)
			IDs1[i] = i;
		///////////////////
		min1 = 0;
		max1 = this->num_rows - 1;
		/////////////////
	}

	if( num2 == 0 && num1 != 0)
	{
		num2 = this->num_cols;
		IDs2.resize(num2);
		for( int i = 0; i < this->num_cols; i++)
			IDs2[i] = i;
		///////////////////
		min2 = 0;
		max2 = this->num_cols - 1;
		/////////////////
	}

	for(int i = 0; i<num1; i++)
		for(int j = 0; j<num2; j++)
			getcellids->InsertNextValue( (IDs1[i])*(this->num_cols) + IDs2[j]);

	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(getcellids);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);
	 
	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInput(0, this->plane->GetOutput());
	extractSelection->SetInput(1, selection);
	extractSelection->Update();
	 
	vtkSmartPointer<vtkUnstructuredGrid> selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extractSelection->GetOutput());
	
	vtkSmartPointer<vtkDataSetMapper> selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selectedMapper->SetInputConnection(selected->GetProducerPort());
	 
	vtkSmartPointer<vtkActor> selectedActor = vtkSmartPointer<vtkActor>::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->GetProperty()->EdgeVisibilityOn();
	selectedActor->GetProperty()->SetEdgeColor(0.6 ,0.19, 0.8);
	selectedActor->GetProperty()->SetLineWidth(0.5);

/////////////////////////////////////////////////////////////////////////////
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];

	p1[0] = -0.5 + min2 * (1.0/this->num_cols);
	p1[1] = -0.5 + (max1 + 1) * (1.0/this->num_rows);
	p1[2] = 0;

	p2[0] = -0.5 + (max2 + 1) * (1.0/this->num_cols);
	p2[1] = -0.5 + (max1 + 1) * (1.0/this->num_rows);
	p2[2] = 0;

	p3[0] = -0.5 + min2 * (1.0/this->num_cols);
	p3[1] = -0.5 + min1 * (1.0/this->num_rows);
	p3[2] = 0;

	p4[0] = -0.5 + (max2 + 1) * (1.0/this->num_cols);
	p4[1] = -0.5 + min1 * (1.0/this->num_rows);
	p4[2] = 0;
	vtkSmartPointer<vtkLineSource> lineSource1 = vtkSmartPointer<vtkLineSource>::New();
	lineSource1->SetPoint1(p1);
	lineSource1->SetPoint2(p2);
	vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper1->SetInputConnection(lineSource1->GetOutputPort());
	vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
	actor1->GetProperty()->SetColor(0.6,0.19,0.8); 
	actor1->GetProperty()->SetLineWidth(3);
	actor1->SetMapper(mapper1);

	vtkSmartPointer<vtkLineSource> lineSource2 = vtkSmartPointer<vtkLineSource>::New();
	lineSource2->SetPoint1(p2);
	lineSource2->SetPoint2(p4);
	vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper2->SetInputConnection(lineSource2->GetOutputPort());
	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->GetProperty()->SetColor(0.6,0.19,0.8); 
	actor2->GetProperty()->SetLineWidth(3);
	actor2->SetMapper(mapper2);

	vtkSmartPointer<vtkLineSource> lineSource3 = vtkSmartPointer<vtkLineSource>::New();
	lineSource3->SetPoint1(p3);
	lineSource3->SetPoint2(p4);
	vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper3->SetInputConnection(lineSource3->GetOutputPort());
	vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
	actor3->GetProperty()->SetColor(0.6,0.19,0.8); 
	actor3->GetProperty()->SetLineWidth(3);
	actor3->SetMapper(mapper3);

	vtkSmartPointer<vtkLineSource> lineSource4 = vtkSmartPointer<vtkLineSource>::New();
	lineSource4->SetPoint1(p1);
	lineSource4->SetPoint2(p3);
	vtkSmartPointer<vtkPolyDataMapper> mapper4 = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper4->SetInputConnection(lineSource4->GetOutputPort());
	vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
	actor4->GetProperty()->SetColor(0.6,0.19,0.8); 
	actor4->GetProperty()->SetLineWidth(3);
	actor4->SetMapper(mapper4);
////////////////////////////////////////////////////////////////////////

	if(this->dragLineFlag1)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor1);
	if(this->dragLineFlag2)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor2);

	if(this->removeActorflag == true)
	{
		for(int i = 0; i < this->num_selection_area; i++)
			this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor (this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastActor());
	}
	
	if(this->localselection == false)
	{
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
		this->num_selection_area = 1;
	}
	else
	{
		this->view->GetRenderer()->AddActor(actor1);
		this->view->GetRenderer()->AddActor(actor2);
		this->view->GetRenderer()->AddActor(actor3);
		this->view->GetRenderer()->AddActor(actor4);
		this->num_selection_area = 4;
	}
	this->removeActorflag = true;
	this->localselection = false;

	if(this->dragLineFlag1)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor1);
	if(this->dragLineFlag2)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(dragLineActor2);
	this->view->Render();
}

void BiHeatmap::drawPoints()
{
	int num_nodes1 = this->tree1.coordinates.size();
	int num_nodes2 = this->tree2.coordinates.size();
	int num_nodes = num_nodes1 + num_nodes2;

	this->graph_Layout = vtkSmartPointer<vtkMutableUndirectedGraph>::New();
	this->points = vtkSmartPointer<vtkPoints>::New();
	this->vertex = vtkSmartPointer<vtkIdTypeArray>::New();
	this->vertex->SetNumberOfValues (num_nodes);

	for(int i = 0; i < num_nodes1; i++)
    {
		this->vertex->SetValue (i,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->tree1.coordinates[i][0], this->tree1.coordinates[i][1], this->tree1.coordinates[i][2]);
	}
	
	for(int i = 0; i < num_nodes2; i++)
    {
		this->vertex->SetValue (i + num_nodes1,graph_Layout->AddVertex());
		this->points->InsertNextPoint(this->tree2.coordinates[i][0], this->tree2.coordinates[i][1], this->tree2.coordinates[i][2]);
	}
    this->graph_Layout->SetPoints(this->points);

	this->vertexColors = vtkSmartPointer<vtkIntArray>::New();   
    this->vertexColors->SetNumberOfComponents(1);
    this->vertexColors->SetName("pointColor1");	

	this->vetexlut = vtkSmartPointer<vtkLookupTable>::New();
	vetexlut->SetNumberOfTableValues(num_nodes);
    for(int i = 0; i < num_nodes; i++)
		this->vetexlut->SetTableValue(i, 0.5, 0.5,0.5);
    this->vetexlut->Build();
   
	this->pointscales1 = vtkSmartPointer<vtkFloatArray>::New(); /// scales for vertex size
    this->pointscales1->SetNumberOfComponents(1);
	this->pointscales1->SetName("pointScales1");

    for(int i = 0; i < num_nodes; i++)
    {
		this->vertexColors->InsertNextValue(i);
		this->pointscales1->InsertNextValue(1);
    }

	this->graph_Layout->GetVertexData()->AddArray(this->pointscales1);
	this->graph_Layout->GetVertexData()->AddArray(this->vertexColors);

    this->view->AddRepresentationFromInput(graph_Layout);
	this->view->SetLayoutStrategy("Pass Through");	
    this->view->ScaledGlyphsOn();
	this->view->SetScalingArrayName("pointScales1");
	this->view->ColorVerticesOn();
	this->view->SetVertexColorArrayName("pointColor1");
	vtkRenderedGraphRepresentation::SafeDownCast(this->view->GetRepresentation()) ->SetGlyphType(vtkGraphToGlyphs::CIRCLE);

	this->theme = vtkSmartPointer<vtkViewTheme>::New();	
	this->theme->SetPointLookupTable(this->vetexlut);
	this->view->ApplyViewTheme(theme);  

	this->selectionCallbackleft =vtkSmartPointer<vtkCallbackCommand>::New();
	this->selectionCallbackleft->SetClientData(this);
	this->selectionCallbackleft->SetCallback(SelectionCallbackFunctionLeftButton);
    this->view->GetRepresentation()->GetAnnotationLink()->AddObserver("AnnotationChangedEvent", this->selectionCallbackleft);
}
void BiHeatmap::SelectionCallbackFunctionLeftButton(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	BiHeatmap* heatmap = (BiHeatmap*)clientData;

	vtkSelectionNode* vertices = NULL;
	vtkSelectionNode* edges = NULL;
	vtkSelectionNode* cells = NULL;

    if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(0);
        }
    else if(selection->GetNode(0)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(0);
        }
 
    if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
        {
        vertices = selection->GetNode(1);
        }
    else if(selection->GetNode(1)->GetFieldType() == vtkSelectionNode::EDGE)
        {
        edges = selection->GetNode(1);
        }

	if( vertices != NULL)
	{
		vtkIdTypeArray* vertexList = vtkIdTypeArray::SafeDownCast(vertices->GetSelectionList());
	
		std::set<long int> IDs;
		if(vertexList->GetNumberOfTuples() > 0)
		{

			for( vtkIdType i = 0; i < vertexList->GetNumberOfTuples(); i++)
			{
				long int value = vertexList->GetValue(i);
				IDs.insert(value);
			}
		}
		heatmap->setSelectIds(IDs);
	}
}

void BiHeatmap::SelectionCallbackFunctionRightButtonDown(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	BiHeatmap* heatmap = (BiHeatmap*)clientData;
	int* pos = heatmap->view->GetInteractor()->GetEventPosition();
 
	vtkCellPicker *cell_picker = (vtkCellPicker *)heatmap->view->GetInteractor()->GetPicker();

	cell_picker->Pick(pos[0], pos[1], 0, heatmap->view->GetRenderer());
	double* worldPosition = cell_picker->GetPickPosition();
 
	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, heatmap->view->GetRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();

		if(cellpicker->GetCellId() != -1)
			heatmap->rightbuttonid1 = cellpicker->GetCellId();
	}
	else if(worldPosition[0]<-0.5 && worldPosition[1]<0.5)
	{
		//heatmap->addDragLine1(worldPosition);
	}	
	else if(worldPosition[0]>-0.5 && worldPosition[1]>0.5)
	{
		//heatmap->addDragLine2(worldPosition);
	}
}

void BiHeatmap::SelectionCallbackFunctionRightButtonUp(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	BiHeatmap* heatmap = (BiHeatmap*)clientData;
	int* pos = heatmap->view->GetInteractor()->GetEventPosition();
 
	vtkCellPicker *cell_picker = (vtkCellPicker *)heatmap->view->GetInteractor()->GetPicker();
 
	// Pick from this location.
	cell_picker->Pick(pos[0], pos[1], 0, heatmap->view->GetRenderer());
	double* worldPosition = cell_picker->GetPickPosition();

	if((worldPosition[0]<=0.5) && (worldPosition[0]>=-0.5) && (worldPosition[1]<=0.5) && (worldPosition[0]>=-0.5))
	{
		vtkSmartPointer<vtkCellPicker> cellpicker = vtkSmartPointer<vtkCellPicker>::New();
		cellpicker->SetTolerance(0.0005);
 
		// Pick from this location.
		cellpicker->Pick(pos[0], pos[1], 0, heatmap->view->GetRenderer());
 
		double* worldPosition = cellpicker->GetPickPosition();
		if(cellpicker->GetCellId() != -1)
		{
			heatmap->rightbuttonid2 = cellpicker->GetCellId();
			heatmap->computeselectedcells();
			heatmap->setselectedCellIds();
		}
	}
}

void BiHeatmap::HandleKeyPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	BiHeatmap* heatmapWin = (BiHeatmap*)clientData;
	char key = heatmapWin->view->GetInteractor()->GetKeyCode();
	int size = 0;
	switch (key)
	{
	case 'c':
		heatmapWin->continueselect = true;
		break;
	case 'i':
		heatmapWin->intersectionselect = true;
		break;
	case 'r':
		heatmapWin->continueselect = false;
		heatmapWin->intersectionselect = false;
		break;
	case 'd':
		size = heatmapWin->Selection->getSelections().size();
		if( size > 0)
		{
			std::cout<< size<< " items deleted! Please rerun SPD analysis!"<<endl;
			heatmapWin->Selection->DeleteCurrentSelectionInTable();
		}
		else
		{
			std::cout<< "No items have been selected!"<<endl;
		}
		break;
	default:
		break;
	}
}

void BiHeatmap::setSelectIds(std::set<long int>& IDs)
{
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;
	std::set<long int>::iterator it;
	long int id;

	this->resetTree1();
	this->resetTree2();

	if(IDs.size() > 0)
	{
		this->localselection = true;
		for(it = IDs.begin(); it != IDs.end(); it++)
		{
			id = *it;
			if(id < this->tree1.coordinates.size())
				reselectIds1(selectedIDs1, id);
			else
				reselectIds2(selectedIDs2, id - this->tree1.coordinates.size());
		}
	}

	if(selectedIDs1.size() > 0)
	{
		for(int i = 0; i < this->num_cols; i++)
			selectedIDs2.insert(i);

		std::cout<<"I selected from tree1"<<std::endl;
		this->updataTree1();
		this->Selection2->select(selectedIDs2);
		this->Selection->select(selectedIDs1);
	}
	else if(selectedIDs2.size() > 0)
	{
		for(int i = 0; i<this->num_rows; i++)
			selectedIDs1.insert( indMapFromIndToVertex[i]);
		std::cout<<"I selected from tree2"<<std::endl;
		this->updataTree2();
		this->Selection2->select(selectedIDs2);
		this->Selection->select(selectedIDs1);
	}
	else
	{
		std::cout<<"I selected nothing"<<std::endl;
		this->Selection2->clear();
		this->Selection->clear();
	}
}

void BiHeatmap::reselectIds1(std::set<long int>& selectedIDs, long int id)
{
	if(id < this->num_rows)
	{
		std::cout<<indMapFromIndToVertex[id]<<"..."<<std::endl;
		selectedIDs.insert(indMapFromIndToVertex[id]);
	}
	else
	{
		int num_line = 0;
		for(int i = 0; i < this->tree1.treeid.size(); i++)
		{
			if(this->tree1.treeid [i][0] == (int)id)
			{
				for(int j = 1; j < this->tree1.treeid[i].size(); j++)	
				{
					this->treecolors1->SetValue((num_line + j - 1)*3, 150);
					this->treecolors1->SetValue((num_line + j - 1)*3 + 1, 50);
					this->treecolors1->SetValue((num_line + j - 1)*3 + 2, 204);
					this->reselectIds1(selectedIDs, this->tree1.treeid[i][j]);
				}
				break;
			}
			num_line += (this->tree1.treeid[i].size() - 1);
		}
	}
}


void BiHeatmap::reselectIds2(std::set<long int>& selectedIDs2, long int id)
{
	if(id < this->num_cols)
	{
		std::cout<<id<<"..."<<std::endl;
		selectedIDs2.insert( id);
	}
	else
	{
		int num_line = 0;
		for(int i = 0; i < this->tree2.treeid.size(); i++)
		{
			if(this->tree2.treeid [i][0] == (int)id)
			{
				for(int j = 1; j < this->tree2.treeid[i].size(); j++)
				{
					this->treecolors2->SetValue((num_line + j - 1)*3, 150);
					this->treecolors2->SetValue((num_line + j - 1)*3 + 1, 50);
					this->treecolors2->SetValue((num_line + j - 1)*3 + 2, 204);
					this->reselectIds2(selectedIDs2, this->tree2.treeid[i][j]);
				}
				break;
			}
			num_line += (this->tree2.treeid[i].size() - 1);
		}
	}
}

void BiHeatmap::computeselectedcells()
{
	this->r1 = this->rightbuttonid1 / this->num_cols;
	this->r2 = this->rightbuttonid2 / this->num_cols;
	this->c1 = this->rightbuttonid1 % this->num_cols;
	this->c2 = this->rightbuttonid2 % this->num_cols;

	std::cout<<r1<<endl;
	std::cout<<r2<<endl;
	this->cellids = vtkSmartPointer<vtkIdTypeArray>::New();
	this->cellids->SetNumberOfComponents(1);
	for(int i = 0; i <= r1 - r2; i++)
	{
		for(int j = 0; j <= c2 - c1; j++)
			this->cellids->InsertNextValue(this->rightbuttonid2 - j + this->num_cols * i);
	}
}

void BiHeatmap::setselectedCellIds()
{
	if(this->treedata1 == true && this->treedata2 == true)
	{
		this->resetTree1();
		this->resetTree2();
	}
	std::set<long int> selectedIDs1;
	std::set<long int> selectedIDs2;
	for(int i = r2; i <= r1; i++)
		selectedIDs1.insert( indMapFromIndToVertex[ this->order1[i]]);
	for(int j = c1; j <= c2; j++)	
		selectedIDs2.insert(this->order2[j]);
		
	this->localselection = true;
	this->Selection2->select(selectedIDs2);
	this->Selection->select(selectedIDs1);
}

void BiHeatmap::updataTree1()
{
	this->treemapper1->ScalarVisibilityOn();
	this->treelinesPolyData1->Modified();
	this->treelinesPolyData1->Update();
	this->treemapper1->Modified();
	this->treemapper1->Update();
	this->treeactor1->Modified();
	this->view->Render();
}

void BiHeatmap::updataTree2()
{
	this->treemapper2->ScalarVisibilityOn();
	this->treelinesPolyData2->Modified();
	this->treelinesPolyData2->Update();
	this->treemapper2->Modified();
	this->treemapper2->Update();
	this->treeactor2->Modified();
	this->view->Render();
}

void BiHeatmap::resetTree1()
{
	for(int i = 0; i < this->treecolors1->GetSize(); i++)
		this->treecolors1->SetValue(i, 0);
	this->updataTree1();
}

void BiHeatmap::resetTree2()
{
	for(int i = 0; i < this->treecolors2->GetSize(); i++)
		this->treecolors2->SetValue(i, 0);
	this->updataTree2();
}

void BiHeatmap::addDragLine1(double* worldPosition)
{
	if(this->dragLineFlag1)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor1);

	double p1[3];
	double p2[3];
	p1[0]=worldPosition[0];
	p1[1]=-0.75;
	p1[2]=0;
	p2[0]=worldPosition[0];
	p2[1]=0.75;
	p2[2]=0;

	dragLineSource1 = vtkSmartPointer<vtkLineSource>::New();
	dragLineSource1->SetPoint1(p1);
	dragLineSource1->SetPoint2(p2);

	dragLineMapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
	dragLineMapper1->SetInputConnection(dragLineSource1->GetOutputPort());

	dragLineActor1 = vtkSmartPointer<vtkActor>::New();
	dragLineActor1->SetMapper(dragLineMapper1);
	dragLineActor1->GetProperty()->SetColor(0.5,0.7,0);
	dragLineActor1->GetProperty()->SetLineWidth(1.5);
	dragLineActor1->DragableOn();
	this->view->GetRenderer()->AddActor(dragLineActor1);
	this->dragLineFlag1 = true;
	this->view->Render();
}

void BiHeatmap::addDragLine2(double* worldPosition)
{
	if(this->dragLineFlag2)
		this->view->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(dragLineActor2);

	double p1[3];
	double p2[3];
	p1[0]=-0.75;
	p1[1]=worldPosition[1];
	p1[2]=0;
	p2[0]=0.75;
	p2[1]=worldPosition[1];;
	p2[2]=0;

	dragLineSource2 = vtkSmartPointer<vtkLineSource>::New();
	dragLineSource2->SetPoint1(p1);
	dragLineSource2->SetPoint2(p2);

	dragLineMapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
	dragLineMapper2->SetInputConnection(dragLineSource2->GetOutputPort());

	dragLineActor2 = vtkSmartPointer<vtkActor>::New();
	dragLineActor2->SetMapper(dragLineMapper2);
	dragLineActor2->GetProperty()->SetColor(0.5,0.7,0);
	dragLineActor2->GetProperty()->SetLineWidth(1.5);
	dragLineActor2->DragableOn();
	this->view->GetRenderer()->AddActor(dragLineActor2);
	this->dragLineFlag2 = true;
	this->view->Render();
}

void BiHeatmap::showFeatureNames()
{
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->SetNumberOfCells(1);

	vtkSmartPointer<vtkDoubleArray> orient = vtkSmartPointer<vtkDoubleArray>::New();
	orient->SetNumberOfComponents(1);
	orient->SetName("orientation");

	vtkSmartPointer<vtkStringArray> label = vtkSmartPointer<vtkStringArray>::New();
	label->SetNumberOfComponents(1);
	label->SetName("label");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for(int i = 0; i < this->num_cols; i++)
    {
		pts->InsertNextPoint(this->tree2.coordinates[i][0], -0.5, 0);
	}

	for(int i=0; i<this->num_cols;i++)
	{
		verts->InsertNextCell(1);
		verts->InsertCellPoint(i);
		orient->InsertNextValue(315.0);
		vtkIdType id = i+1  ;
		label->InsertNextValue(this->table->GetColumn(id)->GetName());
	}

	pd->SetPoints(pts);
	pd->SetVerts(verts);
	pd->GetPointData()->AddArray(label);
	pd->GetPointData()->AddArray(orient);

	vtkSmartPointer<vtkPointSetToLabelHierarchy> hier = vtkSmartPointer<vtkPointSetToLabelHierarchy>::New();
	hier->SetInput(pd);
	hier->SetOrientationArrayName("orientation");
	hier->SetLabelArrayName("label");
	hier->GetTextProperty()->SetColor(0.0, 0.0, 0.0);
  
	vtkSmartPointer<vtkLabelPlacementMapper> lmapper = vtkSmartPointer<vtkLabelPlacementMapper>::New();
	lmapper->SetInputConnection(hier->GetOutputPort());

	vtkSmartPointer<vtkQtLabelRenderStrategy> strategy = vtkSmartPointer<vtkQtLabelRenderStrategy>::New();
	lmapper->SetRenderStrategy(strategy);
	lmapper->SetShapeToNone();
	lmapper->SetBackgroundOpacity(0.0);
	lmapper->SetMargin(0);

	vtkSmartPointer<vtkActor2D> lactor = vtkSmartPointer<vtkActor2D>::New();
	lactor->SetMapper(lmapper);

	vtkSmartPointer<vtkPolyDataMapper> rmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	rmapper->SetInput(pd);

	vtkSmartPointer<vtkActor> ractor = vtkSmartPointer<vtkActor>::New();
	ractor->SetMapper(rmapper);

	this->view->GetRenderer()->AddActor(lactor);
	this->view->GetRenderer()->AddActor(ractor);
}

void BiHeatmap::WriteFile(const char *filename1)
{
	FILE *fp1 = fopen(filename1,"w");
	for(int i=0; i<this->num_rows; i++)
	{
		for(int j=0; j<this->num_cols; j++)
		{
			fprintf(fp1,"%f",this->data[i][j]);
			fprintf(fp1,"\t");
		}				
		fprintf(fp1,"\n");
	}
	fclose(fp1);
}