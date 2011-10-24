#include "HeatmapWindow.h"

Heatmap::Heatmap(QWidget *parent)
: QMainWindow(parent)
{
	this->mainQTRenderWidget;
	this->aPlane = vtkSmartPointer<vtkPlaneSource>::New();
	this->cellData = vtkSmartPointer<vtkFloatArray>::New();
	this->lookuptable = vtkSmartPointer<vtkLookupTable>::New();
	this->mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	this->actor = vtkSmartPointer<vtkActor>::New();
	this->view = vtkSmartPointer<vtkRenderView>::New();
	this->theme = vtkSmartPointer<vtkViewTheme>::New();
	this->cellColors = vtkSmartPointer<vtkIntArray>::New();

}

Heatmap::~Heatmap()
{
}

void Heatmap::setDataForHeatmap(double** features, int* optimalleaforder1, int* optimalleaforder2,int num_samples, int num_features)
{
	this->mapdata = features;
	this->Optimal_Leaf_Order1 = optimalleaforder1;
	this->Optimal_Leaf_Order2 = optimalleaforder2;
	this->num_samples = num_samples;
	this->num_features = num_features;
}

void Heatmap::creatDataForHeatmap()
{
	vector<double > temp;
	temp.resize(num_features);
	///////
	double** tempdata;
	tempdata = new double*[this->num_samples];
	for(int i = 0; i < this->num_samples; i++)
		tempdata[i] = new double[this->num_features];
	/////

	for(int i = 0; i < this->num_samples; i++)
	{
		double mean = 0.0; 
		double std = 0.0;
		double sum = 0.0;
		for(int j = 0; j < this->num_features; j++)
		{
			temp[j] = mapdata[i][Optimal_Leaf_Order2[j]];
			mean += temp[j];
		}

		mean = mean / this->num_features;

		for(int j = 0; j < num_features; j++)
			sum += (temp[j] - mean)*(temp[j] - mean);

		std = sqrt(sum/this->num_features);

		for(int j = 0; j < num_features; j++)
			tempdata[i][j] = (temp[j] - mean)/std;
	}
	for(int i = 0; i < this->num_samples; i++)
		mapdata[this->num_samples - i - 1] = tempdata[Optimal_Leaf_Order1[i]]; 
//////////////////////////////////////////////////////////for debug
	const char* filename = "heatmapdata";
	FILE *fp = fopen(filename,"w");
	for(int i=0; i<this->num_samples; i++)
	{
		for(int j=0; j<this->num_features; j++)
			fprintf(fp,"%14.7e\t",mapdata[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
//////////////////////////////////////////////
}

void Heatmap::setModels(vtkSmartPointer<vtkTable> table = NULL, ObjectSelection * sels = NULL, ObjectSelection * sels2 = NULL)
{
}

void Heatmap::showGraph()
{	
    this->aPlane->SetXResolution(this->num_features);
    this->aPlane->SetYResolution(this->num_samples);

	int index = 0;
	//cellColors->SetNumberOfComponents(1);
    //cellColors->SetName("Color");

	for (int i = 0; i < this->num_samples; i++)
    {
		for(int j = 0; j < this->num_features; j++)
		{
			cellData->InsertNextValue(index++);
			//cellColors->InsertNextValue(index++);
		}
    }
	
	this->lookuptable->SetNumberOfTableValues(this->num_samples*this->num_features);
	this->lookuptable->Build();

	int k = 0;
	for(int i = 0; i < this->num_samples; i++)
	{
		for(int j = 0; j < this->num_features; j++)
		{
			if(mapdata[num_samples - i - 1][j] <= 0)
				lookuptable->SetTableValue(k++, 0, (1 - mapdata[num_samples - i - 1][j])/5, 0);
			else
				lookuptable->SetTableValue(k++, (mapdata[num_samples - i - 1][j])/5, 0, 0);
		}
	}

	this->aPlane->Update(); // Force an update so we can set cell data
	this->aPlane->GetOutput()->GetCellData()->SetScalars(cellData);
	//this->aPlane->GetOutput()->GetCellData()->AddArray(cellColors);
 
	/*this->view->AddRepresentationFromInput(this->aPlane->GetOutput());
	this->theme->SetCellLookupTable(lookuptable);
	this->view->ApplyViewTheme(theme);*/
	// Setup actor and mapper
	this->mapper->SetInputConnection(aPlane->GetOutputPort());
	this->mapper->SetScalarRange(0, this->num_samples*this->num_features - 1);
	this->mapper->SetLookupTable(lookuptable);
	this->actor->SetMapper(mapper);
	
 
	// Setup render window, renderer, and interactor
	//vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	//vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	//renderWindow->AddRenderer(renderer);
	this->mainQTRenderWidget.SetRenderWindow(view->GetRenderWindow());
	this->mainQTRenderWidget.resize(600, 600);
	this->mainQTRenderWidget.show();

	//vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//renderWindowInteractor->SetRenderWindow(renderWindow);
	//renderer->AddActor(actor);
	this->view->GetRenderer()->AddActor(actor);
	this->view->GetRenderer()->SetBackground(1,1,1);
	//renderWindow->Render();
	//renderWindowInteractor->Start();
	this->view->Render();
	this->view->GetInteractor()->Start();
}