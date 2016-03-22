
#include "FTKRenderWindow.h"
#include <QInputDialog>

FTKRenderWindow::FTKRenderWindow(QWidget *parent)
: QMainWindow(parent)
{
	this->QVTK = 0;
	table = NULL;
	selection = NULL;
	selectedFeature = 1;
	selectedColorCode = 1;
	radius = 30;
	discreteColorMap.clear();	
	this->setupUI();
	createDiscreteColorMap();
}

FTKRenderWindow::~FTKRenderWindow()
{
	if(this->QVTK)
	{
		delete this->QVTK;
	}
}

void FTKRenderWindow::setupUI(void)
{
	//QWidget *centralWidget = new QWidget();

	featureMenu = new QMenu(tr("Select Feature"));
	connect(featureMenu, SIGNAL(triggered(QAction *)), this, SLOT(selectFeature(QAction *)));
	menuBar()->addMenu(featureMenu);

	colorMenu = new QMenu(tr("Select Color Code"));
	connect(colorMenu, SIGNAL(triggered(QAction *)), this, SLOT(selectColorCode(QAction *)));
	menuBar()->addMenu(colorMenu);

	settingMenu = new QMenu(tr("Settings"));
	menuBar()->addMenu(settingMenu);
	radiusAction = new QAction("Set Sphere Radius", this);
	connect (radiusAction, SIGNAL(triggered()), this, SLOT(setSphereRadius()));
	settingMenu->addAction(radiusAction);

	setupColorCodeMenu();

	QVTK = new QVTKWidget(this);
	Renderer = vtkSmartPointer<vtkRenderer>::New();
	Renderer->SetBackground(0,0,0);
	QVTK->GetRenderWindow()->AddRenderer(Renderer);
	this->setCentralWidget(QVTK);

	setWindowTitle(tr("Render Window"));
}

void FTKRenderWindow::setModels(vtkSmartPointer<vtkTable> tbl, ObjectSelection * sels)
{
	table = tbl;
	if(!sels)
		selection = new ObjectSelection();
	else
		selection = sels;
	updateFeatureMenu();
	//updateRenderView();
}

void FTKRenderWindow::update()
{
	//!detects changes in the linked space 
	this->updateFeatureMenu();
	QWidget::update();
}

void FTKRenderWindow::closeEvent(QCloseEvent *event)
{
	//! Detects window closing and emits signal 
	emit closing(this);
	event->accept();
}

void FTKRenderWindow::selectFeature(QAction *action)
{
	//! Sets the feature for the color coding of the objects
	if(action->toolTip().toInt() != selectedFeature)
	{
		selectedFeature = action->toolTip().toInt();
		action->setChecked(true);
		selectedColorCode = 1;
		updateRenderView();
	}
}

void FTKRenderWindow::selectColorCode(QAction *action)
{
	int previous_code = selectedColorCode;
	selectedColorCode = action->toolTip().toInt();
	//discrete color code dialog comes here
	if(selectedColorCode == 1)
	{
		classColorMap.clear();
		action->setChecked(true);
	}
	else if(selectedColorCode == 2)
	{
		std::string current_column = table->GetColumnName(selectedFeature+3);
		if(current_column.find("prediction") == std::string::npos )
		{
			selectedColorCode = previous_code;
			std::cout << "Discrete color coding is only for classification results \n";
			return;
		}

		int max_class = 0;
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			if(table->GetValue(row,selectedFeature+3).ToInt() > max_class)
				max_class = table->GetValue(row, selectedFeature+3).ToInt();
		}

		QVector<QString> colors;
		std::map< std::string, std::vector< double > >::iterator it;
		for(it=discreteColorMap.begin() ; it != discreteColorMap.end(); it++)
		{
			colors.push_back(QString::fromStdString((*it).first));
		}

		DiscreteColorDialog *dialog = new DiscreteColorDialog(max_class, colors, this);
		if( dialog->exec() )
		{
			classColorMap = dialog->GetClassColorMap();
			action->setChecked(true);
		}
		else
		{
			selectedColorCode = previous_code;
			return;
		}
		delete dialog;
	}	

	if(selectedColorCode != previous_code)
		updateRenderView();
}

void FTKRenderWindow::setSphereRadius()
{
	bool ok1;
	double tmp = QInputDialog::getDouble(this, tr("Sphere Radius"),tr("Sphere Radius:"), this->radius, 1, 1000, 1, &ok1);
	if(ok1)
	{
		if( abs(this->radius - tmp) > 1e-6)
		{
			this->radius = tmp;
			updateRenderView();
		}
	}
}

void FTKRenderWindow::updateRenderView(void)
{
	vnl_vector<double> feat_col((int)table->GetNumberOfRows());

	if(selectedColorCode == 1)
	{
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			feat_col[row] = table->GetValue(row, selectedFeature+3).ToDouble();
		}	
		double min = feat_col.min_value();
		double max = feat_col.max_value();
		//std::cout << min << "_" << max << "\n";
		for(int i=0; i<feat_col.size(); ++i)
		{
			feat_col[i] = (((feat_col[i] - min)/(max - min))-0.5)*(2);
		}
		//std::cout << feat_col.min_value() << "_" << feat_col.max_value() << "\n";
	}

	if(selectedColorCode == 2)
	{
		for(int row=0; row<(int)table->GetNumberOfRows(); ++row)
		{
			feat_col[row] = table->GetValue(row, selectedFeature+3).ToDouble();
		}		
	}

	std::set<long int> sels = selection->getSelections();
	std::vector<int> ids(sels.begin(), sels.end());

	Renderer->RemoveAllViewProps();
	for(int row=0; row<(int)table->GetNumberOfRows(); row++)
	{
		//std::cout << "myCentroid : " << row <<"\r";
		if(ids.size() == 0);
		else
		{
			std::vector<int>::iterator posn1 = std::find(ids.begin(), ids.end(), table->GetValue(row, 0).ToInt());
			if(posn1 == ids.end())
				return;
		}

		vtkSmartPointer<vtkSphereSource> centroidSphere =  vtkSmartPointer<vtkSphereSource>::New();
		centroidSphere->SetRadius(this->radius);	
		vtkSmartPointer<vtkPolyDataMapper> centroidMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
		centroidMapper->SetInputConnection(centroidSphere->GetOutputPort()); 
		vtkSmartPointer<vtkActor> centroidactor = vtkSmartPointer<vtkActor>::New();
		centroidactor->SetMapper(centroidMapper);
		centroidactor->GetProperty()->SetOpacity(1);
		r_g_b rgbscalar = GetRGBValue(feat_col[row]);
		centroidactor->GetProperty()->SetColor(rgbscalar.r, rgbscalar.g, rgbscalar.b);
		centroidactor->SetPosition(table->GetValue(row, 1).ToDouble(),table->GetValue(row, 2).ToDouble(),table->GetValue(row, 3).ToDouble());
		Renderer->AddActor(centroidactor);
	}
	 
	QVTK->GetRenderWindow()->Render();
	//Interactor->Start();
	
}

r_g_b FTKRenderWindow::GetRGBValue(double val)
{
	if(selectedColorCode == 1)
	{
		int index = 64 * (val+1) - 1;   // when val = 1; index should be the max index
		if( index >= NUM_COLORS)
		{
			index = NUM_COLORS - 1;
		}
		else if( index < 0)
		{
			index = 0;
		}
		return REN_COLOR_MAP[index];
	}

	else
	{
		std::vector<double> disc_color = discreteColorMap[classColorMap[(int)val]];
		return r_g_b(disc_color[0], disc_color[1], disc_color[2]);
	}
}

void FTKRenderWindow::updateFeatureMenu(void)
{
	if(!table) return;

	featureMenu->clear();

	QActionGroup *featureGroup = new QActionGroup(this);
	
	for (int c=4; c<table->GetNumberOfColumns(); ++c)
	{
		QString name = QString(table->GetColumnName(c));
		QAction *featureAct = new QAction( name, this );
		featureAct->setToolTip( QString::number(c-3) );
		featureAct->setCheckable(true);

		featureMenu->addAction(featureAct);
		featureGroup->addAction(featureAct);

		if((c-3) == selectedFeature)
			featureAct->setChecked(true);

	}
	updateRenderView();
}

void FTKRenderWindow::setupColorCodeMenu(void)
{
	colorMenu->clear();

	QActionGroup *colorGroup = new QActionGroup(this);
	
	QAction *continuousColorAct = new QAction( QString("Continuous (Color changes with feature value)"), this );
	continuousColorAct->setToolTip( QString::number(1) );
	continuousColorAct->setCheckable(true);
	colorMenu->addAction(continuousColorAct);
	colorGroup->addAction(continuousColorAct);

	QAction *discreteColorAct = new QAction( QString("Discrete (Choose a color for each class)"), this );
	discreteColorAct->setToolTip( QString::number(2) );
	discreteColorAct->setCheckable(true);
	colorMenu->addAction(discreteColorAct);
	colorGroup->addAction(discreteColorAct);

	if(selectedColorCode == 1)
		continuousColorAct->setChecked(true);
	else
		discreteColorAct->setChecked(true);

}

void FTKRenderWindow::createDiscreteColorMap(void)
{
	std::vector< double > color;
	color.push_back(1.0); color.push_back(0.0); color.push_back(0.0);
	discreteColorMap["Red"] = color;
	color[0] = 0.0; color[1] = 1.0; color[2] = 0.0;
	discreteColorMap["Green"] = color;
	color[0] = 0.0; color[1] = 0.0; color[2] = 1.0;
	discreteColorMap["Blue"] = color;
	color[0] = 0.0; color[1] = 1.0; color[2] = 1.0;
	discreteColorMap["Cyan"] = color;
	color[0] = 1.0; color[1] = 0.6471; color[2] = 0.0;
	discreteColorMap["Orange"] = color;
	color[0] = 0.9333; color[1] = 0.5098; color[2] = 0.9333;
	discreteColorMap["Violet"] = color;
	color[0] = 1.0; color[1] = 1.0; color[2] = 0.0;
	discreteColorMap["Yellow"] = color;
	color[0] = 0.0; color[1] = 0.3922; color[2] = 0.0;
	discreteColorMap["Dark Green"] = color;
	color[0] = 0.2549; color[1] = 0.4118; color[2] = 0.8039;
	discreteColorMap["Royal Blue"] = color;
	color[0] = 0.7451; color[1] = 0.7451; color[2] = 0.7451;
	discreteColorMap["Gray"] = color;
}


//void FTKRenderWindow::BuildNetwork(void)
//{
//	std::cout << "Started building network from feature table\n";
//
//	data_network = new DataNetwork();
//	for(int r=0; r<(int)table->GetNumberOfRows(); ++r)
//	{
//		DataNode d_node;
//		d_node.id = table->GetValue(r,0).ToInt();
//		d_node.x = table->GetValue(r,1).ToInt();
//		d_node.y = table->GetValue(r,2).ToInt();
//		d_node.z = table->GetValue(r,3).ToInt();
//		data_network->GetDataNodesPointer().push_back(d_node);
//	}	
//
//	std::cout << "Finished Building Network\n";
//}


//############################################################################################
//############################################################################################
// A DIALOG TO GET COLOR MAPPING FOR CLASSIFICATION RESULTS
//############################################################################################

DiscreteColorDialog::DiscreteColorDialog(int max_class, QVector<QString> cols, QWidget *parent)
: QDialog(parent)
{	
	colors = cols;

	layout = new QVBoxLayout;
	for(int i=0; i<max_class; ++i)
	{
		std::stringstream ss;
		ss << i+1;
		QLabel * classLabel = new QLabel(QString::fromStdString("Class " + ss.str() + ":"));
		
		QComboBox * classCombo = new QComboBox();
		for(int v = 0; v<colors.size(); ++v)
		{
			classCombo->addItem(colors[v]);
		}
		classCombos.push_back(classCombo);

		QHBoxLayout * classLayout = new QHBoxLayout;
		classLayout->addWidget(classLabel);
		classLayout->addWidget(classCombo);

		layout->addLayout(classLayout);
	}

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	bLayout = new QHBoxLayout;
	bLayout->addStretch(20);
	bLayout->addWidget(okButton);
	layout->addLayout(bLayout);

	this->setLayout(layout);

	this->setWindowTitle(tr("Choose Class Colors"));
	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
}


std::map< int, std::string > DiscreteColorDialog::GetClassColorMap(void)
{
	std::map< int, std::string > colorMap;
	for(int c=0; c<classCombos.size(); ++c)
	{
		colorMap[c+1] = colors[classCombos[c]->currentIndex()].toStdString();
	}

	return colorMap;
}


//############################################################################################
//############################################################################################
// DATA STRUCTURES USED TO RENDER THE NETWORK
//############################################################################################

//DataNetwork::DataNetwork()
//{
//}
//
//DataNetwork::~DataNetwork()
//{
//	//int counter = 0;
//	//for(int n=0; n<data_nodes.size(); n++)
//	//	delete data_nodes[counter];
//	//counter = 0;
//	//for(int l=0; l<data_links.size(); l++)
//	//	delete data_links[counter];	
//}