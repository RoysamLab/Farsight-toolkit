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

#include "ResultLoadWidget.h"

//********************************************************************************
//  ResultLoadWidget
//
//  This widget is a convenience class that loads a result directly from an XML
//  file based on the XML tags.  
//
//  Currently the tags are hard-coded into the model, but in the future this widget
//  can request their names from the user to make the GUI more dynamic and 
//  adaptive.
//********************************************************************************

//Constructor
ResultLoadWidget::ResultLoadWidget(FTKItemModel* mod,QItemSelectionModel* selMod)
{
	model = mod;
	selectionModel = selMod;

	segWin = NULL;
	tblWin = NULL;
	pltWin.clear();

	dataImg = NULL;
	labelImg = NULL;

	filename  = QFileDialog::getOpenFileName(0,"Choose Result File (XML)",".",
			tr("XML Files (*.xml)\n"));

	if( filename.contains(".xml", Qt::CaseInsensitive) )
    {
		extractPath();
		loadFromXML();
	}
}

void ResultLoadWidget::loadFromXML()
{

		loadXML();
		loadOutliers();

		dataImg = new ftkImage();
		dataImg->load(QString(datafname).toStdString());
		labelImg = new ftkImage();
		labelImg->load(QString(labelfname).toStdString());

		CreateNewSegmentationWindow();
		
		CreateNewTableWindow();

		CreateNewPlotWindow();
}

//******************************************************************************
// Load XML simply loads all of the information in the XML file with filename
// into the model m and headerInfo.
// 
//******************************************************************************
void ResultLoadWidget::loadXML()
{
	if (filename == "")
		return;

	vector< vector< string > > features(0);
	vector< string > FeatureNames(0);
	bool namesFilled = false;

	TiXmlDocument doc( filename.toStdString().c_str() );
	doc.LoadFile();

	TiXmlElement* topElement = doc.FirstChildElement();
	if( strstr( topElement->Value(), "nuclei") )
	{
		//get the filenames
		datafname = QString::fromStdString(topElement->Attribute("DataImage"));
		labelfname = QString::fromStdString(topElement->Attribute("LabelImage"));

		TiXmlElement* nucElement = topElement->FirstChildElement();

		//Should now be into the list of nuclei
		while (nucElement)
		{
			vector< string > temp;
			TiXmlElement* featElement = nucElement->FirstChildElement();
			while(featElement)
			{
				if(namesFilled == false)
				{
					FeatureNames.push_back(featElement->Value());
				}
				temp.push_back(featElement->GetText());
				featElement = featElement->NextSiblingElement();
			}
			namesFilled = true;

			features.push_back(temp);
			nucElement = nucElement->NextSiblingElement();
		}
	}

	model->setColumnCount(features[0].size());
	for (int c=0; c<features[0].size(); ++c)
	{
		model->setHeaderData(c, Qt::Horizontal, QString::fromStdString(FeatureNames[c]));
	}

	model->removeRows(0, model->rowCount(QModelIndex()), QModelIndex());
	model->setRowCount(features.size());
	for (int r=0; r<features.size(); ++r)
	{
		for (int c=0; c<features[0].size(); ++c)
		{
			float fval = QString::fromStdString(features[r][c]).toFloat();
			model->setData(model->index(r, c, QModelIndex()), fval);
		}
	}

	model->setColForID(0);
	model->setColForSeedX(1);
	model->setColForSeedY(2);
	model->setColForSeedZ(3);

}

void ResultLoadWidget::extractPath(void)
{
	//Extract path to this XML file
	vector<string> dir;
	dir.clear();

	string ff = filename.toStdString();

	char * pch;
	pch = strtok ( (char*)ff.c_str() , "/" );
	while (pch != NULL)
	{
		dir.push_back(pch);
		pch = strtok (NULL, "/");
	}
	
	//Now we have the directory stucture we need.
	string tpath;
	tpath.clear();
	for (int i=0; i<dir.size()-1; ++i)
	{
		tpath.append(dir[i]);
		tpath.append("/");
	}

	path = QString::fromStdString(tpath);
}

//This function attempts to load outliers from a file called outliers.txt
//It just looks in the directory for it!
void ResultLoadWidget::loadOutliers(void)
{
	string outlierfilename = path.toStdString();
	outlierfilename.append("outliers.txt");

	ifstream inFile; 
	inFile.open(outlierfilename.c_str() );
	if ( !inFile.is_open() )
	{
		std::cerr << "Failed to Load Document: " << outlierfilename << std::endl;
		return;
	}

	const int MAXLINESIZE = 30;	//Numbers could be in scientific notation in this file
	char line[MAXLINESIZE];
	vector< int > outliers(0);

	inFile.getline(line, MAXLINESIZE);
	while ( !inFile.eof() )
	{
		int v = (int)atof(line);
		outliers.push_back(v);

		//Get the next line
		inFile.getline(line, MAXLINESIZE);
	}

	inFile.close();

	//Add to model
	int cols = model->columnCount();
	model->setColumnCount(cols+1);
	model->setHeaderData(cols, Qt::Horizontal, QString::fromStdString("outlier?"));
	
	//Now go through each item and check it in outlier list.
	for(int r=0; r<model->rowCount(); ++r)
	{
		int c = 0;	//column for ID
		int id = ( model->data( model->index(r, c, QModelIndex()) ) ).toInt();
		for (int i=0; i< outliers.size(); ++i)
		{
			if ( id == outliers[i] )
			{
				model->setData(model->index(r, cols, QModelIndex()), 1);
				break;
			}
			else
			{
				model->setData(model->index(r, cols, QModelIndex()), 0);
			}
		}
	}
	model->setColForColoring(cols);
}

//******************************************************************************
// Create a new Plot window and give it the provided model and selection model
//******************************************************************************
void ResultLoadWidget::CreateNewPlotWindow(void)
{
	pltWin.push_back(new PlotWindow());
	connect(pltWin.back(), SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	pltWin.back()->SetDataModel(model);
	pltWin.back()->SetDataSelectionModel(selectionModel);
	pltWin.back()->show();
}

//******************************************************************************
// Create a new table window
//******************************************************************************
void ResultLoadWidget::CreateNewTableWindow(void)
{
	tblWin = new TableWindow();
	connect(tblWin, SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	tblWin->SetDataModel(model);
	tblWin->SetDataSelectionModel(selectionModel);
	tblWin->ResizeToOptimalSize();
	tblWin->show();
}

//******************************************************************************
// Create a new segmentation window.  Shows Images and Segmentation results.
//******************************************************************************
void ResultLoadWidget::CreateNewSegmentationWindow(void)
{
	segWin = new SegmentationWindow();
	connect(segWin, SIGNAL(closing(QWidget*)), this, SLOT(closeWidget(QWidget*)));
	segWin->SetDataModel(model);
	segWin->SetDataSelectionModel(selectionModel);
	segWin->SetChannelImage(dataImg);
	segWin->SetLabelImage(labelImg);
	segWin->show();
}

//*************************************************************************
// THIS SLOT IS CALLED WHEN A WINDOW IS CLOSED
//*************************************************************************
void ResultLoadWidget::closeWidget(QWidget *widget)
{
	std::vector<PlotWindow*>::iterator startIterator = pltWin.begin();
	for (int i=0; i < pltWin.size(); ++i)
	{
		if( widget == pltWin.at(i) )
		{
			pltWin.erase(startIterator+i);
		}
	}
	if(widget == tblWin)
	{
		tblWin = NULL;
	}
	if(widget == segWin)
	{
		segWin = NULL;
	}
}

//******************************************************************************
//Reimplement closeEvent to also close other module windows
//******************************************************************************
void ResultLoadWidget::closeEvent(QCloseEvent *event)
{
	//First Close other widgets
	for (int i=0; i<pltWin.size(); ++i)
	{
		pltWin.at(i)->close();
	}
	pltWin.resize(0);

	if(tblWin)
	{
		tblWin->close();
		tblWin = NULL;
	}
	if(segWin)
	{
		segWin->close();
		segWin = NULL;
	}
	//Then clean up memory
	if(dataImg)
		delete dataImg;
	if(labelImg)
		delete labelImg;
	//Then close myself
	emit closing(this);
	event->accept();
} 



