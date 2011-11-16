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

#include "HistoWindow.h"
#include <vtkAnnotationLink.h>
#include <vtkDataRepresentation.h>
#include <vtkCommand.h>
#include <vtkSelectionNode.h>
#include <vtkIdTypeArray.h>
#include <vtkSelection.h>
#include <vtkQtChartMouseSelectionHandler.h>
#include <vtkQtChartArea.h>
#include <vtkQtChartArea.h>
//Constructor
HistoWindow::HistoWindow(QWidget *parent)
: QMainWindow(parent)
{
  columnNum=0;
  columnName = "Feature";
  numofbins=5;
  status=true; //we assume there will be more than 2 bins (false when less then 2 bins!!)
  distanceToUpperBound = 0.0001; //Change this for more precision
  normalized=false; //data is not normalized initially

  selection = NULL;
  m_table = NULL;
  hisTable = NULL;

  //Setup menu:
  optionsMenu = menuBar()->addMenu(tr("&Options"));
  columnMenu = new QMenu(tr("Set Column"));
  connect(columnMenu, SIGNAL(triggered(QAction *)), this, SLOT(columnChange(QAction *)));
  optionsMenu->addMenu(columnMenu);
  binsMenu = new QMenu(tr("Set Number of Bins"));
  connect(binsMenu, SIGNAL(triggered(QAction *)), this, SLOT(binsChange(QAction *)));
  optionsMenu->addMenu(binsMenu);

  //Setup ChartView:
  chartView = vtkSmartPointer<vtkQtBarChartView>::New();
  chartView->SetupDefaultInteractor();
  //Set the title of the histogram
  chartView->SetTitle("  ");
  chartView->SetAxisTitle(0,"Frequency");
  chartView->SetAxisTitle(1, columnName.c_str() );

#if(CHART_IS_WIDGET)
		setCentralWidget( chartView->GetWidget() );
#else
		QLabel *label = new QLabel(tr("Use this window to select options for the histogram"));
		label->setWordWrap(true);
		label->setAlignment(Qt::AlignCenter);
		setCentralWidget( label );
#endif

  setWindowTitle(tr("Histogram"));
  this->resize(700,500);
}

void HistoWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
} 

void HistoWindow::setModels(vtkSmartPointer<vtkTable> table, ObjectSelection * sels)
{
	m_table = table;
	selection = sels;

	this->update();
	this->selectionCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	this->selectionCallback->SetClientData(this);
    this->selectionCallback->SetCallback ( SelectionCallbackFunction);
	vtkAnnotationLink *link = chartView->GetRepresentation()->GetAnnotationLink();
	link->AddObserver(vtkCommand::AnnotationChangedEvent, this->selectionCallback);
	//chartView->SetupDefaultInteractor();
}

void HistoWindow::update(void)
{
	this->updateOptionMenus();
	this->SyncModel();
	//this->Normalize(); //What the heck why would you normalize a histogram
	this->SetNumofBins(numofbins);
	this->findFrequencies();    
	this->setBucketNames(); 
	this->ConstructBarChart();
}

void HistoWindow::updateOptionMenus()
{
	if(!m_table) return;

	//Add a new Action for each column for each menu item:
	columnMenu->clear();

	QActionGroup *columnGroup = new QActionGroup(this);
	for (int c=0; c<m_table->GetNumberOfColumns(); ++c)
	{
		QString name = QString( m_table->GetColumnName(c) );
		QAction *xAct = new QAction( name, this );
		xAct->setToolTip( QString::number(c) );
		xAct->setCheckable(true);
		columnMenu->addAction(xAct);
		columnGroup->addAction(xAct);

		if(c == columnNum)
		{
			xAct->setChecked(true);
		}
	}

	//Setup bins menu:
	binsMenu->clear();

	QActionGroup *binsGroup = new QActionGroup(this);
	for(int i=2; i<=100; ++i)
	{
		QAction *xAct = new QAction( QString::number(i), this );
		xAct->setCheckable(true);
		binsMenu->addAction(xAct);
		binsGroup->addAction(xAct);

		if(i == numofbins)
		{
		  xAct->setChecked(true);
		}
	}
}

void HistoWindow::columnChange(QAction *action)
{
  columnNum = action->toolTip().toInt();
  columnName = action->text().toStdString();
  chartView->SetAxisTitle(1, columnName.c_str() );
  this->SyncModel();
  //this->Normalize();
  this->SetNumofBins(numofbins);
  this->findFrequencies();    
  this->setBucketNames(); 
  this->ConstructBarChart();
  action->setChecked(true);
}

void HistoWindow::binsChange(QAction *action)
{
  numofbins = action->text().toInt();
  this->SetNumofBins(numofbins);
  this->findFrequencies();    
  this->setBucketNames(); 
  this->ConstructBarChart();
}

void HistoWindow::SyncModel()
{
	if(!m_table) return;
	if(columnNum >= m_table->GetNumberOfColumns())
		return;

	data.clear();
  
	for(int r=0; r<m_table->GetNumberOfRows(); ++r)
	{
		double val = m_table->GetValue(r,columnNum).ToDouble();
		data.insert(val);
	}
}

void HistoWindow::Normalize() 
{ 
  std::multiset<double> dataTmp;
  double sum=0,mean=0,sumVariance=0,variance=0,stdDeviation=0;
  //Find the sum first
  std::multiset<double>::iterator pos;
  for(pos=data.begin(); pos!=data.end(); ++pos)
    sum +=*pos;

  if (data.size() > 0) {
    //Compute the mean
    mean = sum/data.size();
    //Compute the variance
    for(pos=data.begin(); pos!=data.end(); ++pos)
      sumVariance+= ((*pos - mean) * (*pos - mean));
    variance = sumVariance / data.size();
    //Compute the standard Deviation
    stdDeviation = sqrt(variance);
    // We have everything to normalize the data
    if (stdDeviation != 0) {
      //Now normalize the data.
      for(pos=data.begin(); pos!=data.end(); ++pos)
        dataTmp.insert((*pos - mean)/stdDeviation);
      data = dataTmp;
      //Set the signal
      normalized = true;
    } else {
      cerr<<"The standard deviation is 0. Cannot Normalize!"<<endl;
      normalized = false;
    }
  } else {
    cerr<<"No data...Cannot Normalize!"<<endl;
    normalized = false;
  }  
}

void HistoWindow::SetNumofBins(int n)
{
  if((n>1) && (n<1001))
    numofbins=n;
  else 
  {
    numofbins=0;
    status=false;
  }
}

// Store the frequencies in v
bool HistoWindow::findFrequencies() 
{
  int i;
  //initialize the buckets first
  initMap(result_fq,numofbins);
  std::map<int, int>::iterator number;
  std::multiset<double>::iterator pos;
  double bucketSize;

  if (status) 
  {    
    //range of data, size of buckets, and how much each bucket needs to shift
    double range= *(max_element(data.begin(), data.end())) - *(min_element(data.begin(), data.end()));
    bucketSize = range / numofbins;    	
	double left_offset = *(min_element(data.begin(), data.end()));

    double n,n1,n2;
	int count=0;
    for (pos=data.begin(); pos!=data.end();pos++) {
	  count++;
      //find the bucket for the number
      n=*pos;     
      for (int k=0;k<numofbins;k++) {
        n1 = k*bucketSize + left_offset;
        n2 = (k+1)*bucketSize + left_offset;
        //cout<<"n1 : "<<n1<<" n2: "<<n2<<endl;
        if (k == (numofbins - 1)) {
			if ((n >= n1)  && (n <= n2)) {
            result_fq[k]++;
			}
        } else
          //n1 and n2 have to be double but vc++ cannot compare doubles and flips true and false
          // for some reason. Check the first if statement below
		  if ((n >= n1)  && (n < n2)){
              result_fq[k]++;       
		  }
			  
      }
    }
	cout << count<<" elements binned"<<endl;

    for (i=0;i<numofbins;i++ ) 
      cout<<"i : "<<i<<" Frequency "<<result_fq[i]<<endl;
    return true;
  } else //num of bins is less than 2!
    return false;
}

void HistoWindow::setBucketNames() 
{
  std::string bucketName;
  if (status) { 
    //range of data, size of buckets, and how much each bucket needs to shift
    double range= *(max_element(data.begin(), data.end())) - *(min_element(data.begin(), data.end()));
    double bucketSize = range / numofbins;    	
	double left_offset = *(min_element(data.begin(), data.end()));

    for (int i=0; i<numofbins; i++) {   
      std::stringstream p1,p2;      
      if (i == (numofbins - 1)){
        p1 << (bucketSize * i) + left_offset;
        p2 << *(max_element(data.begin(), data.end()));
	  }else{
        p1 << (bucketSize * i) + left_offset;     
        p2 << (bucketSize * (i+1)) + left_offset;
      }

      bucketName = p1.str() + " - " + p2.str();
	  cout<<"Bin "<<i<<": "<<bucketName<<endl;
      names.push_back(bucketName);
    }
  }
}

void HistoWindow::ConstructBarChart() 
{
	if (status) 
	{
		vtkSmartPointer<vtkDoubleArray> column1 = vtkSmartPointer<vtkDoubleArray>::New();
		//vtkSmartPointer<vtkDoubleArray> column2 = vtkSmartPointer<vtkDoubleArray>::New();
		column1->SetName("C1");
		//column2->SetName("C2");

		for (unsigned int i=0; i<(unsigned int)numofbins; ++i) 
		{
			column1->InsertNextValue(result_fq[i]);
			//column2->InsertNextValue(col2[i]);
		}

		// Create a table
		hisTable = vtkSmartPointer<vtkTable>::New();
		// Add the data to the table
		hisTable->AddColumn(column1);
		//table->AddColumn(column2);

		// Set the chart title
		//chartView->SetTitle("HI EVERYONE");
		//chartView->SetAxisTitle(0,"Frequency");
		//chartView->SetAxisTitle(1,"Feature");
		//chartView->SetAxisTitle(1,feature);
		chartView->RemoveAllRepresentations();
		chartView->AddRepresentationFromInput(hisTable);
		chartView->Update();

		#if(!CHART_IS_WIDGET)
			chartView->Show();
		#endif
	}
}


void HistoWindow::initMap(std::map<int,int> &v, int n) {
  for (int i=0;i<n; i++) 
    v[i]=0;
}

/*
bool comp(double n, double n1,double n2) {
  if ((n >= n1)  && (n < n2))
    return true;
  else return false;
};
*/

void HistoWindow::SelectionCallbackFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData )
{
	vtkAnnotationLink* annotationLink = static_cast<vtkAnnotationLink*>(caller);
	vtkSelection* selection = annotationLink->GetCurrentSelection();
	HistoWindow* histWin = (HistoWindow*)clientData;

	vtkSelectionNode* vertices = NULL;


	if( selection->GetNode(0))
	{
		if( selection->GetNode(0)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(0);
		}
	}

	if( selection->GetNode(1))
	{
		if( selection->GetNode(1)->GetFieldType() == vtkSelectionNode::VERTEX)
		{
			vertices = selection->GetNode(1);
		}
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
		//histWin->SetSelectedIds( IDs);
		//histWin->SetSelectedIds2();  // only select the selected feature columns
	}
	histWin->chartView->Update();
}