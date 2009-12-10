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
}

void HistoWindow::update(void)
{
	this->updateOptionMenus();
	this->SyncModel();
	this->Normalize();
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
  this->Normalize();
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

  // Here are the normalized numbers
  /*
  int k=1;
  for(pos=data.begin(); pos!=data.end(); pos++) {
    cout<<"k: "<<k<<" "<<*pos<<endl;
    k++;
  };
  */
  
}

void HistoWindow::SetNumofBins(int n)
{
  //if ((n>1) && (n<11))
	if((n>1) && (n<1000))
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
    pos = max_element(data.begin(), data.end());
    if (*pos != 0) {
      //if (normalized) bucketSize=100/numofbins;
      bucketSize = *pos / numofbins;
      //If the size is big, ignore the decimal
      //We need to work on other cases
      /*
      if (bucketSize >= 2) bucketSize = bucketSize;
      else if (bucketSize >= 0.2) bucketSize = bucketSize;
      else if (bucketSize >= 0.02) bucketSize = bucketSize;
      else bucketSize = bucketSize;
    */
    }
    else {
      cerr<<"The maximum element is 0"<<endl;
      return false;
    }

    double n,n1,n2;
    for (pos=data.begin(); pos!=data.end();pos++) {
      //find the bucket for the number
      n=*pos;     
      for (int k=0;k<numofbins;k++) {
        n1 = k*bucketSize;
        n2 = (k+1)*bucketSize;
        //cout<<"n1 : "<<n1<<" n2: "<<n2<<endl;
        if (k == (numofbins - 1)) {
          if ((n >= n1)  && (n <= n2)) 
            result_fq[k]++;
        } else
          //n1 and n2 have to be double but vc++ cannot compare doubles and flips true and false
          // for some reason. Check the first if statement below
            if ((n >= n1)  && (n < n2))
              result_fq[k]++;       
      }
    }

    for (i=0;i<numofbins;i++ ) 
      cout<<"i : "<<i<<" Frequency "<<result_fq[i]<<endl;
    return true;
  } else //num of bins is less than 2!
    return false;
}

void HistoWindow::setBucketNames() 
{
  std::string bucketName;
  if (status) { // 2<= numofbins <=10
    double maxNumber = GetMaxNumber();    
    double bucketSize=maxNumber/numofbins;
    //double distanceToUpperBound=0;
    //Upper bound of the bins should be decided as follows:
    /*
    if (bucketSize >= 2) distanceToUpperBound = 1;
    else if (bucketSize >= 0.2) distanceToUpperBound = 0.1;
    else if (bucketSize >= 0.02) distanceToUpperBound = 0.01;
    else distanceToUpperBound = 0.001; //not more than 3 decimals for now
    */

    for (int i=0; i<numofbins; i++) {   
      std::stringstream p1,p2;      
      if (i == (numofbins - 1))
        if (bucketSize >= 2) { //if the bin bounds are large numbers, don't show the decimals
          p1<<(((int)bucketSize) * i);          
          p2<<ceil(maxNumber);          
        } else {
          p1<<(bucketSize * i);
          p2<<maxNumber;
          //p2<<((bucketSize * (i+1)) - distanceToUpperBound);
        }
      else {
        //p2<<((((int)(maxNumber/numofbins)) * (i+1)) - 1);
        if (bucketSize >= 2) { //if the bin bounds are large numbers, don't show the decimals
          p1<<(((int)bucketSize) * i);          
          p2<<((((int)bucketSize) * (i+1)) - distanceToUpperBound);
        } else {
          p1<<(bucketSize * i);     
          p2<<((bucketSize * (i+1)) - distanceToUpperBound);
        }

      }

      bucketName = p1.str() + " - " + p2.str();
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
