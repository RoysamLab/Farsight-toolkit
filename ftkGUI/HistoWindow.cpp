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

//#include "Common.h"
#include "HistoWindow.h"
//Constructor
// Requires a selection model (also requires itemmodel)
HistoWindow::HistoWindow(QItemSelectionModel *mod, QWidget *parent)
: QMainWindow(parent)
{
  model = (QAbstractItemModel *)mod->model();

  columnNum=0;
  columnName = "Feature";
  numofbins=5;
  status=true; //we assume there will be more than 2 bins (false when less then 2 bins!!)
  distanceToUpperBound = 0.0001; //Change this for more precision
  normalized=false; //data is nor normalized initially

  //Setup menu:
  optionsMenu = menuBar()->addMenu(tr("&Options"));
  columnMenu = new QMenu(tr("Set Column"));
  connect(columnMenu, SIGNAL(triggered(QAction *)), this, SLOT(columnChange(QAction *)));
  optionsMenu->addMenu(columnMenu);
  binsMenu = new QMenu(tr("Set Number of Bins"));
  connect(binsMenu, SIGNAL(triggered(QAction *)), this, SLOT(binsChange(QAction *)));
  optionsMenu->addMenu(binsMenu);
  this->updateOptionMenus();

  //Setup ChartView:
  chartView = vtkSmartPointer<vtkQtBarChartView>::New();
  chartView->SetupDefaultInteractor();
  //Set the title of the histogram
  chartView->SetTitle("  ");
  chartView->SetAxisTitle(0,"Frequency");
  chartView->SetAxisTitle(1, columnName.c_str() );
 
  //check if farsight is built with VTK 5.4
  #if VTK_MINOR_VERSION==4
    #define VTK_5_4
  #endif

  #ifdef VTK_5_4
    QLabel *label = new QLabel(tr("Use this window to select options for the histogram"));
    label->setWordWrap(true);
    label->setAlignment(Qt::AlignCenter);
    setCentralWidget( label );
  #else
    setCentralWidget( chartView->GetWidget() );
  #endif

  setWindowTitle(tr("Histogram"));
  this->resize(700,500);

  this->SyncModel();
  this->Normalize();
  this->SetNumofBins(numofbins);
  this->findFrequencies();    
  this->setBucketNames(); 
  this->ConstructBarChart();

  connect(model, SIGNAL(dataChanged(const QModelIndex &, const QModelIndex &)), this, SLOT(modelChange(const QModelIndex &, const QModelIndex &)));
}

void HistoWindow::updateOptionMenus()
{
  //Add a new Action for each column for each menu item:
  columnMenu->clear();

  QActionGroup *columnGroup = new QActionGroup(this);
  for (int c=0; c<model->columnCount(); ++c)
  {
    QString name = model->headerData(c,Qt::Horizontal).toString();
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
  for(int i=2; i<=10; ++i)
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

void HistoWindow::modelChange(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
  this->updateOptionMenus();
  this->SyncModel();
  this->Normalize();
  this->SetNumofBins(numofbins);
  this->findFrequencies();    
  this->setBucketNames(); 
  this->ConstructBarChart();
}

void HistoWindow::SyncModel()
{
  if(columnNum >= (int)model->columnCount())
    return;

  data.clear();
  
  for(int r=0; r<(int)model->rowCount(); ++r)
  {
    double val = model->data( model->index(r,columnNum)).toDouble();
    data.insert(val);
  }
}

/*
bool HistoWindow::ReadHistogramData(const char* fileName)
{  
  TiXmlDocument doc(fileName);
  doc.LoadFile();
  TiXmlHandle docHandle( &doc );
  TiXmlElement* levelOneElement =
  docHandle.FirstChild("histogram_data").Element();
  //docHandle.FirstChild("Trace").FirstChild().Element();
  TiXmlElement *levelTwoElement;
  const char *nodeName;
  double num;

  data.clear();
 
  while(levelOneElement)
    {
    nodeName = levelOneElement->Value();

  //Check if this is a histogram_data file
    if (strcmp(nodeName,"histogram_data") == 0)
      {

    //Construct the bar title from the attribute values
    //strcat(barTitle,levelOneElement->Attribute("Class_Membership"));
    //strcat(barTitle," -  Frequency (Y-axis)  vs  ");
    //strcat(feature,levelOneElement->Attribute("feature"));
    //strcat(barTitle,levelOneElement->Attribute("feature"));
    //strcat(barTitle," (X-axis)     Cell Type: ");
    //strcat(barTitle,levelOneElement->Attribute("Class_Membership"));

      levelTwoElement = levelOneElement->FirstChildElement();
      while(levelTwoElement)
        {
        nodeName = (char*)levelTwoElement->Value();

        //Check if there is tag <d> first then get the data located between <d> and </d>
    if (strcmp(nodeName,"d") == 0){
      // Data coming from the XML file is char*. Convert it to a numeric value 
      num=atof(levelTwoElement->GetText());
      data.insert(num);     
    }     
        else if (strcmp(nodeName,"text") == 0) 
          {
          //Do Nothing
          }
        else
          {
          cerr << "XML File contains a tag that cannot be identified! " << nodeName << endl;
          return false;
          }
    levelTwoElement = levelTwoElement->NextSiblingElement("d");
    //cout<<"level two "<<levelTwoElement<<endl;
    }
      }
    else
      {
      cerr << "Incorrect Histogram Data format! " << nodeName << endl;
      return false;
      }
  levelOneElement=levelOneElement->NextSiblingElement();
  }

  // Check if we read any data from the XML file.
  // If not, we cannot construct a histogram
  if (data.size() == 0) {status=false;
             return false;} 

  return true;
}
*/

// Description:
// Normalizes the data using this formula:
// (Xi - Xmean)/StdDeviation
// StdDeviation = sqrt(variance)
// variance = (sum (sqr(Xi-Xmean)))/N
// where N is the number of elements

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
  if ((n>1) && (n<11)) 
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
    table = vtkSmartPointer<vtkTable>::New();
    // Add the data to the table
    table->AddColumn(column1);
    //table->AddColumn(column2);

    // Set the chart title
    //chartView->SetTitle("HI EVERYONE");
    //chartView->SetAxisTitle(0,"Frequency");
    //chartView->SetAxisTitle(1,"Feature");
    //chartView->SetAxisTitle(1,feature);
    chartView->RemoveAllRepresentations();
    chartView->AddRepresentationFromInput(table);
    chartView->Update();

    #ifndef VTK_5_4
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
