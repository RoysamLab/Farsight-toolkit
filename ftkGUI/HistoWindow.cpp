#include "HistoWindow.h"
//Constructor
// Requires a selection model (also requires itemmodel)
HistoWindow::HistoWindow(QItemSelectionModel *mod, QWidget *parent)
: QMainWindow(parent)
{
	model = (QAbstractItemModel *)mod->model();

	columnNum=0;
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
	chartView->SetTitle("HI EVERYONE");
	chartView->SetAxisTitle(0,"Frequency");
	chartView->SetAxisTitle(1,"Feature");
	setCentralWidget( chartView->GetWidget() );
	setWindowTitle(tr("Histogram"));
	this->resize(700,500);

	this->SyncModel();
	this->Normalize();
	this->SetNumofBins(numofbins);
	this->findFrequencies();		
	this->setBucketNames();	
	this->ConstructBarChart();
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
	if(columnNum >= (int)model->columnCount())
		return;

	data.clear();
	
	for(int r=0; r<(int)model->rowCount(); ++r)
	{
		double val = model->data( model->index(r,columnNum)).toDouble();
		data.insert(val);
	}
}

void HistoWindow::Normalize() 
{
	std::multiset<double> dataTmp;
	//double minNumber= *min_element(data.begin(), data.end());
	//double maxNumber= *max_element(data.begin(), data.end());
	//double diff = maxNumber - minNumber;

	double sum=0;
	//Find the sum first
	std::multiset<double>::iterator pos;
	for(pos=data.begin(); pos!=data.end(); ++pos)
		sum +=*pos;

	//Now normalize the data.
	for(pos=data.begin(); pos!=data.end(); ++pos)
    {
		//inefficient hack to workaround the fact that multiset iterators are
    //immutable on gcc
    dataTmp.insert(*pos/sum);
    //*pos = *pos/sum;
    }
  data = dataTmp;

	//Set the signal
	normalized = true;

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



/*
			for (int k=0;k<numofbins;k++) {
				//cout<<"ALLLm KKKK"<<k<<endl;
				if (k == (numofbins - 1)) 
					result_fq[k]++; 
				else {

					//n1 and n2 have to be double but vc++ cannot compare doubles and flips true and false
					// for some reason. Check the first if statement below
					float n1 = k*bucketSize;
					float n2 = (k+1)*bucketSize;
					double n=*pos;
					if ((n >= n1)  && (n < n2)) {
						//cout<<"result "<<k<<" : "<<result_fq[k]<<endl;
						result_fq[k]++;
					}
				}
			}
			
*/


/*
			loc = *pos / bucketSize;
			if (loc == numofbins) result_fq[loc-1]++;   
			else {
				result_fq[loc]=result_fq[loc]+1;
				cout<<"Number is "<<*pos<<" loc: "<<loc<<endl;
			}
*/
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