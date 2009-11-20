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


#include "histogram.h"

#ifdef true
#  error "true is defined as a macro"
#endif
#ifdef false
#  error "false is defined as a macro"
#endif

bool Histogram::ReadHistogramData(char* fileName)
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
 
  while(levelOneElement)
    {
    nodeName = levelOneElement->Value();

	//Check if this is a histogram_data file
    if (strcmp(nodeName,"histogram_data") == 0)
      {

		  /*
		  TiXmlAttribute* node = attributeSet.Find( name );
	if ( node )
	{
		attributeSet.Remove( node );
		delete node;
	}
	*/

	  //Construct the bar title from the attribute values
	  //strcat(barTitle,levelOneElement->Attribute("Class_Membership"));
	  strcat(barTitle," -  Frequency (Y-axis)  vs  ");
	  strcat(barTitle,levelOneElement->Attribute("feature"));
	  strcat(barTitle," (X-axis)     Cell Type: ");
	  strcat(barTitle,levelOneElement->Attribute("Class_Membership"));

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
/*
  int k=1;
  multiset<double>::iterator i;
  for (i=data.begin();i!=data.end();i++ ) { 
	  cout<<"k - "<<k<<" datai : "<<*i<<endl;
		  k++;
  };
*/
  return true;
}

void Histogram::Normalize() {
	multiset<double> dataTmp;
	//double minNumber= *min_element(data.begin(), data.end());
	//double maxNumber= *max_element(data.begin(), data.end());
	//double diff = maxNumber - minNumber;

	double sum=0;
	//Find the sum first
	for(pos=data.begin(); pos!=data.end(); pos++)
		sum +=*pos;

	//Now normalize the data.
	for(pos=data.begin(); pos!=data.end(); pos++)
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
	
};



/*
double Histogram::SetDiffToUpperBound(double n) {
	double stepSize;
	if ((10000 * n) <= 1) status=false; //change this if more than 3 decimals are required
	else if ((1000 * n) <= 1) stepSize=0.00001;
	else if ((100 * n) <= 1) stepSize=0.0001;
	else if ((10 * n) <= 1) stepSize=0.001;
	else if ((1 * n) <= 1) stepSize=0.01;

		((100 * n) <= 1)  
		((10 * n) <= 1))

			stepSize = n/10;
	else if ((1<n) && (n<2)) stepSize=1;
	return stepSize;
}
*/

void Histogram::setBucketNames() {
	string bucketName;
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
			stringstream p1,p2;			
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

void Histogram::initMap(map<int,int> &v, int n) {
	for (int i=0;i<n; i++) 
		v[i]=0;
}

bool comp(double n, double n1,double n2) {
	if ((n >= n1)  && (n < n2))
		return true;
	else return false;
};


// Store the frequencies in v
bool Histogram::findFrequencies() {
	int i;
	//initialize the buckets first
	initMap(result_fq,numofbins);
	map<int, int>::iterator number;
	//multiset<double>::iterator pos;
	double bucketSize;

	if (status) {
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

};


//Records the histogram
void Histogram::RecordImage (char* imgFileName) {

	vtkWindowToImageFilter* img = vtkWindowToImageFilter::New();
	img->SetInput(renWin);
	vtkTIFFWriter *writer = vtkTIFFWriter::New();
	writer->SetInput(img->GetOutput());
	char buff[1024];
	strcat(imgFileName,".tiff"); //Image will be saved as a tiff file
	strcpy(buff,imgFileName);
	writer->SetFileName(buff);
	writer->Update();	
}



	//Define the constructor
Histogram::Histogram() {
		frequencies = vtkIntArray::New();
		dobj = vtkDataObject::New();
		actor = vtkBarChartActor::New();
		//barTitle=new char();
		numofbins=0;
		status=true; //we assume there will be more than 2 bins
		distanceToUpperBound = 0.0001; //Change this for more precision
		normalized=false; //data is nor normalized initially
}

	//Define the destructor
Histogram::~Histogram() {
		frequencies->Delete();
		dobj->Delete();
		actor->Delete();
		ren1->Delete();
		renWin->Delete();
		iren->Delete();
		
	};


void Histogram::ConstructBarChart() {

	if (status) {	
		frequencies->SetNumberOfTuples(numofbins);
 
		for (int i=0; i<numofbins; i++)
			frequencies->SetTuple1(i,result_fq[i]);
    
		dobj->GetFieldData()->AddArray(frequencies);
 
		actor->SetInput(dobj);
		actor->GetProperty()->SetColor(1,1,0);
		
		//vtkTextProperty* p=vtkTextProperty::New();
		//p->SetFontSize(200);
		//p->S
		//actor->SetTitleTextProperty(p);

		actor->SetTitle(barTitle);
		actor->SetYTitle("Frequency");		
		actor->GetPositionCoordinate()->SetValue(0.15,0.05,0.0);
		actor->GetPosition2Coordinate()->SetValue(0.95,0.85,0.0);
		actor->GetProperty()->SetColor(1,1,0);
		actor->GetLegendActor()->SetNumberOfEntries(numofbins);
		for (int i=0; i<numofbins; i++)	{
			double red=vtkMath::Random(0,1);
			double green=vtkMath::Random(0,1);
			double blue=vtkMath::Random(0,1);
			actor->SetBarColor(i,red,green,blue);
		}

		for (int i=0; i<numofbins; i++) 
			actor->SetBarLabel(i, names[i].c_str());
  
		actor->LegendVisibilityOn();

		// Set text colors (same as actor for backward compat with test)
		actor->GetTitleTextProperty()->SetColor(1,1,0);
		actor->GetLabelTextProperty()->SetColor(1,1,0);
	}
	
}

void Histogram::RenderWin(){

	if (status) {
		// Create the RenderWindow, Renderer and both Actors
		
		ren1 = vtkRenderer::New();
		renWin = vtkRenderWindow::New();
		renWin->AddRenderer(ren1);
		iren = vtkRenderWindowInteractor::New();
		iren->SetRenderWindow(renWin);

		ren1->AddActor(actor);
		ren1->SetBackground(0,0,0);
		renWin->SetSize(1000,400);

	
		// render the image
		//ren1->ResetCameraClippingRange();
		renWin->Render();
		RecordImage(imgFile);	
		
		iren->Initialize();
		iren->Start();
	} else {cerr<<"The number of bins has to be between 2 and 10"<<endl;
			return;
	}

}

//----------------------------------------------------------------------------
int main( int argc, char * argv [] )
{
	if (argc <4) {
			cout<<"Incorrect Usage! XML file name, number of bins and a" <<endl;
			cout<<"filename for the histogram to be saved are required!" <<endl;
			return 0;
	};
	
	Histogram* his = new Histogram();
	his->imgFile=argv[3]; //File name for the histogram to be saved
	
	if (! his->ReadHistogramData(argv[1])) {
		cerr<<"No data!"<<endl;
		return 0;
	} else {		
		his->Normalize();
		his->SetNumofBins(atoi(argv[2]));
		//his->SetNumofBins(3);
		//Find the frequencies now
		his->findFrequencies();		
		his->setBucketNames();
		his->ConstructBarChart();
		his->RenderWin();
		//his->RecordImage(argv[3]);
	}	
  return 0;  
}
