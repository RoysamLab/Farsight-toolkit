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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>

#include "vtkBarChartActor.h"
#include "vtkFloatArray.h"
#include "vtkDataObject.h"
#include "vtkFieldData.h"
#include "vtkMath.h"
#include "vtkTextProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkProperty2D.h"
#include "vtkLegendBoxActor.h"
#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"
#include "vtkWindowToImageFilter.h"
#include "vtkTIFFWriter.h"
#include <vtkTextProperty.h>

#include "tinyxml.h"
#include <set> //store histogram data in a multiset
#include <algorithm> //for finding maximum_element
#include <map>
#include <sstream> //convert an int to a string
#include <string>

#include "vtkQtBarChartView.h"
#include "vtkQtChartRepresentation.h"

#include "vtkTable.h"
#include "vtkDoubleArray.h"
#include "vtkSmartPointer.h"

using namespace std;

// Aytekin Vargun:
// This class is used for constructing a histogram from a set of numbers in an XML file.
//
// Notes:
// The number of bins should be more than 1 and less than 11. 
// Otherwise the algorithm will not produce the best results.
// If more than 10 bins is required, update the frequencies and 
// setBucketNames methods.
class Histogram : public QWidget
{
    Q_OBJECT;

public:
	//Define the constructor
	Histogram(std::string filename, int numBins, QWidget *parent = 0);

private: //METHODS:
	bool ReadHistogramData(const char* fileName);
	void Normalize();
	void SetNumofBins(int n);
	bool findFrequencies();
	void setBucketNames();
	void ConstructBarChart();

	void initMap(map<int,int> &v, int n);	
	double GetMaxNumber(){ multiset<double>::iterator pos=max_element(data.begin(), data.end());
						   return *pos;};

	//VARIABLES:

	//histogram data will be read from an XMl file to variable data
	multiset<double> data;

	int numofbins;				//
	bool status;				//will be false if num of bins is less than 2 or larger than 10	
	bool normalized;			//false if not normalized
	double distanceToUpperBound;//change this if you need more precision for bar charts and frequencies

	//return the frequencies here.
	//result_fq[0]=25   stores frequency which is 25 for bin-0
	//result_fq[1]=30   stores frequency which is 30 for bin 1 and so on
	map<int,int> result_fq;
	vector<string> names; //Names of bins will be stored here

	//char barTitle[1024];
	//char feature[1024];

	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkQtBarChartView> chartView;

	QVBoxLayout *layout;

	//vtkSmartPointer<vtkIntArray> frequencies;
	//vtkSmartPointer<vtkDataObject> dobj;

	//vtkBarChartActor *actor;

	//These are required for rendering the result
	//vtkRenderer *ren1;
	//vtkRenderWindow *renWin;
	//vtkRenderWindowInteractor *iren;

	//void RenderWin();
	//void RecordImage (char* imgFileName);
	//double SetDiffToUpperBound(double n);

};

#endif
