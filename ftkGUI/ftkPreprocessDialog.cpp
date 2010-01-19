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
#include "ftkPreprocessDialog.h"
#include "ftkPreprocess.h"

//******************************************************************************************
//******************************************************************************************
// A dialog to get the paramaters file for the preprocessing to use and specify the channel 
// if image has more than one:
//******************************************************************************************
ftkPreprocessDialog::ftkPreprocessDialog(QVector<QString> channels, std::string id, ftk::Image::Pointer img, QWidget *parent)
: QDialog(parent)
{
	InitializeFilters();
	filtername = id;
	myImg = img;
	channelLabel = new QLabel("Choose Channel: ");
	channelCombo = new QComboBox();
	
	for(int v = 0; v<channels.size(); ++v)
	{
		channelCombo->addItem(channels.at(v));
	}
	
	QGridLayout *layout = new QGridLayout;
	this->setLayout(layout);
	this->setWindowTitle(tr("Parameters"));

	layout->addWidget(channelLabel,0,0);
	layout->addWidget(channelCombo,0,1);

	switch(FilterValue[id])
	{
	//Mean and Median Filter
	case Filter1:case Filter2:  
		QTParamLabel1 = new QLabel("Window Size - X");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Window Size - Y");
		QTParam2 = new QLineEdit();
		QTParam2->setText(tr("3"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Window Size - Z");
		QTParam3 = new QLineEdit();
		QTParam3->setText(tr("1"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		 break;

	//G and C Anisotropic Diffusion
	case Filter3:case Filter4: 
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit(); 
		QTParam1->setText(tr("0.125"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit();
		QTParam2->setText(tr("5"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Conductance");
		QTParam3 = new QLineEdit();
		QTParam3->setText(tr("2"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		 break;
		 
	//Sigmoid Filters		 
	case Filter5: 
		QTParamLabel1 = new QLabel("Alpha");
		QTParam1 = new QLineEdit(); 
		QTParam1->setText(tr("10"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Beta");
		QTParam2 = new QLineEdit(); 
		QTParam2->setText(tr("150"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Output Minimum");
		QTParam3 = new QLineEdit();
		QTParam3->setText(tr("30"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		QTParamLabel4 = new QLabel("Output Maximum");
		QTParam4 = new QLineEdit(); 
		QTParam4->setText(tr("240"));
		layout->addWidget(QTParamLabel4,4,0);
		layout->addWidget(QTParam4,4,1);

		 break;
		 
	//Grayscale Morphological Filters	 
	case Filter6:case Filter7:case Filter8:case Filter9: 
		QTParamLabel1 = new QLabel("Radius");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		 
		 break;

	//Smoothing Recursive Gaussian	 
	case Filter10:
		
		QTParamLabel1 = new QLabel("Sigma");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number of Threads");
		QTParam2 = new QLineEdit();
		QTParam2->setText(tr("3"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Normalize across Scale ?");
		QTParam5 = new QCheckBox();
		QTParam5->setText(tr("0"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam5,3,1);
		 break;
	
	//	Opening/Closing by Reconstruction 	 	 
	case Filter11:case Filter12:
		
		QTParamLabel1 = new QLabel("Radius");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel3 = new QLabel("Preserve Intensities ?");
		QTParam5 = new QCheckBox();  
		QTParam1->setText(tr("0"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Fully Connected ?");
		QTParam6 = new QCheckBox();
		QTParam1->setText(tr("1"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam6,3,1);
		 break;
	 
	 case Filter13:
	 
	 //Normalize
	 break;
	 
	 //Shift and Scale
	 case Filter14:
		QTParamLabel1 = new QLabel("Shift Parameter");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Scale Parameter");
		QTParam2 = new QLineEdit();
		QTParam2->setText(tr("1"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
		
	//Sobel Edge Detection	
	case Filter15:	
		QTParamLabel1 = new QLabel("Number of Threads");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("3"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		break;
	
	//Curvature Flow Filter
	case Filter16:
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("0.125"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit();
		QTParam1->setText(tr("5"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
		 
	//Min-Max Flow Filter
	case Filter17:
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit();
		QTParam1->setText(tr("0.125"));
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit();
		QTParam2->setText(tr("5"));
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
	 
		QTParamLabel3 = new QLabel("Stencil Radius");
		QTParam3 = new QLineEdit();
		QTParam3->setText(tr("1"));
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1); 
		 break;
		 
	case Filter18:
	//Resample
	break;
	/* case 18:
		
			QTParamLabel1 = new QLabel("Sigma Min");
			QTParam1 = new QLineEdit(); 
			
			layout->addWidget(QTParamLabel1,1,0);
			layout->addWidget(QTParam1,1,1);
			
			QTParamLabel2 = new QLabel("Sigma Max");
			QTParam2 = new QLineEdit(); 
			layout->addWidget(QTParamLabel2,2,0);
			layout->addWidget(QTParam2,2,1);
			
			QTParamLabel3 = new QLabel("Number of Sigma Steps");
			QTParam3 = new QLineEdit(); 
			layout->addWidget(QTParamLabel3,3,0);
			layout->addWidget(QTParam3,3,1);
			
			QTParamLabel4 = new QLabel("Number Of Iterations");
			QTParam4 = new QLineEdit(); 
			layout->addWidget(QTParamLabel4,4,0);
			layout->addWidget(QTParam4,4,1); 
			break;
	 */		
	case Filter19:
	//Laplacian	
	break;
	case Filter20:
		break;
	} 

	cancelButton = new QPushButton(tr("Cancel"),this);
	connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));
	layout->addWidget(cancelButton,10,1);
	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	layout->addWidget(okButton,10,2);

	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);	
}

int ftkPreprocessDialog::getChannelNumber()
{
	return channelCombo->currentIndex();
}

void ftkPreprocessDialog::accept()
{
	this->doPreprocess();
	QDialog::accept();
}

void ftkPreprocessDialog::doPreprocess(void) 
{		
	std::vector<double> fParams = this->getParams(filtername); 

	ftkPreprocess *ftkpp = new ftkPreprocess();
	ftkpp->myImg = this->myImg;
	ftkpp->channelNumber = this->getChannelNumber();
	ftkpp->filterParams = fParams;
		
	switch( FilterValue[filtername] )
	{
		case Filter1:
			ftkpp->MedianFilter();	
			break;
		case Filter2:
			ftkpp->MeanFilter();	
			break;			
		case Filter3:
			ftkpp->GADiffusion();	
			break;			
		case Filter4:
			ftkpp->CurvAnisotropicDiffusion();	
			break;
		case Filter5:
			ftkpp->SigmoidFilter();	
			break;
		case Filter6:
			ftkpp->GrayscaleErode();	
			break;		
		case Filter7:
			ftkpp->GrayscaleDilate();	
			break;				
		case Filter8:
			ftkpp->GrayscaleOpen();	
			break;
		case Filter9:
			ftkpp->GrayscaleClose();	
			break;
		case Filter10:
			ftkpp->ThreeDSmoothingRGFilter();	
			break;
		case Filter11:
			ftkpp->OpeningbyReconstruction();	
			break;				
		case Filter12:
			ftkpp->ClosingbyReconstruction();	
			break;
		case Filter13:
			ftkpp->NormalizeImage();	
			break;
		case Filter14:
			ftkpp->ShiftScale();	
			break;						
		case Filter15:
			ftkpp->SobelEdgeDetection();	
			break;
		case Filter16:
			ftkpp->CurvatureFlow();	
			break;
		case Filter17:
			ftkpp->MinMaxCurvatureFlow();	
			break;
		case Filter18:
			ftkpp->Resample();	
			break;
		case Filter19:
			ftkpp->LaplacianFilter();	
			break;	
		case Filter20:
			ftkpp->InvertIntensity();
			break;
		default:
			std::cout<<"Something went wrong. Please contact the System Administrator"<<std::endl;
			break;
	}
	delete ftkpp;
}	
	
void ftkPreprocessDialog::InitializeFilters(void) 
{
	FilterValue["Median"] = Filter1;
	FilterValue["Mean"] = Filter2;
	FilterValue["GAnisotropicDiffusion"] = Filter3;
	FilterValue["CAnisotropicDiffusion"] = Filter4;
	FilterValue["Sigmoid"] = Filter5;
	FilterValue["GSErode"] = Filter6;
	FilterValue["GSDilate"] = Filter7;
	FilterValue["GSOpen"] = Filter8;
	FilterValue["GSClose"] = Filter9;
	FilterValue["SRGaussian"] = Filter10;//Needs work
	FilterValue["OpeningbyReconstruction"] = Filter11;
	FilterValue["ClosingbyReconstruction"] = Filter12;
	FilterValue["Normalize"] = Filter13;
	FilterValue["ShiftandScale"] = Filter14;		
	FilterValue["Sobel"] = Filter15;
	FilterValue["CurvatureFlow"] = Filter16;
	FilterValue["MinMax"] = Filter17;
	FilterValue["Resample"] = Filter18;
	FilterValue["Laplacian"] = Filter19;
	FilterValue["Invert"] = Filter20;
}

std::vector<double> ftkPreprocessDialog::getParams(std::string id)
{
	paramVal1 = 0;
	paramVal2 = 0;
	paramVal3 = 0;
	paramVal4 = 0; 
	paramVal5 = 0; 
	
	switch(FilterValue[id])
	{
	case Filter1:case Filter2:case Filter3:case Filter4:case Filter10:case Filter11:case Filter12:case Filter17: 
		paramVal1 =  QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		paramVal3  = QTParam3->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		parameters.push_back(paramVal3);
		break;


	case Filter5: 
		paramVal1 = QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		paramVal3  = QTParam3->text().toDouble();
		paramVal4  = QTParam4->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		parameters.push_back(paramVal3);
		parameters.push_back(paramVal4);
		break;
		
		
	case Filter6: case Filter7:case Filter8:case Filter9:case Filter15: 
		paramVal1 = QTParam1->text().toDouble();
		parameters.push_back(paramVal1);
		break;
		
	case Filter14:case Filter16: 
		
		paramVal1 = QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		break;											

	default:
		break;
	}

	return this->parameters;
}




