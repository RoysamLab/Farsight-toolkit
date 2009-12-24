/*
 *  ftkPreprocessDialog.cpp
 *  Farsight
 *
 *  Created by RAGHAV on 12/3/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "ftkPreprocessDialog.h"
#include "ftkPreprocess.h"


//******************************************************************************************
//******************************************************************************************
// A dialog to get the paramaters file for the preprocessing to use and specify the channel 
// if image has more than one:
//******************************************************************************************

ftkPreprocessDialog::ftkPreprocessDialog(QVector<QString> channels, unsigned char id, QWidget *parent)
: QDialog(parent)
{

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

	 switch(id)
	{
	
	//Mean and Median Filter
	case 0:case 1:  
		QTParamLabel1 = new QLabel("Window Size - X");
		QTParam1 = new QLineEdit(); 
		
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Window Size - Y");
		QTParam2 = new QLineEdit(); 
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Window Size - Z");
		QTParam3 = new QLineEdit(); 
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		 break;


	//G and C Anisotropic Diffusion
	case 2:case 3: 
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit(); 
		
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit(); 
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Conductance");
		QTParam3 = new QLineEdit(); 
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		 break;
		 
	//Sigmoid Filters		 
	case 4: 
		QTParamLabel1 = new QLabel("Alpha");
		QTParam1 = new QLineEdit(); 
		
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Beta");
		QTParam2 = new QLineEdit(); 
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Output Minimum");
		QTParam3 = new QLineEdit(); 
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1);
		
		QTParamLabel4 = new QLabel("Output Maximum");
		QTParam4 = new QLineEdit(); 
		layout->addWidget(QTParamLabel4,4,0);
		layout->addWidget(QTParam4,4,1);

		 break;
		 
	//Grayscale Morphological Filters	 
	case 5:case 6:case 7:case 8: 
		QTParamLabel1 = new QLabel("Radius");
		QTParam1 = new QLineEdit(); 
		
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		 
		 break;

	//Smoothing Recursive Gaussian	 
	case 9:
		
		QTParamLabel1 = new QLabel("Sigma");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number of Threads");
		QTParam2 = new QLineEdit();  
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Normalize across Scale ?");
		QTParam5 = new QCheckBox();  
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam5,3,1);
		 break;
	
	//	Opening/Closing by Reconstruction 	 	 
	case 10:case 11:
		
		QTParamLabel1 = new QLabel("Radius");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel3 = new QLabel("Preserve Intensities ?");
		QTParam5 = new QCheckBox();  
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		
		QTParamLabel3 = new QLabel("Fully Connected ?");
		QTParam6 = new QCheckBox();  
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam6,3,1);
		 break;
	 
	 case 12:
	 
	 //Normalize
	 break;
	 
	 //Shift and Scale
	 case 13:
		QTParamLabel1 = new QLabel("Shift Parameter");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Scale Parameter");
		QTParam2 = new QLineEdit();  
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
		
	//Sobel Edge Detection	
	case 14:	
		QTParamLabel1 = new QLabel("Number of Threads");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		break;
	
	//Curvature Flow Filter
	case 15:
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit(); 
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
		 
	//Min-Max Flow Filter
	case 16:
		QTParamLabel1 = new QLabel("TimeStep");
		QTParam1 = new QLineEdit(); 
		layout->addWidget(QTParamLabel1,1,0);
		layout->addWidget(QTParam1,1,1);
		
		QTParamLabel2 = new QLabel("Number Of Iterations");
		QTParam2 = new QLineEdit(); 
		layout->addWidget(QTParamLabel2,2,0);
		layout->addWidget(QTParam2,2,1);
		break;
	 
		QTParamLabel3 = new QLabel("Stencil Radius");
		QTParam3 = new QLineEdit(); 
		layout->addWidget(QTParamLabel3,3,0);
		layout->addWidget(QTParam3,3,1); 
		 break;
		 
	case 17:
	//Resample
	
	case 18:
	
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
		
	case 19:
	//Laplacian	
	break;	
	 
	} 

	okButton = new QPushButton(tr("OK"),this);
	connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
	layout->addWidget(okButton,10,1);	
	Qt::WindowFlags flags = this->windowFlags();
	flags &= ~Qt::WindowContextHelpButtonHint;
	this->setWindowFlags(flags);
			
}



int ftkPreprocessDialog::getChannelNumber()
{
	return channelCombo->currentIndex();
}


ftk::Image::Pointer ftkPreprocessDialog::getImage() 
{

		vector<double> fParams = this->getParams(16); 
		ftkPreprocess *ftkpp = new ftkPreprocess();
		ftkpp->myImg = this->myImg;
		ftkpp->filterParams = fParams;
		myImg = ftkpp->MinMaxCurvatureFlow();
		return myImg;

}	
	

vector<double> ftkPreprocessDialog::getParams(unsigned char id)
{

	paramVal1 = 0;
	paramVal2 = 0;
	paramVal3 = 0;
	paramVal4 = 0; 
	paramVal5 = 0; 

	
	 switch(id)
	{
	case 0:case 1:case 2:case 3:case 9:case 10:case 11:case 16: 
		paramVal1 = QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		paramVal3  = QTParam3->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		parameters.push_back(paramVal3);
		break;


	case 4:case 18: 
		paramVal1 = QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		paramVal3  = QTParam3->text().toDouble();
		paramVal4  = QTParam4->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		parameters.push_back(paramVal3);
		parameters.push_back(paramVal4);
		break;
		
		
	case 5: case 6:case 7:case 8:case 14: 
		paramVal1 = QTParam1->text().toDouble();
		
		parameters.push_back(paramVal1);

		break;
		
	case 13:case 15: 
		paramVal1 = QTParam1->text().toDouble();
		paramVal2  = QTParam2->text().toDouble();
		
		parameters.push_back(paramVal1);
		parameters.push_back(paramVal2);
		
		break;											

	}

return this->parameters;

}

