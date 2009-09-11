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

//**************************************************************************************
// This class provides the widgets around segmentationview to allow for multiple 
// dimension images to be displayed.  It passes the current Z,T, and channel info
// to the view, so theview can display the currect image
//**************************************************************************************

#include "SegmentationWindow.h"

//Constructor
SegmentationWindow::SegmentationWindow(QWidget *parent)
  : QWidget(parent)
{
	//Setup Vertical slider widgets
    vSlider = new QSlider();
    vSlider->setOrientation(Qt::Vertical);
	vSlider->setDisabled(true);
	vSlider->setRange(0,0);
	vSlider->setValue(0);

	vSpin = new QSpinBox();
	vSpin->setDisabled(true);
	vSpin->setRange(0,0);
	vSpin->setValue(0);
	vSpin->resize( vSpin->minimumSizeHint() );

	vLabel = new QLabel("z");
	vLabel->setDisabled(true);

	QGridLayout *vsliderLayout = new QGridLayout;
	vsliderLayout->addWidget(vSpin,0,0,1,2);
	vsliderLayout->addWidget(vSlider,1,0,1,1);
	vsliderLayout->addWidget(vLabel,1,1,1,1);

	//Setup Horizontal slider widgets
    hSlider = new QSlider();
    hSlider->setOrientation(Qt::Horizontal);
	hSlider->setDisabled(true);
	hSlider->setRange(0,0);
	hSlider->setValue(0);

	hSpin = new QSpinBox();
	hSpin->setDisabled(true);
	hSpin->setRange(0,0);
	hSpin->setValue(0);

	hLabel = new QLabel("t");
	hLabel->setDisabled(true);

	QHBoxLayout *hsliderLayout = new QHBoxLayout;
	hsliderLayout->addWidget(hSpin);
	hsliderLayout->addWidget(hLabel);
	hsliderLayout->addWidget(hSlider);

	//Setup viewer window
    segview = new SegmentationView();

	QGridLayout *viewerLayout = new QGridLayout();
    viewerLayout->addWidget(segview, 0, 0);
    viewerLayout->addLayout(vsliderLayout, 0, 1);
	viewerLayout->addLayout(hsliderLayout, 1, 0);

	QGridLayout *allLayout = new QGridLayout;
	allLayout->addLayout(viewerLayout,1,0);

	setLayout(allLayout);

	//setup connections
	connect(vSlider, SIGNAL(valueChanged(int)), segview, SLOT(setZ(int)));
	connect(hSlider, SIGNAL(valueChanged(int)), segview, SLOT(setT(int)));
	connect(vSlider, SIGNAL(valueChanged(int)), vSpin, SLOT(setValue(int)));
	connect(hSlider, SIGNAL(valueChanged(int)), hSpin, SLOT(setValue(int)));
	connect(vSpin, SIGNAL(valueChanged(int)), vSlider, SLOT(setValue(int)));
	connect(hSpin, SIGNAL(valueChanged(int)), hSlider, SLOT(setValue(int)));

	connect(segview, SIGNAL(goToZ(int)),vSlider,SLOT(setValue(int)));

	setAttribute ( Qt::WA_DeleteOnClose );
	setWindowTitle(tr("Image Browser"));

	numZSlices = 0;
	numTSlices = 0;
	channelWidget = NULL;
}

void SegmentationWindow::SetBoundsVisible(bool val)
{
	segview->setBoundsVisible(val);
}

void SegmentationWindow::SetIDsVisible(bool val)
{
	segview->setIDsVisible(val);
}

void SegmentationWindow::closeEvent(QCloseEvent *event)
{
	emit closing(this);
	event->accept();
} 

//***********************************************************************************
// channelWidget is the channel selection widget.  We don't show it until we show
// this window.
//***********************************************************************************
void SegmentationWindow::showEvent ( QShowEvent * event )
{
	QWidget::showEvent(event);
	if (channelWidget)
	{
		channelWidget->show();
	}
}

//****************************************************************************************
// Reimplemented moveEvent so that when the image window is moved, the dock widget is 
//  moved with it
//****************************************************************************************
void SegmentationWindow::moveEvent ( QMoveEvent * event )
{
	QWidget::moveEvent ( event );
	if ( channelWidget )
	{
		int dx = event->pos().x() - event->oldPos().x();
		int dy = event->pos().y() - event->oldPos().y();
		channelWidget->move(channelWidget->x() + dx, channelWidget->y() + dy );
	}
}

//***********************************************************************************
// This creates the channel window (widget) that allows for channels to be turn on 
// and off
//***********************************************************************************
void SegmentationWindow::createChannelWindow(ftk::Image::Pointer image)
{
	closeChannelWindow();

	channelWidget = new QWidget(this);
	channelWidget->setWindowTitle(tr("Channels"));
	QVBoxLayout *chLayout = new QVBoxLayout;

	const ftk::Image::Info *info = image->GetImageInfo();

	chBoxes = new QCheckBox * [numChannels];
	for (int ch=0; ch < numChannels; ++ch)
	{
		std::string chName = (*info).channelNames[ch];
		chBoxes[ch] = new QCheckBox( tr( chName.c_str() ) );
		chBoxes[ch]->setChecked(true);
		connect(chBoxes[ch], SIGNAL(clicked(bool)), this, SLOT(updateChannels(bool)));
		chLayout->addWidget(chBoxes[ch]);
	}
	chLayout->setSizeConstraint(QLayout::SetFixedSize);
	channelWidget->setLayout(chLayout);
	channelWidget->setWindowFlags( Qt::Tool );
	channelWidget->move(this->x()+20,this->y()+100);
	channelWidget->show();
}

void SegmentationWindow::closeChannelWindow(void)
{
	if(channelWidget)
	{
		channelWidget->close();
		channelWidget = NULL;
		for (int ch=0; ch < numChannels; ++ch)
		{
			delete chBoxes[ch]; 
		}
		delete[] chBoxes;
		chBoxes = NULL;
	}
}

void SegmentationWindow::SetModels(SegmentationModel *sModel)
{
	segview->setModels(sModel);
}

//***********************************************************************************
// ChannelImage is the data image (ftkImage) that is was segmented
//***********************************************************************************
//void SegmentationWindow::AddChannelImage( QString &fileName )
void SegmentationWindow::SetChannelImage(ftk::Image::Pointer image)
{
	ftk::Image::Info info;
	if(!image)
	{
		info.numZSlices = 0;
		info.numTSlices = 0;
		info.numChannels = 0;
	}
	else
	{
		info = *(image->GetImageInfo());
	}

	if( numZSlices != info.numZSlices )
	{
		//this->resize((*info).numColumns+60,(*info).numRows+85); 
		numZSlices = info.numZSlices;
		updateVSlider();
	}
	
	if( numTSlices != info.numTSlices )
	{
		numTSlices = info.numTSlices;
		updateHSlider();
	}

	if (numChannels != info.numChannels )
	{
		numChannels = info.numChannels;
		if (numChannels > 1)
		{
			createChannelWindow(image);
			updateChannels();
		}
		else
		{
			closeChannelWindow();
			std::vector<bool> flags;
			flags.clear(); 
			flags.push_back(true);
			segview->setChannelFlags(flags);
		}
	}
	segview->setChannelImage(image);
}

//***********************************************************************************
// The label image is the segmentation result
//***********************************************************************************
//void SegmentationWindow::AddLabelImage( QString &fileName )
void SegmentationWindow::SetLabelImage( ftk::Image::Pointer image)
{
	ftk::Image::Info info;
	if(!image)
	{
		info.numZSlices = 0;
		info.numTSlices = 0;
		info.numChannels = 0;
	}
	else
	{
		info = *(image->GetImageInfo());
	}

	if(!(segview->getChannelImage()))
	{
		//this->resize( (*info).numColumns+60, (*info).numRows+85);
		
		if ( numZSlices != info.numZSlices )
		{
			numZSlices = info.numZSlices;
			updateVSlider();
		}
		if (numTSlices != info.numTSlices )
		{
			numTSlices = info.numTSlices;
			updateHSlider();
		}
	}
	segview->setLabelImage(image);
}

void SegmentationWindow::updateVSlider(void)
{
	vSlider->setRange(0,numZSlices-1);
	vSpin->setRange(0,numZSlices-1);
	vSlider->setValue(0);
	vSpin->setValue(0);
	if (numZSlices > 1)
	{
		vSlider->setEnabled(true);
		vSpin->setEnabled(true);
		vLabel->setEnabled(true);
	}
	else
	{
		vSlider->setEnabled(false);
		vSpin->setEnabled(false);
		vLabel->setEnabled(false);
	}
}

void SegmentationWindow::updateHSlider(void)
{
	hSlider->setRange(0,numTSlices-1);
	hSlider->setValue(0);
	hSpin->setRange(0,numTSlices-1);
	hSpin->setValue(0);
	if (numTSlices > 1)
	{
		hSlider->setEnabled(true);
		hSpin->setEnabled(true);
		hLabel->setEnabled(true);
	}
	else
	{
		hSlider->setEnabled(false);
		hSpin->setEnabled(false);
		hLabel->setEnabled(false);
	}
}

void SegmentationWindow::updateChannels(void)
{
	std::vector<bool> flags;
	flags.clear();
	for (int ch=0; ch<numChannels; ++ch)
	{
		flags.push_back(chBoxes[ch]->isChecked());
	}
	segview->setChannelFlags(flags);
}
