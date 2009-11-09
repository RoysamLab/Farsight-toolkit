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

//*****************************************************************************************
// Use the label image to find out if a click occurred on an object, and emit the pixel value there. 
//*****************************************************************************************
#include "LabelImageViewQT.h"

LabelImageViewQT::LabelImageViewQT(QWidget *parent) 
  : QWidget(parent)
{
	this->setupUI();

	channelImg = NULL;
	labelImg = NULL;

	showBounds = true;
	//showIDs = false;
	
	currentScale = 1;					//Image scaling and zooming variables:
	ZoomInFactor = 1.25f;
	ZoomOutFactor = 1 / ZoomInFactor;	
	initGrayscaleColorTable();			//Create a color table of gray values
	backgroundThreshold = 1;			//When adjusting intensities, only change values >= this
	foregroundOffset = 0;				//Offset to ADD to intensity values.

	channelFlags.clear();				//Is channel visible or not?
	channelFlags.push_back(true);
	channelWidget = NULL;

	colorForSelections = Qt::yellow;	//Color used for boundary of selected object
	colorForNormal = Qt::cyan;			//Color for unselected boundaries
	colorForIDs = Qt::green;
	selection = NULL;					//pointer to ObjectSelections class

	setMouseTracking(true);				//Will emit the current mouse position!
	rubberBand = NULL;					//For drawing a box!!
	pointsMode = false;					//For collecting points

}

void LabelImageViewQT::setupUI(void)
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

	imageLabel = new QLabel();
	imageLabel->setBackgroundRole(QPalette::Base);
	imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	imageLabel->setMouseTracking(true);
	imageLabel->setScaledContents(true);	//Image will be stretched to fill label
	imageLabel->resize(0,0);
	scrollArea = new QScrollArea();
	scrollArea->setMouseTracking(true);
	scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
	scrollArea->horizontalScrollBar()->setRange(0,0);
    scrollArea->verticalScrollBar()->setRange(0,0);		
	scrollArea->horizontalScrollBar()->setValue(0);
	scrollArea->verticalScrollBar()->setValue(0);

	QGridLayout *viewerLayout = new QGridLayout();
    viewerLayout->addWidget(scrollArea, 0, 0);
    viewerLayout->addLayout(vsliderLayout, 0, 1);
	viewerLayout->addLayout(hsliderLayout, 1, 0);

	QGridLayout *allLayout = new QGridLayout;
	allLayout->addLayout(viewerLayout,1,0);

	setLayout(allLayout);

	//setup connections
	connect(vSlider, SIGNAL(valueChanged(int)), this, SLOT(sliderChange(int)));
	connect(hSlider, SIGNAL(valueChanged(int)), this, SLOT(sliderChange(int)));
	connect(vSpin, SIGNAL(valueChanged(int)), this, SLOT(spinChange(int)));
	connect(vSpin, SIGNAL(valueChanged(int)), this, SLOT(spinChange(int)));

	setAttribute ( Qt::WA_DeleteOnClose );
	setWindowTitle(tr("Image Browser"));
}

//***************************************************************************************
// This is the original image that the segmentation came from
//***************************************************************************************
void LabelImageViewQT::SetChannelImage(ftk::Image::Pointer img)
{	
	if(!img)	//If img is NULL then I'm trying to remove the channel image
	{
		channelImg = NULL;
		removeChannelWidget();
		refreshDisplayImage();
		return;
	}

	channelImg = img;	//Save the pointer
	if(labelImg)	//If already a label image and sizes don't match, remove it
	{
		if(  labelImg->Size() != channelImg->Size() )
		{
			labelImg = NULL;
			//labelGeometries.clear();
			selection = NULL;
		}
	}
	updateVSlider();
	updateHSlider();
	createChannelWidget();
	//refreshDisplayImage();	//Called after creating the channel widget;
}

//***************************************************************************************
// This is the label image from the segmentation result
//***************************************************************************************
void LabelImageViewQT::SetLabelImage(ftk::Image::Pointer img, ObjectSelection * sels)
{
	if(!img)
	{
		labelImg = NULL;
		refreshDisplayImage();
		return;
	}

	labelImg = img;
	if(channelImg)
	{
		if(  labelImg->Size() != channelImg->Size() )
		{
			removeChannelWidget();
			channelImg = NULL;
		}
	}

	if(sels)
	{
		selection = sels;
		connect(selection, SIGNAL(changed()), this, SLOT(refreshDisplayImage()));
	}

	updateVSlider();
	updateHSlider();
	//refreshFeatures();
	refreshDisplayImage();
}

void LabelImageViewQT::SetBoundsVisible(bool val)
{
	this->showBounds = val;
	refreshDisplayImage();
}

/*
void LabelImageViewQT::SetIDsVisible(bool val)
{
	this->showIDs = val;
	refreshDisplayImage();
}
*/

void LabelImageViewQT::GetBox(void)
{
	if(pointsMode) return;

	if(!rubberBand)	
		rubberBand = new MyRubberBand(this);
}

void LabelImageViewQT::Get2Points(void)
{
	if(rubberBand) return;

	pointsMode = true;
}

void LabelImageViewQT::updateVSlider(void)
{
	const ftk::Image::Info *info;
	if(channelImg)
	{
		info = channelImg->GetImageInfo();	//Get info of new image
	}
	else if(labelImg)
	{
		info = labelImg->GetImageInfo();
	}
	else
	{
		return;
	}

	int numZSlices = (*info).numZSlices;

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

void LabelImageViewQT::updateHSlider(void)
{
	const ftk::Image::Info *info;
	if(channelImg)
	{
		info = channelImg->GetImageInfo();	//Get info of new image
	}
	else if(labelImg)
	{
		info = labelImg->GetImageInfo();
	}
	else
	{
		return;
	}

	int numTSlices = (*info).numTSlices;

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

//***********************************************************************************
// This creates the channel window (widget) that allows for channels to be turned on 
// and off
//***********************************************************************************
void LabelImageViewQT::createChannelWidget()
{
	if(!channelImg) return;

	removeChannelWidget();

	channelWidget = new QWidget(this);
	channelWidget->setWindowTitle(tr("Channels"));
	QVBoxLayout *chLayout = new QVBoxLayout;

	const ftk::Image::Info *info = channelImg->GetImageInfo();
	numChBoxes = (*info).numChannels;

	chBoxes = new QCheckBox * [numChBoxes];
	for (int ch=0; ch < numChBoxes; ++ch)
	{
		std::string chName = (*info).channelNames[ch];
		chBoxes[ch] = new QCheckBox( tr( chName.c_str() ) );
		chBoxes[ch]->setChecked(true);
		connect(chBoxes[ch], SIGNAL(toggled(bool)), this, SLOT(updateChFlags(bool)));
		chLayout->addWidget(chBoxes[ch]);
	}
	chLayout->setSizeConstraint(QLayout::SetFixedSize);
	channelWidget->setLayout(chLayout);
	channelWidget->setWindowFlags( Qt::Tool );
	channelWidget->move(scrollArea->x()+20,scrollArea->y()+20);

	updateChFlags(true);
	if (numChBoxes > 1 && this->isWindow())
		channelWidget->show();
}

void LabelImageViewQT::removeChannelWidget(void)
{
	if(channelWidget)
	{
		channelWidget->close();
		channelWidget = NULL;
		for (int ch=0; ch < numChBoxes; ++ch)
		{
			delete chBoxes[ch]; 
		}
		delete[] chBoxes;
		chBoxes = NULL;
		numChBoxes = 0;
	}
}

void LabelImageViewQT::updateChFlags(bool b)
{
	if(!channelWidget) return;

	channelFlags.clear();
	for (int ch=0; ch<numChBoxes; ++ch)
	{
		channelFlags.push_back(chBoxes[ch]->isChecked());
	}
	refreshDisplayImage();
}

//*******************************************************************
// SLOTS: sliderChange and spinChange
//
// This slot gets called when a slider changes, simply need to update
// the display image
//
// //DON'T USE V UNLESS YOU KNOW WHAT SLIDER IT CAME FROM
//*******************************************************************
void LabelImageViewQT::sliderChange(int v)
{
	vSpin->setValue(vSlider->value());
	hSpin->setValue(hSlider->value());
	refreshDisplayImage();
}
void LabelImageViewQT::spinChange(int v)
{
	vSlider->setValue(vSpin->value());
	hSlider->setValue(hSpin->value());
	refreshDisplayImage();
}

void LabelImageViewQT::update()
{
	refreshDisplayImage();
	QWidget::update();
}

//****************************************************************************************
// Reimplemented moveEvent so that when the image window is moved, the dock widget is 
//  moved with it
//****************************************************************************************
void LabelImageViewQT::moveEvent ( QMoveEvent * event )
{
	QWidget::moveEvent ( event );
	if ( channelWidget )
	{
		int dx = event->pos().x() - event->oldPos().x();
		int dy = event->pos().y() - event->oldPos().y();
		channelWidget->move(channelWidget->x() + dx, channelWidget->y() + dy );
	}
}

void LabelImageViewQT::keyPressEvent(QKeyEvent *event)
 {
	 int key = event->key();
	 if(key == Qt::Key_Equal)
	 {
		 zoom(ZoomInFactor);
	 }
	 else if(key == Qt::Key_Minus)
	 {
		 zoom(ZoomOutFactor);
	 }
	 else if((key >= Qt::Key_0) && (key <= Qt::Key_9))
	 {
		int num = key - 0x30;
		if(numChBoxes > num)
		{
			chBoxes[num]->toggle();
		}
	 }
	 else
	 {
		 QWidget::keyPressEvent(event);
	 }
 }

void LabelImageViewQT::mousePressEvent(QMouseEvent *event)
{	
	QPoint corner = scrollArea->pos();
	origin = event->pos();				// This is this position in the whole widget

	if(rubberBand)	
	{
		rubberBand->setGeometry(QRect(origin,QSize()));
		rubberBand->show();
		return;
	}

	if(pointsMode) return;

	QPoint v_origin = origin - corner;	// This is a local position (in viewport coordinates)

	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;
	int numChannels = (*info).numChannels;
	int currentT = hSpin->value();
	int currentZ = vSpin->value();

	//Compute value in image coordinates and make sure we click within the image:
	int xx = (v_origin.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
	int yy = (v_origin.y() + scrollArea->verticalScrollBar()->value()) / currentScale;
	if( xx < 0 || yy < 0 || xx >= totalWidth || yy >= totalHeight )
		return;

	//Find out the label where I pressed:
	if(!labelImg) return;
	int labelval = 0;
	for(int ch=0; ch < numChannels; ++ch)
	{
		labelval = (int)labelImg->GetPixel(currentT, ch, currentZ, int(yy), int(xx));
		if(labelval != 0) break;
	}
	if(labelval == 0) return;

	Qt::MouseButton button = event->button();
	Qt::KeyboardModifiers modifiers = event->modifiers();
	if(button == Qt::LeftButton && modifiers == Qt::ControlModifier)
	{
		if(selection)
			selection->toggle(labelval);
	}
	else if(button == Qt::LeftButton && modifiers == Qt::NoModifier)
	{
		if(selection)
			selection->select(labelval);
	}
	else if(button == Qt::RightButton)
	{
		//This shows the tooltip at the global position (screen coordinates)
		QToolTip::showText(event->globalPos(), QString("ID: ") + QString::number(labelval) );
	}
}


void LabelImageViewQT::mouseMoveEvent(QMouseEvent *event)
{
	QPoint pos = event->pos();
	if(rubberBand)
	{
		if(pos.x() <= origin.x() || pos.y() <= origin.y())
			rubberBand->setGeometry(QRect(origin,QSize()));
		else
			rubberBand->setGeometry(QRect(origin, pos).normalized());
	}

	QPoint corner = scrollArea->pos();
	QPoint v_pos = pos - corner;		// This is a local position (in viewport coordinates)
	int xx = (v_pos.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
	int yy = (v_pos.y() + scrollArea->verticalScrollBar()->value()) / currentScale;

	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();	//Get info of new image
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;
	int currentZ = vSpin->value();

	if( xx>=0 && xx<totalWidth && yy>=0 && yy<totalHeight )
		emit mouseAt(xx, yy, currentZ);
}

void LabelImageViewQT::mouseReleaseEvent(QMouseEvent *event)
{
	if(rubberBand || pointsMode)
	{
		QPoint corner = scrollArea->pos();
		QPoint pos = event->pos() - corner;		// This is a local position (in viewport coordinates)
		QPoint org = origin - corner;
		int x1 = (org.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
		int y1 = (org.y() + scrollArea->verticalScrollBar()->value()) / currentScale;
		int x2 = (pos.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
		int y2 = (pos.y() + scrollArea->verticalScrollBar()->value()) / currentScale;

		if(rubberBand)
		{
			delete rubberBand;
			rubberBand = NULL;
			origin = QPoint();
			emit boxDrawn(x1, y1, x2, y2, vSpin->value());
		}
		else //pointsMode
		{
			if( abs(x1-x2) < 5 && abs(y1-y2) < 5 )		//I haven't moved for during this click
			{
				if(origin3.size() == 0)
				{
					origin3.push_back(x2);
					origin3.push_back(y2);
					origin3.push_back(vSpin->value());
				}
				else
				{
					x1 = origin3.at(0);
					y1 = origin3.at(1);
					int z1 = origin3.at(2);
					origin3.clear();
					pointsMode = false;
					emit pointsClicked(x1,y1,z1,x2,y2,vSpin->value());
				}
			}
		}
	}
}

//*****************************************************************************************
// change the currentScale
//*****************************************************************************************
void LabelImageViewQT::zoom(double zf)
{
	double newScale = currentScale * zf;
	if(zf > 1 && newScale > 5)
		return;
	if(zf < 1 && newScale < 0.20)
		return;

	currentScale = newScale;
	refreshDisplayImage();
}

 //***************************************************************************************************
//The function redraws the plot onto the off-screen pixmap and updates the display.
//***************************************************************************************************
void LabelImageViewQT::refreshDisplayImage()
{
	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();	//Get info of new image
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;

	displayImage = QImage(totalWidth, totalHeight, QImage::Format_ARGB32);	
	displayImage.fill(qRgb(0,0,0));
	QPainter painter(&displayImage);
	painter.setCompositionMode(QPainter::CompositionMode_Plus);
	drawImage(&painter);
	painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
	drawBoundaries(&painter);
	//drawObjectIDs(&painter);

	//Do zooming:
	int oldX = scrollArea->horizontalScrollBar()->value();
	int oldY = scrollArea->verticalScrollBar()->value();

	imageLabel->setPixmap(QPixmap::fromImage(displayImage));
	imageLabel->adjustSize();
	if(currentScale != 1)
	{
		QSize newSize = displayImage.size()*abs(currentScale);
		imageLabel->resize(newSize);
	}

	scrollArea->horizontalScrollBar()->setValue(oldX);
	scrollArea->verticalScrollBar()->setValue(oldY);

}

//***************************************************************************************
// The drawImage function draws the iamge
//***************************************************************************************
void LabelImageViewQT::drawImage(QPainter *painter)
{
	if(!channelImg)
		return;

	int currentZ = vSpin->value();
	int currentT = hSpin->value();

	const ftk::Image::Info *info = channelImg->GetImageInfo();
	for (int i=0; i < (*info).numChannels; i++)
	{
		if (channelFlags[i])
		{
			QImage gray((*info).numColumns, (*info).numRows, QImage::Format_ARGB32_Premultiplied);
			std::vector<unsigned char> color = (*info).channelColors[i];
			gray.fill(qRgb(color[0],color[1],color[2]));
			unsigned char * p = channelImg->GetSlicePtr<unsigned char>( currentT, i, currentZ );
			if(p)
			{
				//Get the image:
				QImage img(p, (*info).numColumns, (*info).numRows, (*info).numColumns, QImage::Format_Indexed8);
				if(foregroundOffset == 0)
				{
					gray.setAlphaChannel(img);	//Set it to the alpha channel
				}
				else
				{
					QImage img2 = img.copy();
					scaleIntensity( &img2, backgroundThreshold, foregroundOffset );
					gray.setAlphaChannel( img2 );	//Set it to the alpha channel
				}
			}
			painter->drawImage(0,0,gray);
		}
	}
}

//img must be of QImage::Format_Indexed8
void LabelImageViewQT::scaleIntensity(QImage *img, int threshold, int offset)
{
	if( threshold < 0) threshold = 0;

	img->setColorTable(grayscaleColorTable);
	int old_v, new_v;
	for(int c=0; c<img->width(); ++c)
	{
		for(int r=0; r<img->height(); ++r)
		{
			old_v = img->pixelIndex(r,c);
			if( old_v >= threshold)
			{
				new_v = old_v + offset;
				if(new_v > 255) new_v=255;
				else if(new_v < 0) new_v=0;
				img->setPixel(r,c, new_v);
			}
		}
	}
}

void LabelImageViewQT::AdjustImageIntensity(void)
{
	IntensityDialog *dialog = new IntensityDialog(backgroundThreshold, foregroundOffset, this);
	connect(dialog, SIGNAL(valuesChanged(int,int)), this, SLOT(adjustImageIntensity(int,int)));
	dialog->show();
}

void LabelImageViewQT::adjustImageIntensity(int threshold, int offset)
{
	backgroundThreshold = threshold;
	foregroundOffset = offset;
	this->refreshDisplayImage();
}

void LabelImageViewQT::initGrayscaleColorTable(void)
{
	grayscaleColorTable.clear();
	for(int i=0; i<256; ++i)
	{
		grayscaleColorTable.append(qRgb(i,i,i));
	}
}

//***************************************************************************************
// The drawBoundaries function draws the boundaries of the label image
//***************************************************************************************
void LabelImageViewQT::drawBoundaries(QPainter *painter)
{
	if(!labelImg)
		return;

	if(!showBounds)
		return;

	const ftk::Image::Info *info = labelImg->GetImageInfo();
	int chs = (*info).numChannels;
	int h = (*info).numRows;
	int w = (*info).numColumns;
	int currentZ = vSpin->value();
	int currentT = hSpin->value();

	for(int ch = 0; ch < chs; ++ch)
	{
		std::vector<unsigned char> color = (*info).channelColors[ch];
		QColor qcolor(color[0],color[1],color[2]); 

		int v, v1, v2, v3, v4;
		for(int i=1; i < h-1; i++)
		{
			for(int j=1; j < w-1; j++)
			{
				v = (int)labelImg->GetPixel(currentT, ch, currentZ, i, j);
				if (v > 0)
				{
					v1 = (int)labelImg->GetPixel(currentT, ch, currentZ, i, j+1);
					v2 = (int)labelImg->GetPixel(currentT, ch, currentZ, i+1, j);
					v3 = (int)labelImg->GetPixel(currentT, ch, currentZ, i, j-1);
					v4 = (int)labelImg->GetPixel(currentT, ch, currentZ, i-1, j);
					if(v!=v1 || v!=v2 || v!=v3 || v!=v4)
					{
						painter->setPen(qcolor);
						if(selection) 
							if(selection->isSelected(v))
								painter->setPen(colorForSelections);
						painter->drawPoint(j,i);
					}
				}
			}
		}
	}
}

/* SHOULD REWRITE THIS SO IT DOESN'T NEED CENTROID, JUST DRAWS ID
void LabelImageViewQT::drawObjectIDs(QPainter *painter)
{
	if(!showIDs) return;
	if(!labelImg) return;
	if(labelGeometries.size() == 0) return;

	int currentZ = vSpin->value();

	//Iterate through each object and write its id at its centroid.
	std::map<int,LabelGeometry>::iterator it;
	for ( it = labelGeometries.begin() ; it != labelGeometries.end(); ++it )
	{
		int id = (*it).first;
		LabelGeometry geo = (*it).second;
		if ( (currentZ >= geo.BoundingBox[4]) && (currentZ <= geo.BoundingBox[5]) )
		{
			painter->setPen(colorForIDs);
			painter->drawText(geo.Centroid[0], geo.Centroid[1], QString::number(id));
		}
	}
}
*/
/*
void LabelImageViewQT::refreshFeatures(void)
{
	if(!this->labelImg)
		return;

	labelGeometries.clear();

	//Compute features:
	typedef itk::Image<unsigned char, 3> ChannelImgType;
	typedef itk::Image<unsigned short, 3> LabelImgType;
	typedef itk::LabelGeometryImageFilter< LabelImgType, ChannelImgType > FeatureCalcType;
	FeatureCalcType::Pointer labFilter = FeatureCalcType::New();
	labFilter->SetInput( labelImg->GetItkPtr<unsigned short>(0,0) );
	labFilter->Update();

	//Get the centroid and bounding box:
	std::vector< FeatureCalcType::LabelPixelType > ls = labFilter->GetLabels();
	for (int l = 0; l < (int)ls.size(); ++l)
	{
		int label = (int)ls.at(l);
		LabelGeometry geo;
		
		FeatureCalcType::LabelPointType c = labFilter->GetCentroid( label );
		for (unsigned int i = 0; i < 3; ++i)
		{
			geo.Centroid[i] = float( c[i] );
		}
		FeatureCalcType::BoundingBoxType bbox = labFilter->GetBoundingBox( label );
		for (unsigned int i = 0; i < 3*2; ++i)
		{
			geo.BoundingBox[i] = float( bbox[i] );
		}
		labelGeometries[l] = geo; 
	}
}
*/

MyRubberBand::MyRubberBand(QWidget * p)
: QWidget(p)
{
	this->setWindowFlags( Qt::FramelessWindowHint | this->windowFlags() );
	this->setMouseTracking(true);
}

void MyRubberBand::paintEvent(QPaintEvent *event)
{
	QPainter painter(this);
	painter.setPen(Qt::blue);
	painter.setBrush(Qt::NoBrush);
	QRect r = this->rect();
	r.setWidth(r.width()-1);
	r.setHeight(r.height()-1);
	painter.drawRect(r.normalized());
}

void MyRubberBand::mouseMoveEvent(QMouseEvent * event)
{
	this->setGeometry(QRect(this->pos(), this->pos()+ event->pos()).normalized());
}

IntensityDialog::IntensityDialog(int threshold, int offset, QWidget *parent)
: QDialog(parent)
{
	QLabel * header = new QLabel(tr("Adjust values for Intensity shifting:"));
	QLabel * label1 = new QLabel(tr("Threshold: "));
	QString thresholdMessage(tr("Any intensity value below this threshold with not be changed"));
	label1->setToolTip(thresholdMessage);
	thresholdSpin = new QSpinBox();
	thresholdSpin->setRange(1,255);
	thresholdSpin->setValue(threshold);
	thresholdSpin->setToolTip(thresholdMessage);
	connect(thresholdSpin, SIGNAL(valueChanged(int)), this, SLOT(changeThreshold(int)));

	QLabel * label2 = new QLabel(tr("Offset: "));
	QString offsetMessage(tr("This value will be added to the intensity values of each channel"));
	label2->setToolTip(offsetMessage);
	offsetSpin = new QSpinBox();
	offsetSpin->setRange(-255,255);
	offsetSpin->setSingleStep(5);
	offsetSpin->setValue(offset);
	offsetSpin->setToolTip(offsetMessage);
	connect(offsetSpin, SIGNAL(valueChanged(int)), this, SLOT(changeOffset(int)));

	hideButton = new QPushButton(tr("DONE"));
	connect(hideButton, SIGNAL(clicked()), this, SLOT(close()));

	QGridLayout * layout = new QGridLayout();
	layout->addWidget(header,0,0,1,2);
	layout->addWidget(label1,1,0,1,1);
	layout->addWidget(thresholdSpin,1,1,1,1);
	layout->addWidget(label2,2,0,1,1);
	layout->addWidget(offsetSpin,2,1,1,1);
	layout->addWidget(hideButton,3,1,1,1);
	this->setLayout(layout);
	this->setWindowTitle(tr("Adjust Image Intensity"));
	this->setModal(false);
	this->setAttribute(Qt::WA_DeleteOnClose, true);
}

void IntensityDialog::changeThreshold(int v)
{
	emit valuesChanged(v, offsetSpin->value());
}

void IntensityDialog::changeOffset(int v)
{
	emit valuesChanged(thresholdSpin->value(), v);
}
