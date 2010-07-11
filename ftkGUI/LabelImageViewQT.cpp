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

LabelImageViewQT::LabelImageViewQT(QMap<QString, QColor> * new_colorItemsMap, QWidget *parent) 
  : QWidget(parent)
{
	this->readSettings();
	this->setupUI();

	//Setup all the color needs:
	centroidColorTable = this->CreateColorTable();
	if(!new_colorItemsMap)
		this->SetColorsToDefaults();
	else
		this->colorItemsMap = new_colorItemsMap;

	channelImg = NULL;
	labelImg = NULL;
	centerMap = NULL;
	bBoxMap = NULL;
	
	currentScale = 1;					//Image scaling and zooming variables:
	ZoomInFactor = 1.25f;
	ZoomOutFactor = 1 / ZoomInFactor;	
	initGrayscaleColorTable();			//Create a color table of gray values
	backgroundThreshold = 1;			//When adjusting intensities, only change values >= this
	foregroundOffset = 0;				//Offset to ADD to intensity values.

	channelFlags.clear();				//Is channel visible or not?
	//channelFlags.push_back(true);

	classMap.clear();

	selection = NULL;					//pointer to ObjectSelections class

	setMouseTracking(true);				//Will emit the current mouse position!
	rubberBand = NULL;					//For drawing a box!!
	pointsMode = false;					//For collecting points
	roiMode = false;					//For drawing Region of Interest (ROI)
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
	connect(hSpin, SIGNAL(valueChanged(int)), this, SLOT(spinChange(int)));

	setAttribute ( Qt::WA_DeleteOnClose );
	setWindowTitle(tr("Image Browser"));
}

//******************************************************************************
//Write/read settings functions
//******************************************************************************
void LabelImageViewQT::readSettings()
{
	QSettings settings;
	settings.beginGroup("labelImageViewQT");
	showBounds = settings.value("showBounds", true).toBool();
	showCrosshairs = settings.value("showCrosshairs", true).toBool();
	showIDs = settings.value("showIDs", false).toBool();
	showCentroids = settings.value("showCentroids", false).toBool();
	showROI = settings.value("showROI", false).toBool();
	settings.endGroup();
}

//Destructor:
LabelImageViewQT::~LabelImageViewQT()
{
	writeSettings();
}

void LabelImageViewQT::writeSettings()
{
	QSettings settings;
	settings.beginGroup("labelImageViewQT");
	settings.setValue("showBounds", showBounds);
	settings.setValue("showCrosshairs", showCrosshairs);
	settings.setValue("showIDs", showIDs);
	settings.setValue("showCentroids", showCentroids);
	//settings.setValue("showROI", showROI);
	settings.endGroup();
}

void LabelImageViewQT::SetColorsToDefaults(void)
{
	colorItemsMap = new QMap<QString, QColor>();
	(*colorItemsMap)["Selected Objects"] = Qt::yellow;		//Color used for boundary of selected object
	(*colorItemsMap)["Object Boundaries"] = Qt::cyan;		//Color for unselected boundaries
	(*colorItemsMap)["Object IDs"] = Qt::green;				//Color for ID numbers
	(*colorItemsMap)["ROI Boundary"] = Qt::gray;	
}
//***************************************************************************************
// This is the original image that the segmentation came from
//***************************************************************************************
void LabelImageViewQT::SetChannelImage(ftk::Image::Pointer img)
{	
	if(!img)	//If img is NULL then I'm trying to remove the channel image
	{
		channelImg = NULL;
		refreshBaseImage();
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
			refreshBoundsImage();
		}
	}
	if(!labelImg)
	{
		updateVSlider();
		updateHSlider();
	}

	initChannelFlags();
	refreshBaseImage();
}

//***************************************************************************************
// This is the label image from the segmentation result
//***************************************************************************************
void LabelImageViewQT::SetLabelImage(ftk::Image::Pointer img, ObjectSelection * sels)
{
	if(!img)
	{
		labelImg = NULL;
		refreshBoundsImage();
		return;
	}

	labelImg = img;
	if(channelImg)
	{
		if(  labelImg->Size() != channelImg->Size() )
		{
			channelImg = NULL;
			refreshBaseImage();
		}
	}

	if(!channelImg)
	{
		updateVSlider();
		updateHSlider();
	}

	if(sels)
	{
		selection = sels;
		connect(selection, SIGNAL(changed()), this, SLOT(selectionChange()));
	}

	//refreshFeatures();
	refreshBoundsImage();
}

void LabelImageViewQT::SetClassMap(vtkSmartPointer<vtkTable> table, std::vector<std::string> columns){
	classMap.clear();

	vtkAbstractArray * output = table->GetColumnByName(columns.at(columns.size()-1).c_str());
	if(output == 0)
		return;

	for(int i=0; i<table->GetNumberOfRows(); ++i){
		classMap[table->GetValue(i,0).ToInt()] = table->GetValueByName(i,columns.at(columns.size()-1).c_str()).ToInt();
	}

	refreshBoundsImage();
}

void LabelImageViewQT::SetROIMaskImage(QImage img)
{
	roiImage = img;
	this->SetROIVisible(true);
}

void LabelImageViewQT::SetROIVisible(bool val)
{
	this->showROI = val;
	refreshBoundsImage();
}

void LabelImageViewQT::SetCrosshairsVisible(bool val)
{
	this->showCrosshairs = val;
	refreshBoundsImage();
}

void LabelImageViewQT::SetBoundsVisible(bool val)
{
	this->showBounds = val;
	refreshBoundsImage();
}

void LabelImageViewQT::SetIDsVisible(bool val)
{
	this->showIDs = val;
	refreshBoundsImage();
}

void LabelImageViewQT::SetCentroidsVisible(bool val)
{
	this->showCentroids = val;
	refreshBoundsImage();
}


void LabelImageViewQT::SetCenterMapPointer(std::map<int, ftk::Object::Point> * cMap)
{
	centerMap = cMap;
	this->update();
}

void LabelImageViewQT::SetBoundingBoxMapPointer(std::map<int, ftk::Object::Box> * bMap)
{
	bBoxMap = bMap;
	this->update();
}

void LabelImageViewQT::ClearGets(void)
{
	//stop getting points:
	pointsMode = false;
	origin3.clear();

	//stop getting ROI:
	roiMode = false;
	roiPoints.clear();

	//stop getting box:
	if(rubberBand)
	{
		delete rubberBand;
		rubberBand = NULL;
		origin = QPoint();
	}
}

void LabelImageViewQT::GetBox(void)
{
	if(pointsMode) return;
	if(roiMode) return;

	if(!rubberBand)	
		rubberBand = new MyRubberBand(this);
}

void LabelImageViewQT::Get2Points(void)
{
	if(rubberBand) return;
	if(roiMode) return;

	pointsMode = true;
}

void LabelImageViewQT::GetROI(void)
{
	if(rubberBand) return;
	if(pointsMode) return;

	roiMode = true;
}

void LabelImageViewQT::SaveDisplayImageToFile(QString fileName)
{
	bool ok = displayImage.save( fileName );
	if(!ok)
		QMessageBox::warning(this, tr("Save Failure"), tr("Image Failed to Save"));
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

void LabelImageViewQT::initChannelFlags()
{
	channelFlags.clear();
	std::vector<std::string> channel_names = channelImg->GetChannelNames();
	for (int ch=0; ch<(int)channel_names.size(); ++ch)
	{
		channelFlags.push_back(true);
	}
}

void LabelImageViewQT::SetChannelFlags(std::vector<bool> ch_fg)
{
	channelFlags = ch_fg;
	refreshBaseImage();
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
}
void LabelImageViewQT::spinChange(int v)
{
	vSlider->setValue(vSpin->value());
	hSlider->setValue(hSpin->value());
	refreshBaseImage();		//Only need this in either slider or spin changes!!!
	refreshBoundsImage();
}

void LabelImageViewQT::update()
{
	refreshBaseImage();
	refreshBoundsImage();
	QWidget::update();
}

//****************************************************************************************
// Reimplemented moveEvent so that when the image window is moved, the dock widget is 
//  moved with it
//****************************************************************************************
void LabelImageViewQT::moveEvent ( QMoveEvent * event )
{
	QWidget::moveEvent ( event );
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
		if((int)channelFlags.size() > num)
		{
			channelFlags.at(num) = !channelFlags.at(num);
			refreshBaseImage();
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
	if(roiMode) return;

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
	if(rubberBand || pointsMode || roiMode)
	{
		QPoint corner = scrollArea->pos();
		QPoint pos = event->pos() - corner;		// This is a local position (in viewport coordinates)
		QPoint org = origin - corner;			// This the last click position
		int x1 = (org.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
		int y1 = (org.y() + scrollArea->verticalScrollBar()->value()) / currentScale;
		int x2 = (pos.x() + scrollArea->horizontalScrollBar()->value()) / currentScale;
		int y2 = (pos.y() + scrollArea->verticalScrollBar()->value()) / currentScale;

		const ftk::Image::Info *info;
		if(channelImg)    info = channelImg->GetImageInfo();	//Get info of new image
		else if(labelImg) info = labelImg->GetImageInfo();
		else return;
		int totalWidth = (*info).numColumns;
		int totalHeight = (*info).numRows;
		int currentZ = vSpin->value();

		if(x2 >= totalWidth) x2 = totalWidth-1;
		if(y2 >= totalHeight) y2 = totalHeight-1;

		if(rubberBand)
		{
			delete rubberBand;
			rubberBand = NULL;
			origin = QPoint();
			emit boxDrawn(x1, y1, x2, y2, currentZ);
		}
		else	// pointsMode or roiMode
		{
			if( abs(x1-x2) < 5 && abs(y1-y2) < 5 )		//I haven't moved far during this click
			{
				if(pointsMode)
				{
					if(origin3.size() == 0)				//This is the first click
					{
						origin3.push_back(x2);
						origin3.push_back(y2);
						origin3.push_back(currentZ);
						this->repaint();
					}
					else								//Must be the second click
					{
						x1 = origin3.at(0);
						y1 = origin3.at(1);
						int z1 = origin3.at(2);
						origin3.clear();
						pointsMode = false;
						emit pointsClicked(x1,y1,z1,x2,y2,currentZ);
					}
				}
				else //roiMode
				{
					ftk::Object::Point last;
					last.t = hSpin->value();
					last.x = x2;
					last.y = y2;
					last.z = currentZ;
					roiPoints.push_back(last);
					if(roiPoints.size() >= 4)
					{
						ftk::Object::Point first = roiPoints.at(0);
						//Check to see if I've closed the region and should end:
						if( abs(first.x-last.x) < 5 && abs(first.y-last.y) < 5 )
						{
							//Erase last point:
							roiPoints.erase(roiPoints.end()-1);
							roiMode = false;
							
							createROIMask();	//This will call repaint
							emit roiDrawn();
							return;
						}
					}
					this->repaint();
				}
			}
		}
	}
}

void LabelImageViewQT::createROIMask()
{
	roiPoints.push_back( roiPoints[0] );
	int numPoints = (int)roiPoints.size();

	//Turn my list of points into a path
	QPainterPath path;
	path.moveTo(roiPoints[0].x, roiPoints[0].y);
	for (int i=0; i < numPoints; i++)
	{
		path.lineTo(roiPoints[i].x, roiPoints[i].y);
	}
	//Draw the path in an image
	QRect rect = this->displayImage.rect();
	QImage img(rect.width(),rect.height(),QImage::Format_Mono);
	img.fill(Qt::black);
	QPainter painter(&img);
	painter.setPen(Qt::white);
	painter.setBrush(Qt::white);
	painter.drawPath(path);

	roiPoints.clear();

	roiImage = img;
	this->SetROIVisible(true);

	//img.save(QString("mask_test.png"));
}

QImage * LabelImageViewQT::GetROIMaskImage()
{ 
	if(this->showROI)
		return &roiImage;
	else
		return NULL;
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
	this->repaint();
}

//****************************************************************************************************
// Reimplement paintEvent to set the displayImage!!
//****************************************************************************************************
void LabelImageViewQT::paintEvent(QPaintEvent * event)
{	
	QWidget::paintEvent(event);

	if(baseImage.height() <= 0 || baseImage.width() <= 0)
		return;

	displayImage = baseImage;
	QPainter painter(&displayImage);
	//if(labelImg)
	//{
		//painter.setCompositionMode(QPainter::CompositionMode_SourceOver); //Default
		//painter.setCompositionMode(QPainter::CompositionMode_Overlay); //Only see overlay when color matches (almost none)
		//painter.setCompositionMode(QPainter::CompositionMode_Plus); // Looked like SourceOver
		//painter.setCompositionMode(QPainter::CompositionMode_Source); // Only shows source (alpha blending didn't seem to work)
		//painter.setCompositionMode(QPainter::CompositionMode_DestinationOver); //Don't see overlay at all
		//painter.setCompositionMode(QPainter::CompositionMode_Clear);  //Everything was white!
		painter.drawImage(0,0,boundsImage);
	//}

	if(pointsMode)
	{
		if(origin3.size() > 0)
		{
			if( origin3.at(2) == vSpin->value() )
			{
				painter.setPen(Qt::red);
				painter.drawPoint( origin3.at(0), origin3.at(1) );
			}
		}
	}
	else if(roiMode)
	{
		painter.setPen(Qt::red);
		for(int p=1; p<(int)roiPoints.size(); ++p)
		{
			painter.drawLine( roiPoints.at(p).x, roiPoints.at(p).y, roiPoints.at(p-1).x, roiPoints.at(p-1).y );
		}
	}

	//Do zooming:
	int oldX = scrollArea->horizontalScrollBar()->value();
	int oldY = scrollArea->verticalScrollBar()->value();

	imageLabel->setPixmap(QPixmap::fromImage(displayImage));
	imageLabel->adjustSize();
	if(currentScale != 1)
	{
		QSize newSize = displayImage.size()*currentScale;
		imageLabel->resize(newSize);
	}

	scrollArea->horizontalScrollBar()->setValue(oldX);
	scrollArea->verticalScrollBar()->setValue(oldY);
}

void LabelImageViewQT::refreshBaseImage()
{
	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();	//Get info of new image
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	int totalWidth = (*info).numColumns;
	int totalHeight = (*info).numRows;

	//baseImage = QImage(totalWidth, totalHeight, QImage::Format_ARGB32);	
	baseImage = QImage(totalWidth, totalHeight, QImage::Format_ARGB32_Premultiplied);	
	baseImage.fill(qRgb(0,0,0));

	if(!channelImg)
		return;

	QPainter painter(&baseImage);
	painter.setCompositionMode(QPainter::CompositionMode_Plus);

	int currentZ = vSpin->value();
	int currentT = hSpin->value();

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
			painter.drawImage(0,0,gray);
		}
	}
	this->repaint();
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
			old_v = img->pixelIndex(c,r);
			if( old_v >= threshold)
			{
				new_v = old_v + offset;
				if(new_v > 255) new_v=255;
				else if(new_v < 0) new_v=0;
				img->setPixel(c,r, new_v);
			}
		}
	}
}

void LabelImageViewQT::AdjustImageIntensity(void)
{
	if(!channelImg)
		return;

	IntensityDialog *dialog = new IntensityDialog(backgroundThreshold, foregroundOffset, this);
	connect(dialog, SIGNAL(valuesChanged(int,int)), this, SLOT(adjustImageIntensity(int,int)));
	dialog->show();
}

void LabelImageViewQT::adjustImageIntensity(int threshold, int offset)
{
	backgroundThreshold = threshold;
	foregroundOffset = offset;
	this->refreshBaseImage();
}

void LabelImageViewQT::initGrayscaleColorTable(void)
{
	grayscaleColorTable.clear();
	for(int i=0; i<256; ++i)
	{
		grayscaleColorTable.append(qRgb(i,i,i));
	}
}

void LabelImageViewQT::selectionChange(void)
{
	this->goToSelection();
	this->refreshBoundsImage();
}

void LabelImageViewQT::goToSelection(void)
{
	if(!centerMap || !bBoxMap) return;

	std::set<long> sels = selection->getSelections();
	if(sels.size() != 1) return;

	long id = *(sels.begin());
	int visZ = vSpin->value();

	if( visZ < ((*bBoxMap)[id]).min.z || visZ > ((*bBoxMap)[id]).max.z )
		vSpin->setValue( ((*centerMap)[id]).z );
}

void LabelImageViewQT::refreshBoundsImage(void)
{
	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();	//Get info of new image
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	//int chs = (*info).numChannels;
	int h = (*info).numRows;
	int w = (*info).numColumns;
	//int currentZ = vSpin->value();
	//int currentT = hSpin->value();

	boundsImage = QImage(w, h, QImage::Format_ARGB32_Premultiplied);	
	boundsImage.fill(qRgba(0,0,0,0));

	QPainter painter(&boundsImage);
	this->drawObjectBoundaries(&painter);
	this->drawObjectIDs(&painter);
	this->drawObjectCentroids(&painter);
	this->drawSelectionCrosshairs(&painter);
	this->drawROI(&painter);
	this->repaint();
}

void LabelImageViewQT::drawObjectBoundaries(QPainter *painter)
{
	if(!showBounds) return;
	if(!labelImg) return;

	const ftk::Image::Info *info = labelImg->GetImageInfo();

	int chs = (*info).numChannels;
	int h = (*info).numRows;
	int w = (*info).numColumns;
	int currentZ = vSpin->value();
	int currentT = hSpin->value();

	for(int ch = 0; ch < chs; ++ch)
	{
		QColor qcolor;
		if(chs > 1)
		{
			std::vector<unsigned char> color = (*info).channelColors[ch];
			qcolor = QColor(color[0],color[1],color[2]);
		}
		else
		{
			qcolor = (*colorItemsMap)["Object Boundaries"];
		}

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
						{
							if(selection->isSelected(v))
							{
								painter->setPen( (*colorItemsMap)["Selected Objects"] );
							}
						}
						painter->drawPoint(j,i);
					}
				}
			}
		}
	}
}

// SHOULD REWRITE THIS SO IT DOESN'T NEED CENTROID, JUST DRAWS ID
void LabelImageViewQT::drawObjectIDs(QPainter *painter)
{
	if(!showIDs) return;
	if(!labelImg) return;
	if(!centerMap) return;
	if(!bBoxMap) return;

	int currentZ = vSpin->value();

	//Iterate through each object and write its id at its centroid.
	std::map<int,ftk::Object::Point>::iterator it;
	for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	{
		int id = (*it).first;
		ftk::Object::Point point = (*it).second;
		//if ( (currentZ == point.z) )
		if( currentZ >= ((*bBoxMap)[id]).min.z && currentZ <= ((*bBoxMap)[id]).max.z )
		{
			painter->setPen( (*colorItemsMap)["Object IDs"] );
			painter->drawText(point.x + 2, point.y - 2, QString::number(id));
		}
	}
}

void LabelImageViewQT::drawObjectCentroids(QPainter *painter)
{
	if(!showCentroids) return;
	if(!labelImg) return;
	if(!centerMap) return;

	int currentZ = vSpin->value();

	//Iterate through each object and write its id at its centroid.
	std::map<int,ftk::Object::Point>::iterator it;
	for ( it = centerMap->begin() ; it != centerMap->end(); ++it )
	{
		int id = (*it).first;
		int cls = 1;
		if(classMap.size() > 0)
			cls = classMap[id];

		int numColors = (int)centroidColorTable.size();
		QColor myColor = centroidColorTable.at( (cls-1)%numColors );
		painter->setPen(Qt::black);
		painter->setBrush(myColor);

		ftk::Object::Point point = (*it).second;
		if ( (currentZ == point.z) )
		{
			//painter->drawRect(point.x - 2, point.y - 2, 5, 5);
			painter->drawEllipse(point.x - 2, point.y - 2, 5, 5);
		}
	}
}

void LabelImageViewQT::drawSelectionCrosshairs(QPainter *painter)
{
	if(!showCrosshairs) return;
	if(!centerMap) return;

	const ftk::Image::Info *info;
	if(channelImg)    info = channelImg->GetImageInfo();	//Get info of image
	else if(labelImg) info = labelImg->GetImageInfo();
	else return;
	int h = (*info).numRows;
	int w = (*info).numColumns;

	std::set<long> sels = selection->getSelections();
	if(sels.size() > 0)
	{
		int id = *(sels.begin());
		ftk::Object::Point center = (*centerMap)[id];
		painter->setPen(Qt::gray);
		painter->drawLine(center.x, 0, center.x, h-1);
		painter->drawLine(0, center.y, w-1, center.y);
	}

}

void LabelImageViewQT::drawROI(QPainter *painter)
{
	if(!showROI) return;

	int h = roiImage.rect().height();
	int w = roiImage.rect().width();
	int v, v1, v2, v3, v4;

	for(int i=1; i < w-1; i++)
	{
		for(int j=1; j < h-1; j++)
		{
			v = roiImage.pixelIndex(i, j);
			if (v > 0)
			{
				v1 = roiImage.pixelIndex(i, j+1);
				v2 = roiImage.pixelIndex(i+1, j);
				v3 = roiImage.pixelIndex(i, j-1);
				v4 = roiImage.pixelIndex(i-1, j);
				if(v!=v1 || v!=v2 || v!=v3 || v!=v4)
				{
					painter->setPen( (*colorItemsMap)["ROI Boundary"] );
					painter->drawPoint(i,j);
				}
			}
		}
	}
}

QVector<QColor> LabelImageViewQT::CreateColorTable()
{
	QVector<QColor> colorTable;
	
	//Colors for 10 classes
	colorTable.append(Qt::magenta);
	colorTable.append(Qt::cyan);
	colorTable.append(Qt::yellow);
	colorTable.append(Qt::red);
	colorTable.append(Qt::green);
	colorTable.append(Qt::blue);
	colorTable.append(Qt::white);
	colorTable.append(Qt::gray);
	colorTable.append(Qt::darkGreen);
	colorTable.append(Qt::darkBlue);

	return colorTable;
}

QString LabelImageViewQT::GetColorNameFromTable( int class_num ){
	QString color_name;
	QVector<QColor> colors;
	colors = CreateColorTable();
	if( colors.at( class_num ) == Qt::yellow )
		color_name = "Yellow";
	else if( colors.at( class_num ) == Qt::white )
		color_name = "White";
	else if( colors.at( class_num ) == Qt::green )
		color_name = "Green";
	else if( colors.at( class_num ) == Qt::cyan )
		color_name = "Cyan";
	else if( colors.at( class_num ) == Qt::magenta )
		color_name = "Magenta";
	else if( colors.at( class_num ) == Qt::red )
		color_name = "Red";
	else if( colors.at( class_num ) == Qt::blue )
		color_name = "Blue";
	else if( colors.at( class_num ) == Qt::gray )
		color_name = "Gray";
	else if( colors.at( class_num ) == Qt::darkGreen )
		color_name = "DarkGreen";
	else if( colors.at( class_num ) == Qt::darkBlue )
		color_name = "DarkBlue";
	else color_name = "No more colors";
	return color_name;
}

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
