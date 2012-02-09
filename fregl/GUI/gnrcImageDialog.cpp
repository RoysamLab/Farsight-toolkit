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

#include "gnrcimagedialog.h"
#include <QPainter>
#include <QRectF>
#include <QtGui>
#include <QGraphicsScene>
#include <assert.h>
#include <QList>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkStreamingImageFilter.h"

#include <vil3d/vil3d_load.h>
#include <vil3d/vil3d_save.h>
#include <vil3d/vil3d_image_view.h>
#include <vil3d/vil3d_image_resource.h>

#include <fregl/fregl_util.h>

//Constructor
GNRCImageDialog::GNRCImageDialog(QMainWindow *parent)
{
	action = new QAction(this);
	action->setCheckable(true);
	setupUi(this);

	ImageView->setDragMode(QGraphicsView::NoDrag);
	ImageView->setMouseTracking(false);
	setAttribute(Qt::WA_DeleteOnClose); 
	ITextBrowser->setMaximumHeight(150);
	scaling = 1;
}

// Not in use for now. Replaced by loadXML(.)
void GNRCImageDialog::loadImage( const QString& fileName )
{
	//set range for vertical silder (number of zooms
	VSlider->setRange(0,10);
	//by default, point the slider to the first picture element
	VSlider->setValue(10);

	origImg.load(fileName);
	origSize = origImg.size();

	displayItem.setPixmap(QPixmap::fromImage(origImg));
	scene.addItem(&displayItem);
	this->ImageView->setScene(&scene);
	this->ImageView->repaint();

	connect(VSlider, SIGNAL(valueChanged(int)), this, SLOT(updateScene(int)));
	connect(&scene,SIGNAL(mouseClicked()),this,SLOT(doLastClick()));
	connect(&scene,SIGNAL(mouseDragged()),this,SLOT(doLastDrag()));
}

void GNRCImageDialog::loadXML( const QString& qfileName )
{
	//Initialize space_transformer. 
	//
	std::string xml_fileName, fileName_2d;
	xml_fileName = qfileName.toStdString();
	space_transformer_.read_xml( xml_fileName, montage_3d_default_, fileName_2d );

	//set range for vertical silder (number of zooms). By default, point
	//the slider to the first picture element
	//
	VSlider->setRange(0,10);
	VSlider->setValue(10);

	// Set up the image display
	//
	QString filepath = strippedPath(qfileName);
	QString fileName(fileName_2d.c_str());
	montage_path_ = filepath.toStdString();

	origImg.load(filepath+fileName);
	origSize = origImg.size();
	displayItem.setPixmap(QPixmap::fromImage(origImg));
	scene.addItem(&displayItem);
	this->ImageView->setScene(&scene);
	this->ImageView->repaint();

	//Set up the communication with the scene
	//
	connect(VSlider, SIGNAL(valueChanged(int)), this, SLOT(updateScene(int)));
	connect(&scene,SIGNAL(mouseClicked()),this,SLOT(doLastClick()));
	connect(&scene,SIGNAL(mouseDragged()),this,SLOT(doLastDrag()));
}

void GNRCImageDialog::saveROIAs( const QString& in_fileName, const QString &out_fileName  )
{
	std::string image_name_3d = montage_path_;
	if (in_fileName == "") 
		image_name_3d += montage_3d_default_;
	else image_name_3d = in_fileName.toStdString();
	std::cout<<"Extracing image from "<<image_name_3d<<std::endl;

	// ROI
	int ind = -1;
	for (int i = 0; i<scene.myBoxes.size(); i++) {
		if (scene.myBoxes[i] == scene.lastBox) {
			ind  = i;
		}
	}
	if (ind <0) {
		QString str = "\nERROR => No ROI selected! To select, click on the ROI.\n";
		ITextBrowser->append(str);
		return;
	}

	vnl_vector_fixed<float,3> start,size;
	start[0] = scene.lastBox->rect().left()/scene.zooms[ind]; // first index on X
	start[1] = scene.lastBox->rect().top()/scene.zooms[ind]; // first index on Y
	start[2] = 0; // first index on Z

	//size[0] = 430;
	//size[1] = 430;
	size[0] = scene.lastBox->rect().right()/scene.zooms[ind]-start[0]+1 ; // size along X
	size[1] = scene.lastBox->rect().bottom()/scene.zooms[ind]-start[1]+1; // size along Y
	size[2] = space_transformer_.montage_size()[2];

	// Crop the image and write to a file
#if defined(VCL_WIN32) && !defined(__CYGWIN__)
	std::string image_name_pattern = image_name_3d+std::string("\\slice_###.png");
#else
	std::string image_name_pattern = image_name_3d+std::string("/slice_###.png");
#endif

	vil3d_image_resource_sptr image_sptr = vil3d_load_image_resource( image_name_pattern.c_str());
	if (!image_sptr) std::cout<<"Image not loaded!"<<std::endl;

	std::cout<<"cropped region "<<start<<", size "<<size<<std::endl;
	vil3d_image_view<vxl_byte> view = image_sptr->get_view(start[0],size[0],start[1],size[1], 0, size[2]);
	//vil3d_save( view, out_fileName.toStdString().c_str() );

	//convert to itk format to be saved as 3D tiff
	ImageType::Pointer itk_image = fregl_util_convert_vil_to_itk(view);

	typedef itk::ImageFileWriter< ImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( out_fileName.toStdString() );
	writer->SetInput( itk_image );
	writer->Update();

	/*
	typedef itk::Image< unsigned char, 3 > ImageType;
	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( image_name_3d );

	ImageType::IndexType start;
	start[0] = scene.lastBox->rect().top(); // first index on X
	start[1] = scene.lastBox->rect().left(); // first index on Y
	start[2] = 0; // first index on Z
	ImageType::SizeType size;
	size[0] = scene.lastBox->rect().bottom()-start[0]+1; // size along X
	size[1] = scene.lastBox->rect().right()-start[1]+1; // size along Y
	size[2] = space_transformer_.montage_size()[2];

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractFilterType;
	typedef itk::StreamingImageFilter< ImageType, ImageType > StreamingFilterType;

	ExtractFilterType::Pointer efilter = ExtractFilterType::New();
	StreamingFilterType::Pointer sfilter = StreamingFilterType::New();

	sfilter->SetNumberOfStreamDivisions( 30 );
	sfilter->SetInput( reader->GetOutput() );
	efilter->SetInput( sfilter->GetOutput() );
	efilter->SetExtractionRegion( region );
	efilter->Update();

	typedef itk::ImageFileWriter< ImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( out_fileName.toStdString() );
	writer->SetInput( efilter->GetOutput() );
	writer->Update();
	*/
}

void GNRCImageDialog::deleteROI()
{
	// Delete the lastBox
	int ind = -1;
	for (int i = 0; i<scene.myBoxes.size(); i++) {
		if (scene.myBoxes[i] == scene.lastBox) {
			ind  = i;
		}
	}
	if (ind > -1) {
		scene.removeItem(scene.myBoxes[ind]);
		scene.myBoxes.remove(ind);
		scene.zooms.remove(ind);
		delete scene.lastBox;
		scene.lastBox = NULL;
	}
}

void GNRCImageDialog::updateScene(int zoom_power)
{
	int zoom = zoom_power - 10;

	qreal scale; //w.r.t. the original image size
	if (zoom == 0)
		scale = 1;
	else if (zoom < 0)
		scale = 1/qreal(abs(zoom)+1);
	else 
		scale = (zoom+1);

	// Change the image: The image is scaled, and bounding box of the
	// scene redefined
	//
	QSize newSize;
	QImage newImg;
	//scene.removeItem(&displayItem);

	//Scale the image based on location of the slider
	newSize = origSize*scale;
	if ( !zoom )
		newImg = origImg;
	else 
		newImg = origImg.scaled(newSize,Qt::KeepAspectRatio,Qt::FastTransformation);
	displayItem.setPixmap(QPixmap::fromImage(newImg));
	//scene.addItem(&displayItem);
	scene.setSceneRect(0,0,newSize.width(),newSize.height());

	this->ImageView->update();
	this->ImageView->repaint();

	// Change the rectangles
	QVector<QGraphicsRectItem *>::iterator itr = scene.myBoxes.begin();
	QVector<float>::iterator fitr = scene.zooms.begin();
	for (; itr!=scene.myBoxes.end(); itr++, fitr++) {
		(*itr)->scale(scale/scaling, scale/scaling);
		//(*itr)->scale(scale, scale);
		QPen pen = (*itr)->pen();
		pen.setWidth(2*(*fitr)/scale);
		(*itr)->setPen(pen);
	}

	scaling = scale;
	scene.scaling = scale;
}

//extract only the file name
QString GNRCImageDialog::strippedName(const QString &fullFileName)
{
	return QFileInfo(fullFileName).fileName();
}

QString GNRCImageDialog::strippedPath(const QString &fullFileName)
{
	QString filename = QFileInfo(fullFileName).fileName();
	QString fullPath = QFileInfo(fullFileName).filePath();
	return fullPath.remove(filename);
}

void GNRCImageDialog::doLastClick()
{
	qreal x = scene.lastMouseReleaseEvent.x()/scaling;  //x-coordinate
	qreal y = scene.lastMouseReleaseEvent.y()/scaling;  //y-coordinate
	scene.removeItem(scene.lastBox); //the lastBox is just a point, so remove
	delete scene.lastBox;
	scene.lastBox = NULL;

	// Check if the (x,y) are in a box. If so, change lastBox to the one
	// selected and emit the right signal. The lastBox (the one last
	// clicked) is always colored in magenta.
	for (int i = 0; i<scene.myBoxes.size(); i++) {
		QRectF rect = scene.myBoxes[i]->rect();
		if ( x>rect.left()/scene.zooms[i] && x<rect.right()/scene.zooms[i] && 
			y>rect.top()/scene.zooms[i] && y<rect.bottom()/scene.zooms[i] ) {
				scene.lastBox = scene.myBoxes[i];
				QPen pen = scene.lastBox->pen();
				pen.setBrush(Qt::magenta);
				scene.lastBox->setPen(pen);
		}
		else {
			QPen pen = scene.myBoxes[i]->pen();
			pen.setBrush(Qt::yellow);
			scene.myBoxes[i]->setPen(pen);
		}
	}

	// Output the location of the mouse click, and display the mapping
	// with the original images
	QString str;
	str = "Montage location = ( ";
	str += QString::number(x);
	str += " , ";
	str += QString::number(y);
	str += " ) =>";
	ITextBrowser->append(str);

	// look up the position in the montage space to find images that
	// overlap at this point. The original image names and positions are
	// displayed
	float anchor_x = x+space_transformer_.origin()[0];
	float anchor_y = y+space_transformer_.origin()[1];
	float z = space_transformer_.montage_size()[2]/2;

	vnl_vector_fixed< float, 3 > loc(anchor_x,anchor_y,z); 
	vnl_vector_fixed< float, 2 >   xformed_loc;
	std::vector<std::string> images = space_transformer_.image_names();
	for (unsigned int i = 0; i<images.size(); i++) {

		if ( !space_transformer_.in_image_2d(loc, i, xformed_loc) ) {
			continue;
		}    

		QString msg;
		msg = "\t";
		msg += QString(images[i].c_str());
		msg += " : ( ";
		msg += QString::number(xformed_loc[0]);
		msg += " , ";
		msg += QString::number(xformed_loc[1]);
		msg += " )";
		ITextBrowser->append(msg);
	}
}

void GNRCImageDialog::doLastDrag()
{
	qreal x1 = scene.lastMousePressEvent.x()/scaling;
	qreal y1 = scene.lastMousePressEvent.y()/scaling;
	qreal x2 = scene.lastMouseReleaseEvent.x()/scaling;  //x-coordinate
	qreal y2 = scene.lastMouseReleaseEvent.y()/scaling;  //y-coordinate

	QString str;
	str = "Mouse drag from ";
	str += "x = ";
	str += QString::number(x1);
	str += " y = ";
	str += QString::number(y1);
	str += " to x = ";
	str += QString::number(x2);
	str += " y = ";
	str += QString::number(y2);
	str += ", size = (";
	str += QString::number(abs(x2-x1));
	str += ",";
	str += QString::number(abs(y2-y1));
	str += ")";
	ITextBrowser->append(str);
}

//Returns a pointer to the current image
QImage* GNRCImageDialog::currentImage(void)
{
	return &origImg;
}
