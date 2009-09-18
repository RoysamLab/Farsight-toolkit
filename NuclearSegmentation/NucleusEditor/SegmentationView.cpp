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
// This view is intended to be part of the model view framework so that when items (cells)
// are selected in another view, we can easily locate them in the original data and look 
// at there segmentation.
//
// Currently I use the label image from the segmentation to find out if a click occurred on
// an object, and what object it occurred on.  I then use the bounding box information from
// the model to mark the selections in the window.
//*****************************************************************************************
#include "SegmentationView.h"

SegmentationView::SegmentationView(QWidget *parent) 
  : QAbstractItemView(parent)
{
	resultModel = NULL;
	channelImg = NULL;
	labelImg = NULL;
	showBounds = true;
	showIDs = true;

	currentZ = 0;
	currentT = 0;
	
	currentScale = 1;
	ZoomInFactor = 1.25f;
	ZoomOutFactor = 1 / ZoomInFactor;
	initGrayscaleColorTable();
	backgroundThreshold = 1;		//When adjusting intensities, only change values >= this
	foregroundOffset = 0;			//Offset to ADD to intensity values.

	//numObjects = 0;
	channelFlags.clear();
	channelFlags.push_back(true);

    horizontalScrollBar()->setRange(0,0);	//scrollbars for 2D navigation
    verticalScrollBar()->setRange(0,0);	

	//Added by Aytekin
	horizontalScrollBar()->setValue(0);
	verticalScrollBar()->setValue(0);

	totalWidth = 0;		//width of image in pixels
	totalHeight = 0;	//height of image in pixels

	colorForSelections = Qt::yellow;

	//setSelectionMode(QAbstractItemView::SingleSelection);
	setSelectionMode(QAbstractItemView::MultiSelection);

	setMouseTracking(true);
	rubberBand = NULL;
	//These were just for testing:
	//setToolTip("TOOL TIP");
	//setStatusTip("STATUS TIP");
	//setWhatsThis("WHAT'S THIS");
}

void SegmentationView::setBoundsVisible(bool val)
{
	this->showBounds = val;
	refreshDisplayImage();
}

void SegmentationView::setIDsVisible(bool val)
{
	this->showIDs = val;
	refreshDisplayImage();
}

//*********************************************************************************************************
//	This slot is called when items are changed in the model. 
//  The changed items are those from topLeft to bottomRight inclusive. 
//  If just one item is changed topLeft == bottomRight. 
//*********************************************************************************************************
void SegmentationView::dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
    QAbstractItemView::dataChanged(topLeft, bottomRight);
	refreshDisplayImage();
}

//*****************************************************************************************
// Add to the setModel function so that the mapping is updated
//*****************************************************************************************
void SegmentationView::setModels ( SegmentationModel * sModel ) 
{	
	resultModel = sModel;
	QAbstractItemView::setModel( resultModel->GetModel() );
	QAbstractItemView::setSelectionModel( resultModel->GetSelectionModel() );
	connect(resultModel,SIGNAL(modelChanged()),this, SLOT(refreshDisplayImage()));
	viewport()->update();
}

//*********************************************************************************************************
//	Starts editing the item corresponding to the given index if it is editable
//*********************************************************************************************************
bool SegmentationView::edit(const QModelIndex &index, EditTrigger trigger, QEvent *event)
{
	//return QAbstractItemView::edit(index, trigger, event);
	// Model is not editable from the image view
	return false;
}

//*********************************************************************************************************
//Returns the item that covers the coordinate given in the view.
//
//Required function
//*********************************************************************************************************
QModelIndex SegmentationView::indexAt(const QPoint &point) const
{
	if(!resultModel)
		return QModelIndex();

	QModelIndex retval = QModelIndex();
	// Transform the view coordinates into contents viewport coordinates.
    double wx = ( point.x() + horizontalScrollBar()->value() ) / currentScale;
	double wy = ( point.y() + verticalScrollBar()->value() ) / currentScale;

	//Find the extrema of the image area in viewport coordinates
	//QRect rect(0,0,viewport()->width(),viewport()->height());
	QRect rect(0,0,totalWidth,totalHeight);

	//First check to be sure clicked location is inside of plot area
	if ( (wx > rect.left()) && (wx < rect.right()) && (wy > rect.top()) && (wy < rect.bottom()) )
	{
		if (labelImg)
		{
			//Now find out what Label ID this point has (by checking on label image)
			int labelval = (int)labelImg->GetPixel(currentT,0,currentZ,int(wy),int(wx));
			//return index of first item in model at this row
			if(labelval > 0)
				retval = model()->index(resultModel->RowForID(labelval),0,rootIndex());
		}
		else
		{  /*
			//We have no label image, so must use the seed point marker
			int Xcol = ((FTKItemModel *)model())->getColForSeedX();
			int Ycol = ((FTKItemModel *)model())->getColForSeedY();
			int Zcol = ((FTKItemModel *)model())->getColForSeedZ();

			for (int row = 0; row < model()->rowCount(rootIndex()); ++row) 
			{
				QModelIndex indexX = model()->index(row, Xcol, rootIndex());
				QModelIndex indexY = model()->index(row, Ycol, rootIndex());
				QModelIndex indexZ = model()->index(row, Zcol, rootIndex());
				int x = model()->data(indexX).toInt();
				int y = model()->data(indexY).toInt();
				int z = model()->data(indexZ).toInt();

				if( z == currentZ )
				{
					if ( ( wx > (x-2) ) && ( wx < (x+2) ) && ( wy > (y-2) ) && ( wy < (y+2) ) ) 
					{
						//returns the index of the first column which is ID
						retval = model()->index(row, 0, rootIndex());
						break;
					}
				}
			}*/
		}
	}
	return retval;
}

//*********************************************************************************************************
//Returns false because we do not hide any items in our view
//
//Required function
//*********************************************************************************************************
bool SegmentationView::isIndexHidden(const QModelIndex & /*index*/) const
{
	return false;
}

//*********************************************************************************************************
//     Returns the rectangle of the item with index index in the
//     model. The rectangle is in viewport coordinates.
//*********************************************************************************************************
QRect SegmentationView::itemRect(const QModelIndex &index) const
{
	if (!resultModel)
		return QRect();

	if (!index.isValid())
		return QRect();

	// Get current viewport region
	QRect viewportRect(0, 0, viewport()->width(), viewport()->height());

	// If not a valid viewport leave
	if (!viewportRect.isValid())
		return QRect();

	int row = index.row();	
	QModelIndex indexID = model()->index(row, resultModel->ColumnForID(), rootIndex());
	int id = model()->data(indexID).toInt();
	ftk::Object::Box b = resultModel->SegResult()->GetObjectPtr(id)->GetBoundingBox();

	if ( currentZ >= b.min.z && currentZ <= b.max.z )
	{	
		int minX = int( double(b.min.x)*currentScale );
		int minY = int( double(b.min.y)*currentScale );
		int sizeX = int( double(b.max.x-b.min.x+5)*currentScale );
		int sizeY = int( double(b.max.y-b.min.y+5)*currentScale );
		return QRect( minX, minY, sizeX, sizeY );
	}

	return QRect();
}


//*********************************************************************************************************
//     Returns the region of the item with index index in the
//     model. The region is in viewport coordinates.
//*********************************************************************************************************
QRegion SegmentationView::itemRegion(const QModelIndex &index) const
{
    if (!index.isValid())
		 return QRegion();

	//Need to return a the region in the current z slice for the object with this index!!!
	// HARD SO DO THIS FOR NOW:

	return itemRect(index);
}

void SegmentationView::resizeEvent(QResizeEvent * /* event */)
{
	//Need to adjust the viewport scrollbars
    updateGeometries();
}

//*****************************************************************************
// This funciton sets the range and stepsizes used for the viewport scrollbars
// It is called when the window is resized or zoom is performed
//*****************************************************************************
void SegmentationView::updateGeometries()
{
	//Scaled image Size:
	int w = int( double(totalWidth) * currentScale );		//scaled image width
	int h = int( double(totalHeight) * currentScale );		//scaled image height

	int oldrangeH = horizontalScrollBar()->maximum();
	int oldrangeV = verticalScrollBar()->maximum();

	int oldvalH = horizontalScrollBar()->value();
	int oldvalV = verticalScrollBar()->value();

	double ratioH = 0.5;
	double ratioV = 0.5;

	//if( (oldrangeH == 0) && (oldvalH == 0) )
	//	ratioH = 0.5;
	if(oldrangeH != 0)
		ratioH = (double)oldvalH / (double)oldrangeH;

	//if( (oldrangeV == 0) && (oldvalV == 0) )
	//	ratioV = 0.5;
	if(oldrangeV != 0)
		ratioV = (double)oldvalV / (double)oldrangeV;

	int newrangeH = qMax(0, w - viewport()->width());
	int newrangeV = qMax(0, h - viewport()->height());

	horizontalScrollBar()->setPageStep(viewport()->width());
	horizontalScrollBar()->setRange(0, newrangeH);
	verticalScrollBar()->setPageStep(viewport()->height());
	verticalScrollBar()->setRange(0, newrangeV);

	horizontalScrollBar()->setValue(int(ratioH*newrangeH));
	verticalScrollBar()->setValue(int(ratioV*newrangeV));
}

//*********************************************************************************
// This slot is called automatically when the viewport scrollbars are changed.
//*********************************************************************************
void SegmentationView::scrollContentsBy(int dx, int dy)
{
	viewport()->scroll(dx, dy);
}

int SegmentationView::horizontalOffset() const
{
    return horizontalScrollBar()->value();
}

int SegmentationView::verticalOffset() const
{
    return verticalScrollBar()->value();
}

void SegmentationView::mousePressEvent(QMouseEvent *event)
{	
	QAbstractItemView::mousePressEvent(event);
	QPoint oldOrigin = origin;
	origin = event->pos();				// This is a local position (in viewport coordinates)
	//Compute value in image coordinates and make sure we click within the image:
	int xx = (origin.x() + horizontalOffset()) / currentScale;
	int yy = (origin.y() + verticalOffset()) / currentScale;
	if( xx < 0 || yy < 0 || xx >= totalWidth || yy >= totalHeight )
		return;

	Qt::MouseButton button = event->button();
	if(button == Qt::LeftButton)
	{
		if(!resultModel)  return;
		//if we are in spliting mode, then add the point to the splitting list
		if(resultModel->isSplittingMode())
		{
			resultModel->addPointToSplitList(xx, yy, currentZ);
		}

		if(resultModel->isAddMode())
		{
			int sz = resultModel->getSizeOfAddList();
			if(sz == 0)
			{
				resultModel->addPointToAddList(xx,yy,currentZ);
				if(!rubberBand)	
					rubberBand = new MyRubberBand(this);
				rubberBand->setGeometry(QRect(origin,QSize()));
				rubberBand->show();
			}
			else if(sz == 1)
			{
				if(oldOrigin.x() < origin.x() && oldOrigin.y() < origin.y())
				{
					resultModel->addPointToAddList(xx,yy,currentZ);
				}
				else
				{
					resultModel->abortAdd();
				}
				delete rubberBand;
				rubberBand = NULL;
			}
		}
	}
	else if(button == Qt::RightButton)
	{
		if(rubberBand && resultModel)
		{
			resultModel->abortAdd();
			delete rubberBand;
			rubberBand = NULL;
		}

		if(!this->labelImg)	return;
		int labelval = (int)labelImg->GetPixel(currentT,0,currentZ,int(yy),int(xx));
		//if(labelval > 0)
		QToolTip::showText(event->globalPos(), QString("ID: ") + QString::number(labelval) );	//This shows the tooltip at the global position (screen coordinates)
	}
	
}

void SegmentationView::mouseReleaseEvent(QMouseEvent *event)
{
    QAbstractItemView::mouseReleaseEvent(event);
	/*
	QPoint click = event->pos();
	//Check to make sure this was just a click
	if ( (click.x() < origin.x()-2) || (click.x() > origin.x()+2) )
		return;
	if ( (click.y() < origin.y()-2) || (click.y() > origin.y()+2) )
		return;

	//Check for CTRL modifier
	Qt::KeyboardModifiers modifiers = event->modifiers();
	if ( modifiers.testFlag(Qt::ControlModifier) )
	{
		//QModelIndex index = indexAt(click);
		//selectionModel()->select(QItemSelection(index,index), QItemSelectionModel::SelectCurrent);
		//setSelection(QRect(click.x(),click.y(),1,1), QItemSelectionModel::SelectCurrent);
		//setSelection(QRegion(click.x(),click.y(),1,1),QItemSelectionModel::Select);
	}
	*/
}

void SegmentationView::mouseMoveEvent(QMouseEvent *event)
{
    //QAbstractItemView::mouseMoveEvent(event);
	QPoint pos = event->pos();
	int xx = ( pos.x() + horizontalOffset() ) / currentScale;
	int yy = ( pos.y() + verticalOffset() ) / currentScale;

	//Draw a box!
	if(resultModel)
	{
		if(resultModel->isAddMode())
		{	
			if(rubberBand)
			{
				if(pos.x() <= origin.x() || pos.y() <= origin.y())
					rubberBand->setGeometry(QRect(origin,QSize()));
				else
					rubberBand->setGeometry(QRect(origin, pos).normalized());
			}
		}
	}

	if( xx>=0 && xx<totalWidth && yy>=0 && yy<totalHeight )
		emit mouseAt(xx, yy, currentZ);
}

void SegmentationView::wheelEvent ( QWheelEvent *e )
{
	 int numDegrees = e->delta() / 8;
     double numSteps = numDegrees / 15.0f;
	 double zf = pow(ZoomInFactor, numSteps);
	 zoom(zf);

	if(rubberBand && resultModel)
	{
		resultModel->abortAdd();
		delete rubberBand;
		rubberBand = NULL;
	}
}

void SegmentationView::keyPressEvent(QKeyEvent *event)
 {
     switch (event->key()) {
     case Qt::Key_Plus:
         zoom(ZoomInFactor);
         break;
     case Qt::Key_Minus:
         zoom(ZoomOutFactor);
         break;
     case Qt::Key_Left:
		 horizontalScrollBar()->setValue(horizontalScrollBar()->value() - 2);
         break;
     case Qt::Key_Right:
         horizontalScrollBar()->setValue(horizontalScrollBar()->value() + 2);
         break;
     case Qt::Key_Down:
         verticalScrollBar()->setValue(verticalScrollBar()->value() + 2);
         break;
     case Qt::Key_Up:
         verticalScrollBar()->setValue(verticalScrollBar()->value() - 2);
         break;
	 case Qt::Key_D:	//For delete
		 resultModel->deleteTrigger();
		 break;
	 case Qt::Key_M:
		 resultModel->mergeTrigger();
	 case Qt::Key_S:
		 //resultModel->splitTrigger();
     default:
         QWidget::keyPressEvent(event);
     }
 }

//*****************************************************************************************************************
//	Moves the cursor in the view according to the given cursorAction and keyboard modifiers specified by modifiers.
//*****************************************************************************************************************
QModelIndex SegmentationView::moveCursor(QAbstractItemView::CursorAction cursorAction, Qt::KeyboardModifiers /*modifiers*/)
{
     QModelIndex current = currentIndex();
/*
     switch (cursorAction) {
         case MoveLeft:
         case MoveUp:
             if (current.row() > 0)
                 current = model()->index(current.row() - 1, current.column(),
                                          rootIndex());
             else
                 current = model()->index(0, current.column(), rootIndex());
             break;
         case MoveRight:
         case MoveDown:
             if (current.row() < rows(current) - 1)
                 current = model()->index(current.row() + 1, current.column(),
                                          rootIndex());
             else
                 current = model()->index(rows(current) - 1, current.column(),
                                          rootIndex());
             break;
         default:
             break;
     }
*/
	 //refreshDisplayImage();
     return current;
 }

 //*********************************************************************************************************
 //	Scrolls the view if necessary to ensure that the item at index is visible. 
 //	The view will try to position the item according to the given hint
 //*********************************************************************************************************
void SegmentationView::scrollTo(const QModelIndex &index, ScrollHint)
{
	
	//Also need to change the current Z to make sure visible
	int row = index.row();
	QModelIndex indexID = model()->index(row, resultModel->ColumnForID(), rootIndex());
	int id = model()->data(indexID).toInt();
	
	ftk::Object::Box b = resultModel->SegResult()->GetObjectPtr(id)->GetBoundingBox();
	if( currentZ < b.min.z || currentZ > b.max.z )	//If I'm on wrong z:
	{
		ftk::Object::Point p = resultModel->SegResult()->GetObjectPtr(id)->GetCentroid();
		emit goToZ(p.z);
	}

     QRect area = viewport()->rect();
     QRect rect = visualRect(index);

     if (rect.left() < area.left())
         horizontalScrollBar()->setValue(
             horizontalScrollBar()->value() + rect.left() - area.left());
     else if (rect.right() > area.right())
         horizontalScrollBar()->setValue(
             horizontalScrollBar()->value() + qMin(
                 rect.right() - area.right(), rect.left() - area.left()));

     if (rect.top() < area.top())
         verticalScrollBar()->setValue(
             verticalScrollBar()->value() + rect.top() - area.top());
     else if (rect.bottom() > area.bottom())
         verticalScrollBar()->setValue(
             verticalScrollBar()->value() + qMin(
                 rect.bottom() - area.bottom(), rect.top() - area.top()));

	 //refreshDisplayImage();
	 viewport()->update();
}

//***********************************************************************************************************************
//   Find the indices corresponding to the extent of the selection.
//	 Applies the selection flags to the items in or touched by the rectangle, rect. 
//   When implementing your own itemview setSelection should call selectionModel()->select(selection, flags) 
//		where selection is either an empty QModelIndex or a QItemSelection that contains all items that are contained in rect. 
//***********************************************************************************************************************
 void SegmentationView::setSelection(const QRect &rect, QItemSelectionModel::SelectionFlags command)
 {
	if(!resultModel)
	{
		return;
	}

	if(resultModel->isSplittingMode() || resultModel->isAddMode())
	{
		return;
	}

	//This is the selection region
    QRect contentsRect = rect.translated( horizontalScrollBar()->value(), verticalScrollBar()->value()).normalized();

	QItemSelection selection;
	selection.clear();

	/* REMOVED 9-02-2009 Isaac Abbott
	int rows = model()->rowCount(rootIndex());
	for (int row = 0; row < rows; ++row) 
	{
		//get the region of the item
		QModelIndex index1 = model()->index(row, 0, rootIndex());						
		QModelIndex index2 = model()->index(row,(model()->columnCount())-1, rootIndex());

		//if it intersects with the selection region, save the index
        QRegion region = itemRegion(index1);
		if (!region.intersect(contentsRect).isEmpty())
		{
			//indexes.append(index);
			//Add each item's region to the selection
			selection.merge(QItemSelection(index1,index2),command);
		}
    }
	*/
	// ADDED 9-02-2009 Isaac Abbott
	int x = contentsRect.x();
	int y = contentsRect.y();
	x = x/currentScale;
	y = y/currentScale;

	if (labelImg)
	{
		//Now find out what Label ID this point has (by checking on label image)
		int labelval = (int)labelImg->GetPixel(currentT,0,currentZ,int(y),int(x));
		if(labelval > 0)
		{
			int row = resultModel->RowForID(labelval);
			if(row >= 0)
			{
				QModelIndex index1 = model()->index(row, 0, rootIndex());						
				QModelIndex index2 = model()->index(row,(model()->columnCount())-1, rootIndex());
				selection.merge(QItemSelection(index1,index2),command);
			}
		}
	}

	selectionModel()->select(selection, command);
    viewport()->update(); 
 }

//********************************************************************************************************************
//	 Returns the rectangle on the viewport (in viewport coordinates) occupied by the item at index. 
//   If your item is displayed in several areas then visualRect should return the primary area that contains index 
//		and not the complete area that index might encompasses, touch or cause drawing. 
//   In the base class this is a pure virtual function. 
//********************************************************************************************************************
QRect SegmentationView::visualRect(const QModelIndex &index) const
{
	QRect rect = itemRect(index);
    if (rect.isValid())
        return QRect(rect.left() - horizontalScrollBar()->value(), rect.top() - verticalScrollBar()->value(), rect.width(), rect.height());
    else
        return rect;
}

//********************************************************************************************************************
//   Returns a region corresponding to the selection in viewport coordinates.
//	 Returns the region from the viewport of the items in the given selection.
//********************************************************************************************************************
QRegion SegmentationView::visualRegionForSelection(const QItemSelection &selection) const
{
    int ranges = selection.count();

    if (ranges == 0)
        return QRect();

	QRegion region;
    for (int i = 0; i < ranges; ++i) {
		QItemSelectionRange range = selection.at(i);
        for (int row = range.top(); row <= range.bottom(); ++row) 
		{
            for (int col = range.left(); col <= range.right(); ++col) 
			{
				QModelIndex index = model()->index(row, col, rootIndex());
                region += visualRect(index);
            }
        }
    }
    return region;
} 


//************************************************************************************
// paintEvent updates the offscreen pixmap and then paints the image.
//************************************************************************************
void SegmentationView::paintEvent(QPaintEvent *event)
{
	if (displayImage.width() == 0 || displayImage.height() == 0)
		return;

	//First draw the image and boundaries
	QImage tmpImg(displayImage.width(),displayImage.height(),QImage::Format_ARGB32);
	QPainter p(&tmpImg);
	p.drawImage(0,0,displayImage);
	//Then draw the selection markers
	if(selectionModel())
		drawSelectionMarkers(&p);

	QSize newSize = tmpImg.size()*abs(currentScale);
	
	//Paint the image on the viewport
	QStylePainter painter(viewport());

	//Want to move to center
	//int xshift = qMax(0, ( viewport()->width() - newSize.width() ) / 2 );
	//int yshift = qMax(0, ( viewport()->height() - newSize.height() ) / 2 );

	

	int xshift=0;int yshift=0;
	painter.drawImage(xshift - horizontalScrollBar()->value(), 
		                yshift - verticalScrollBar()->value(),
						tmpImg.scaled(newSize,Qt::KeepAspectRatio,Qt::FastTransformation));

	//horizontalScrollBar()->setValue(0);
	//verticalScrollBar()->setValue(0);
	/*
	int h=horizontalScrollBar()->value();
	int v=verticalScrollBar()->value();
	int *l,*t,*r,*b;
	l=new int();t=new int();r=new int();b=new int();
	*l=0;*t=0;*r=0;*b=0;
	viewport()->getContentsMargins(l,t,r,b);
	int ll,tt,rr,bb;
	ll=*l;
	tt=*t;
	rr=*r;
	bb=*b;
*/
}


 //***************************************************************************************************
//The refreshPixmap() function redraws the plot onto the off-screen pixmap and updates the display.
//We resize the pixmap to have the same size as the widget and fill it with the widget's erase color.
//This color is the "dark" component of the palette, because of the call to setBackgroundRole() in the
//Plotter constructor. If the background is a non-solid brush, QPixmap::fill() needs to know the
//offset in the widget where the pixmap will end up to align the brush pattern correctly. Here, the
//pixmap corresponds to the entire widget, so we specify position (0, 0).
//Then we create a QPainter to draw on the pixmap. The initFrom() call sets the painter's pen,
//background, and font to the same ones as the Plotter widget. Next we call drawGrid() and
//drawCurves() to perform the drawing. At the end, we call update() to schedule a paint event for the
//whole widget. The pixmap is copied to the widget in the paintEvent() function (p. 123).
//***************************************************************************************************
void SegmentationView::refreshDisplayImage()
{
	if (totalWidth == 0 || totalHeight == 0)
	{
		displayImage.fill(0);
		return;
	}

	displayImage = QImage(totalWidth,totalHeight,QImage::Format_ARGB32);
	displayImage.fill(qRgb(0,0,0));
	QPainter painter(&displayImage);
	painter.setCompositionMode(QPainter::CompositionMode_Plus);
	painter.initFrom(this);
	if(channelImg)
	{
		drawImage(&painter);
	}
	painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
	if(labelImg)
	{
		drawBoundaries(&painter);
	}
	if(resultModel)
	{
		drawObjectIDs(&painter);
		drawExclusionMargin(&painter);
	}

	//drawObjects(&painter);

	viewport()->update();
}

//***************************************************************************************
// The drawImage function draws the iamge
//***************************************************************************************
void SegmentationView::drawImage(QPainter *painter)
{
	if(!channelImg)
		return;

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
void SegmentationView::scaleIntensity(QImage *img, int threshold, int offset)
{
	if( threshold < 0) threshold = 0;

	img->setColorTable(grayscaleColorTable);
	int old_v, new_v;
	//int max = 0;
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
				//if(old_v > max) max = old_v;
			}
		}
	}
	//std::cerr << "max: " << max << std::endl;
}

void SegmentationView::AdjustImageIntensity(void)
{
	IntensityDialog *dialog = new IntensityDialog(backgroundThreshold, foregroundOffset, this);
	connect(dialog, SIGNAL(valuesChanged(int,int)), this, SLOT(AdjustImageIntensity(int,int)));
	dialog->show();
}

void SegmentationView::AdjustImageIntensity(int threshold, int offset)
{
	backgroundThreshold = threshold;
	foregroundOffset = offset;
	this->refreshDisplayImage();
}

void SegmentationView::initGrayscaleColorTable(void)
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
void SegmentationView::drawBoundaries(QPainter *painter)
{
	if(!labelImg)
		return;

	if(!showBounds)
		return;

	const ftk::Image::Info *info = labelImg->GetImageInfo();

	int h = (*info).numRows;
	int w = (*info).numColumns;

	QImage borders(w,h,QImage::Format_ARGB32);
	borders.fill(qRgba(0,0,0,0));
	int v, v1, v2, v3, v4;
	for(int i=1; i < h-1; i++)
	{
		for(int j=1; j < w-1; j++)
		{
			v = (int)labelImg->GetPixel(0, 0, currentZ, i, j);
			if (v > 0)
			{
				v1 = (int)labelImg->GetPixel(0, 0, currentZ, i, j+1);
				v2 = (int)labelImg->GetPixel(0, 0, currentZ, i+1, j);
				v3 = (int)labelImg->GetPixel(0, 0, currentZ, i, j-1);
				v4 = (int)labelImg->GetPixel(0, 0, currentZ, i-1, j);
				if(v!=v1 || v!=v2 || v!=v3 || v!=v4)
				{
					QColor myColor = Qt::cyan;
					/*
					if(resultModel)
					{
						int row = resultModel->RowForID(v);
						int columnNumForColor = resultModel->ColumnForColor();
						QModelIndex indexC = model()->index(row, columnNumForColor, rootIndex());
						int c = model()->data(indexC).toInt();
						myColor = resultModel->ColorMap().value(c);
						borders.setPixel(j,i,qRgb(0,255,0));
					}
					else
					{
						myColor = Qt::cyan;
					}
					*/
					borders.setPixel(j,i,myColor.rgb());
				}
				else
				{
					borders.setPixel(j,i,qRgba(0,0,0,0));
				}
			}
		}
	}
	painter->drawImage(0,0,borders);
}

void SegmentationView::drawObjectIDs(QPainter *painter)
{
	if(!resultModel)
		return;

	if(!showIDs)
		return;

	QColor myIdColor = Qt::magenta;
	QColor myBadColor = Qt::red;

	//Iterate through each object and write its id at its centroid.
	std::vector<ftk::Object> * objects = resultModel->SegResult()->GetObjectsPtr();
	for(int i=0; i<(int)objects->size(); ++i)
	{
		ftk::Object * obj = &(objects->at(i));
		int id = obj->GetId();

		//Find z extremes:
		int max_z = 0;
		int min_z = 1000;
		ftk::Object::Box b = obj->GetBoundingBox();
		ftk::Object::Point center = obj->GetCentroid();

		if ( (currentZ >= b.min.z ) && (currentZ <= b.max.z) )
		{
			QFont f = painter->font();
			f.setBold(true);
			painter->setFont(f);
			char val = obj->GetValidity();
			if( val == ftk::Object::VALID )
			{
				painter->setPen(myIdColor);
				painter->drawText(center.x, center.y, QString::number(id));
			}
			else if( val == ftk::Object::EXCLUDED )
			{
				painter->setPen(myBadColor);
				painter->drawText(center.x, center.y, QString("X"));
			}
		}
	}
}

void SegmentationView::drawExclusionMargin(QPainter *painter)
{
	if(!resultModel)
		return;

	if(!channelImg && !labelImg)
		return;

	if(resultModel->xyMargin == 0 && resultModel->zMargin == 0)
		return;

	QColor myMarginColor = Qt::white;

	const ftk::Image::Info *info;
	if(channelImg)
		info = channelImg->GetImageInfo();
	else
		info = labelImg->GetImageInfo();

	int zSlices = (*info).numZSlices;

	int x = 0 + resultModel->xyMargin;
	int y = 0 + resultModel->xyMargin;
	int w = totalWidth - (2*x);
	int h = totalHeight - (2*y);
	int min_z = 0 + resultModel->zMargin;
	int max_z = zSlices - min_z - 1;

	if ( (currentZ >= min_z ) && (currentZ <= max_z) )
	{
		painter->setPen(myMarginColor);
		painter->drawRect(x, y, w, h);
	}
}

//***************************************************************************************
// The drawSelectionMarkers function paints the boundary and the ID in the appropriate color
//***************************************************************************************
void SegmentationView::drawSelectionMarkers(QPainter *painter)
{
	if(!resultModel)
		return;

	QItemSelectionModel *selections = selectionModel();
	QModelIndexList selIndices = selections->selectedIndexes();

	painter->setPen(colorForSelections);
	QFont f = painter->font();
	f.setBold(true);
	painter->setFont(f);

	for (int selIndex = 0; selIndex < selIndices.size(); ++selIndex) 
	{
		int row = selIndices.at(selIndex).row();
		QModelIndex indexID = model()->index(row, resultModel->ColumnForID(), rootIndex());
		int id = model()->data(indexID).toInt();
		ftk::Object::Box b = resultModel->SegResult()->GetObjectPtr(id)->GetBoundingBox();

		if ( (currentZ >= b.min.z ) && (currentZ <= b.max.z) )
		{
			//currently in a slice that should show the box
			/* This code is to draw a box around selected objects REMOVED 9/02/2009 IMA
			painter->drawRect( b.min.x, b.min.y, (b.max.x-b.min.x+1), (b.max.y-b.min.y+1) );
			//*******************************************************************************/
			//Replaced with this code to redraw the boundary of the object in this slice:
			int v, v1, v2, v3, v4;
			for(int r = b.min.y; r <= b.max.y; r++)
			{
				for(int c = b.min.x; c <= b.max.x; c++)
				{
					v = (int)labelImg->GetPixel(0, 0, currentZ, r, c);
					if (v > 0 && v == id)
					{
						v1 = (int)labelImg->GetPixel(0, 0, currentZ, r, c+1);
						v2 = (int)labelImg->GetPixel(0, 0, currentZ, r+1, c);
						v3 = (int)labelImg->GetPixel(0, 0, currentZ, r, c-1);
						v4 = (int)labelImg->GetPixel(0, 0, currentZ, r-1, c);
						if(v!=v1 || v!=v2 || v!=v3 || v!=v4)
						{
							painter->drawPoint(c, r);
						}
					}
				}		
			}
			//*******************************************************************************/
		}
		ftk::Object::Point center = resultModel->SegResult()->GetObjectPtr(id)->GetCentroid();
		if ( (currentZ >= b.min.z ) && (currentZ <= b.max.z) )
		{
			painter->drawText(center.x, center.y, QString::number(id));
		}
	}

	//Added by Yousef 7-31-2009, modified by Isaac 9-02-2009 to only draw points inside objects
	//if we are selecting points for splitting, mark the points
	if(resultModel->isSplittingMode())
	{
		painter->setPen(Qt::red);
		for(int i=0; i<resultModel->getSizeOfSplittingList(); i++)		
		{
			std::vector<int> xx = resultModel->getPointFromSplittingList(i);
			if ( currentZ == xx.at(2) )
				painter->drawPoint(xx.at(0), xx.at(1));
		}
	}
}

void SegmentationView::drawObjects(QPainter *painter)
{/*
	int columnNumForClass = ((FTKItemModel *)model())->getColForColoring();
	int Xcol = ((FTKItemModel *)model())->getColForSeedX();
	int Ycol = ((FTKItemModel *)model())->getColForSeedY();
	int Zcol = ((FTKItemModel *)model())->getColForSeedZ();

	for (int row = 0; row < model()->rowCount(rootIndex()); ++row) 
	{
		QModelIndex indexX = model()->index(row, Xcol, rootIndex());
        QModelIndex indexY = model()->index(row, Ycol, rootIndex());
		QModelIndex indexZ = model()->index(row, Zcol, rootIndex());
		QModelIndex indexC = model()->index(row, columnNumForClass, rootIndex());
		int c = model()->data(indexC).toInt();
		int x = model()->data(indexX).toInt();
        int y = model()->data(indexY).toInt();
		int z = model()->data(indexZ).toInt();

		if( z == currentZ )
		{
			if ( itemInRowIsSelected(row) )
			{
				painter->setPen((((FTKItemModel*)model())->colors)[uint(c-1) % 6]);
				painter->setBrush(QBrush(Qt::black,Qt::SolidPattern));
			}
			else
			{
				painter->setBrush(QBrush((((FTKItemModel*)model())->colors)[uint(c-1) % 6],Qt::SolidPattern));
				painter->setPen(Qt::black);
			}
			painter->drawRect(x-1,y-1,3,3);
		}
	}*/
}
//**************************************************************************************
//Checks to see if any of the items in row are selected
//**************************************************************************************
bool SegmentationView::itemInRowIsSelected(int row)
{
	if( !resultModel )
		return false;

	QItemSelectionModel *selections = selectionModel();
	QModelIndex index;
	int columns = model()->columnCount();
	for (int column=0; column<columns;column++)
	{
		index = model()->index(row, column, rootIndex());
		if (selections->isSelected(index))
			return true;
	}
	return false;
}

//***************************************************************************************
// This is the original image that the segmentation came from
//***************************************************************************************
void SegmentationView::setChannelImage(ftk::Image::Pointer img)
{	
	if(!img)
	{
		channelImg = NULL;
		if(!labelImg)
		{
			totalWidth = 0;
			totalHeight = 0;
		}
		update();
		refreshDisplayImage();
		return;
	}

	const ftk::Image::Info *info = img->GetImageInfo();

	channelImg = img;
	totalWidth = (*info).numColumns;
	totalHeight = (*info).numRows;
	update();

	if(labelImg)
	{
		if(  labelImg->Size() != channelImg->Size() )
		{
			labelImg = NULL;
			//Sizes don't match so put up message
			/*QMessageBox::critical(0, tr("QMessageBox::critical()"),
				tr("Image Sizes do not match \r\n Clearing labelImg"));*/
		}
	}
	refreshDisplayImage();
}

//***************************************************************************************
// This is the label image from the segmentation result
//***************************************************************************************
void SegmentationView::setLabelImage(ftk::Image::Pointer img)
{
	if(!img)
	{
		labelImg = NULL;
		if(!channelImg)
		{
			totalWidth = 0;
			totalHeight = 0;
		}
		update();
		refreshDisplayImage();
		return;
	}

	const ftk::Image::Info *info = img->GetImageInfo();

	labelImg = img;
	totalWidth = (*info).numColumns;
	totalHeight = (*info).numRows;
	update();

	if(channelImg)
	{
		if(  labelImg->Size() != channelImg->Size() )
		{
			channelImg = NULL;
			//Sizes don't match so put up message
			/*QMessageBox::critical(0, tr("QMessageBox::critical()"),
				tr("Image Sizes do not match \r\n Clearing channelImg"));*/
		}
	}
	refreshDisplayImage();
}

//***************************************************************************************
// Change the zslice currently displayed
//***************************************************************************************
void SegmentationView::setZ(int z)
{
	currentZ = z;
	refreshDisplayImage();
}

//***************************************************************************************
// Change the tslice currently displayed
//***************************************************************************************
void SegmentationView::setT(int t)
{
	currentT = t;
	refreshDisplayImage();
}

//*****************************************************************************************
// change the currentScale
//*****************************************************************************************
void SegmentationView::zoom(double zf)
{
	double newScale = currentScale * zf;
	if(zf > 1 && newScale > 5)
		return;
	if(zf < 1 && newScale < 0.20)
		return;

	currentScale = newScale;
	refreshDisplayImage();
	updateGeometries();
}

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