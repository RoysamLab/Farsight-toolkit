//****************************************************************************************
// Adapted from Qt Example pieview.cpp
//
// This view inherits from QAbstractItemView and is therefore part of the model/view 
// framework of QT.  It allows for data in a model to be represented by scatterplots.
// The selection method is slightly adjusted so that a plot point is shown as selected
// if any of the items in its row are actually selected.
//****************************************************************************************
#include <math.h>

#include "ScatterView.h"

ScatterView::ScatterView(QWidget *parent) : QAbstractItemView(parent)
{
    horizontalScrollBar()->setRange(0, 0);
    verticalScrollBar()->setRange(0, 0);

	//setSelectionMode(QAbstractItemView::ExtendedSelection);
	setSelectionMode(QAbstractItemView::MultiSelection);
	//setSelectionMode(QAbstractItemView::SingleSelection);

	//Set Defaults
	columnNumForX = 1;
	columnNumForY = 2;
	columnNumForColoring = -1;

	mySettings = new PlotSettings();

	selMode = 0;
    rubberBand = 0;

	x = 0;
	y = 0;
}

//********************************************************************************************
// SLOT
//********************************************************************************************
void ScatterView::SetColForX(int x)
{ 
	if (x > 0 && x < model()->columnCount())
	{
		columnNumForX = x; 
		updateAxis();
		viewport()->update();
	}
}
//********************************************************************************************
// SLOT
//********************************************************************************************
void ScatterView::SetColForY(int y)
{ 
	if (y > 0 && y < model()->columnCount())
	{
		columnNumForY = y; 
		updateAxis();
		viewport()->update();
	}
}
//********************************************************************************************
// SLOT
//********************************************************************************************
void ScatterView::SetColForColor(int c, QMap<int, QColor>  newMap)
{ 
	if (c > 0 && c < model()->columnCount())
	{
		columnNumForColoring = c;
		colorMap = newMap;
		viewport()->update();
	}
}
void ScatterView::SetColForColor(int c)
{
	if (c > 0 && c < model()->columnCount())
	{
		columnNumForColoring = c;
		colorMap = GetDefaultColors();
		viewport()->update();
	}
}

//*******************************************************************************************
// SLOT: We have two selection modes. One is the standard Single Selection mode 
//    Implemented by QT.
//      The second mode allows for the use of the select button to select all points
//		within a region with boundary points defined by selectionRegion.
//**************************************************************************************
void ScatterView::selModeChanged(int s)
{
	if (s == 0)
	{
		//setSelectionMode(QAbstractItemView::SingleSelection);
		setSelectionMode(QAbstractItemView::MultiSelection);
		selectionRegion.clear();
		selMode = 0;
	}
	else
	{
		setSelectionMode(QAbstractItemView::NoSelection);
		selectionRegion.clear();
		selMode = 1;
	}
	viewport()->update();
}

//**************************************************************************************
// SLOT: call this slot to update the selection.  All objects inside the region defined
//  for the points in selectionRegion will be selected.
//**************************************************************************************
void ScatterView::selectClicked(void)
{
	if (selMode != 1)
		return;

	if(selectionRegion.size() <= 0)
		return;

	//I replace the last point with the first point to create a closed loop
	int numPoints = selectionRegion.size();
	selectionRegion[numPoints-1] = selectionRegion[0];

	//I turn my list of clicks into a path
	QPainterPath path;
	path.moveTo(selectionRegion[0]);
	for (int i=1; i < numPoints; i++)
	{
		path.lineTo(selectionRegion[i]);
	}

	//Draw the path in an image
	//QRect rect(Margin, Margin, viewport()->width() - 2*Margin, viewport()->height() - 2*Margin);
	QRect rect(LMargin, TMargin, viewport()->width() - (LMargin+RMargin), viewport()->height() - (BMargin+TMargin));
	QImage img(rect.width(),rect.height(),QImage::Format_Mono);
	img.fill(Qt::white);
	QPainter painter(&img);
	painter.setPen(Qt::black);
	painter.setBrush(Qt::black);
	painter.drawPath(path);

	//Turn the image into a Region
	QBitmap bitmap = QBitmap::fromImage(img, Qt::MonoOnly | Qt::ThresholdDither);
	QRegion region(bitmap);
	//QItemSelectionModel::SelectionFlags command = QItemSelectionModel::Select;
	//setSelection(region,QItemSelectionModel::ClearAndSelect);
	setSelection(region,QItemSelectionModel::Select);
	viewport()->update();
	//selectionRegions.append(selectionRegion);
	selectionRegion.clear();
}

//**************************************************************************************
// All objects in region are selected/de according to selection flags command
//**************************************************************************************
void ScatterView::setSelection(const QRegion &region, QItemSelectionModel::SelectionFlags command)
{
	int rows = model()->rowCount(rootIndex());
	QItemSelection selection;
	selection.clear();
    //QModelIndexList indexes;
	for (int row = 0; row < rows; ++row) 
	{
		//get the region of the item
		QModelIndex index1 = model()->index(row, 0, rootIndex());
		QModelIndex index2 = model()->index(row, model()->columnCount()-1, rootIndex());
        QRegion point = itemRegion(index1);
		//if it intersects with the selection region, save the index
		if (!point.intersect(region).isEmpty())
		{
			//indexes.append(index);
			selection.merge(QItemSelection(index1,index2),command);
		}
    }
	selectionModel()->select(selection, command);

	/*
	if (indexes.size() > 0) 
	{	
		QItemSelection selection;
		selection.clear();
        for (int i = 0; i < indexes.size(); ++i) 
		{
			//Add each item's region to the selection
			selection.merge(QItemSelection(indexes[i],indexes[i]),command);
        }
        selectionModel()->select(selection, command);
		selectionModel()->setCurrentIndex(indexes.last(),QItemSelectionModel::NoUpdate);
    } 
	else 
	{
        QModelIndex noIndex;
        QItemSelection selection(noIndex, noIndex);
        selectionModel()->select(selection, command);
    }
	*/
    update();
}

//**************************************************************************************
// SLOT: clears the selection model and selectionRegion
//**************************************************************************************
void ScatterView::clearClicked(void)
{
	selectionRegion.clear();
	//selectionRegions.clear();
	selectionModel()->clear();
	viewport()->update();
}

//*********************************************************************************************************
//	This slot is called when items are changed in the model. 
//  The changed items are those from topLeft to bottomRight inclusive. 
//  If just one item is changed topLeft == bottomRight. 
//*********************************************************************************************************
void ScatterView::dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight)
{
    QAbstractItemView::dataChanged(topLeft, bottomRight);
	
	//Find max/min x/y and adjust the axis accordingly
	updateAxis();

	viewport()->update();
}

//**********************************************************************************************************
// This is a convenience slot that shows that all data has changed
//**********************************************************************************************************
void ScatterView::dataChanged()
{
	const QModelIndex topLeft = model()->index(0,0);
	const QModelIndex bottomRight = model()->index( model()->rowCount() - 1, model()->columnCount() - 1 );
	this->dataChanged(topLeft,bottomRight);
}

//*****************************************************************************************
// Add to the setModel function so that the axis are updated
//*****************************************************************************************
void ScatterView::setModel(QAbstractItemModel *model)
{
	QAbstractItemView::setModel(model);
	updateAxis();
	viewport()->update();
}

//*********************************************************************************************************
// Find the maximum and minimum values of the data points and update the axis to show all data
//*********************************************************************************************************
void ScatterView::updateAxis(void)
{
	double x1 = 0;
	double x2 = 0;
	double y1 = 0;
	double y2 = 0;
	double v = 0;
	QModelIndex index;
	int rows = model()->rowCount(rootIndex());

	if (rows > 0)
	{	//initialize values
		index = model()->index(0, columnNumForX, rootIndex());
		x1 = model()->data(index).toDouble();
		x2 = x1;
		index = model()->index(0, columnNumForY, rootIndex());
		y1 = model()->data(index).toDouble();
		y2 = y1;
	}
	for (int row = 1; row < rows; ++row) 
	{
		//find min/max of each axis
		index = model()->index(row, columnNumForX, rootIndex());
		v = model()->data(index).toDouble();
		if (v < x1) x1 = v;
		else if (v > x2) x2 = v;
		index = model()->index(row, columnNumForY, rootIndex());
		v = model()->data(index).toDouble();
		if (v < y1) y1 = v;
		else if (v > y2) y2 = v;
    }
	
	//now extend axis slightly so all data points will be visible
	double xrange = x2-x1;
	double yrange = y2-y1;
	double xp = .03*double(xrange);
	double yp = .03*double(yrange);

	mySettings->setRange(x1-xp,x2+xp,y1-yp,y2+yp);
	mySettings->adjust();
}

//*********************************************************************************************************
//	Starts editing the item corresponding to the given index if it is editable
//*********************************************************************************************************
bool ScatterView::edit(const QModelIndex &index, EditTrigger trigger, QEvent *event)
{
   //return QAbstractItemView::edit(index, trigger, event);
	return false;
}

//*********************************************************************************************************
//Returns the item that covers the coordinate given in the view.
//
//Required function
//*********************************************************************************************************
QModelIndex ScatterView::indexAt(const QPoint &point) const
{
	QModelIndex retval = QModelIndex();
	// Transform the view coordinates into contents viewport coordinates.
    double wx = point.x() + horizontalScrollBar()->value();
	double wy = point.y() + verticalScrollBar()->value();

	// Find the extrema of the plot area in viewport coordinates
	//QRect rect(Margin, Margin, viewport()->width() - 2*Margin, viewport()->height() - 2*Margin);
	QRect rect(LMargin, TMargin, viewport()->width() - (LMargin+RMargin), viewport()->height() - (BMargin+TMargin));

	//First check to be sure clicked location is inside of plot area
	if ( (wx > rect.left()) && (wx < rect.right()) && (wy > rect.top()) && (wy < rect.bottom()) )
	{
		//Convert wx and wy to real x and real y values
		PlotSettings settings = *mySettings;

		for (int row = 0; row < model()->rowCount(rootIndex()); ++row) 
		{
			QModelIndex indexX = model()->index(row, columnNumForX, rootIndex());
            QModelIndex indexY = model()->index(row, columnNumForY, rootIndex());
            double valueX = model()->data(indexX).toDouble();
			double valueY = model()->data(indexY).toDouble();

            double dx = valueX - settings.minX;
			double dy = valueY - settings.minY;
			double x = rect.left() + (dx * (rect.width() - 1) / settings.spanX());
			double y = rect.bottom() - (dy * (rect.height() - 1) / settings.spanY());

			if ( ( wx > (x-2) ) && ( wx < (x+2) ) && ( wy > (y-2) ) && ( wy < (y+2) ) ) 
			{
				//returns the index of the first column which is ID
				retval = model()->index(row, 0/*columnNumForX*/, rootIndex());
				break;
			}
		}
	}
	//retval = model()->index(0, 0, rootIndex());
	return retval;
}

//*********************************************************************************************************
//Returns false because we do not hide any items in our view
//
//Required function
//*********************************************************************************************************
bool ScatterView::isIndexHidden(const QModelIndex & /*index*/) const
{
	return false;
}

//*********************************************************************************************************
//     Returns the rectangle of the item with index index in the
//     model. The rectangle is in viewport coordinates.
//*********************************************************************************************************
QRect ScatterView::itemRect(const QModelIndex &index) const
{
	//std::cerr << "In itemRect" << std::endl;

	if (!index.isValid())
		return QRect();

	// Get current viewport region
	//QRect viewportRect(Margin, Margin, viewport()->width() - 2*Margin, viewport()->height() - 2*Margin);
	QRect viewportRect(LMargin, TMargin, viewport()->width() - (LMargin+RMargin), viewport()->height() - (BMargin+TMargin));

	// If not a valid viewport leave
	if (!viewportRect.isValid())
		return QRect();

	//PlotSettings settings = PlotSettings();
	PlotSettings settings = *mySettings;
	int row = index.row();
	QModelIndex indexX = model()->index(row, columnNumForX, rootIndex());
    QModelIndex indexY = model()->index(row, columnNumForY, rootIndex());

    double valueX = model()->data(indexX).toDouble();
	double valueY = model()->data(indexY).toDouble();

    double dx = valueX - settings.minX;
	double dy = valueY - settings.minY;
	double x = viewportRect.left() + (dx * (viewportRect.width() - 1) / settings.spanX());
	double y = viewportRect.bottom() - (dy * (viewportRect.height() - 1) / settings.spanY());

	return QRect(x-3,y-3,7,7);
}


//*********************************************************************************************************
//     Returns the region of the item with index index in the
//     model. The region is in viewport coordinates.
//*********************************************************************************************************
QRegion ScatterView::itemRegion(const QModelIndex &index) const
{
	//std::cerr << "In itemRegion" << std::endl;

    if (!index.isValid())
		 return QRegion();

	
	return itemRect(index);
}
//**************************************************************************************
// If we used scrollbars we would need this
//**************************************************************************************
int ScatterView::horizontalOffset() const
{
    return horizontalScrollBar()->value();
}
//**************************************************************************************
// If we used scrollbars we would need this
//**************************************************************************************
int ScatterView::verticalOffset() const
{
    return verticalScrollBar()->value();
}

void ScatterView::keyPressEvent(QKeyEvent * event)
{
	switch (event->key()) 
	{
		case Qt::Key_Up:
			y--;
			break;
		case Qt::Key_Down:
			y++;
			break;
		case Qt::Key_Left:
			x--;
			break;
		case Qt::Key_Right:
			x++;
			break;
		default:
			QWidget::keyPressEvent(event);
    }
	viewport()->update();
}

void ScatterView::mousePressEvent(QMouseEvent *event)
{
    QAbstractItemView::mousePressEvent(event);
    origin = event->pos();
    /*if (!rubberBand)
	{
        rubberBand = new QRubberBand(QRubberBand::Rectangle, this);
	}
	rubberBand->setGeometry(QRect(origin, QSize()));
	rubberBand->show();*/

}

void ScatterView::mouseMoveEvent(QMouseEvent *event)
{
	QAbstractItemView::mouseMoveEvent(event);
    /*if (rubberBand)
	{
        rubberBand->setGeometry(QRect(origin, event->pos()).normalized());
		//setSelection(QRect(origin, event->pos()).normalized(), QItemSelectionModel::ClearAndSelect);
	}*/
}
//**************************************************************************************
// Upon mouse release, check to make sure I'm close to press (indicating a mouse click)
// Then check my selection mode.  If selection mode is 1 then I add the object under the 
// mouse to the selection
//**************************************************************************************
void ScatterView::mouseReleaseEvent(QMouseEvent *event)
{
    QAbstractItemView::mouseReleaseEvent(event);
	QPoint click = event->pos();

	//Check to make sure this was just a click
	if ( (click.x() < origin.x()-2) || (click.x() > origin.x()+2) )
		return;
	if ( (click.y() < origin.y()-2) || (click.y() > origin.y()+2) )
		return;

	//Check to make sure we are in region selection mode
	if ( selMode != 1 )
		return;

	if ( indexAt(click) == QModelIndex() )
	{
		selectionRegion.append(click);
	}
	if ( indexAt(click) != QModelIndex() )
	{
		setSelection(QRegion(click.x(),click.y(),1,1),QItemSelectionModel::Select);
	}
    /*if (rubberBand)
	{
        rubberBand->hide();
	}	*/
    viewport()->update();
}

//**************************************************************************************
// This event happens when the window is resized.  For simplicity I clear the selection
// region so that I do not have to scale it accordingly and redraw it.
//**************************************************************************************
void ScatterView::resizeEvent(QResizeEvent *event)
{
	selectionRegion.clear();
	//selectionRegions.clear();
	QAbstractItemView::resizeEvent(event);
}

//*****************************************************************************************************************
//	Moves the cursor in the view according to the given cursorAction and keyboard modifiers specified by modifiers.
//*****************************************************************************************************************
QModelIndex ScatterView::moveCursor(QAbstractItemView::CursorAction cursorAction, Qt::KeyboardModifiers /*modifiers*/)
{
     QModelIndex current = currentIndex();
	 //std::cerr << "In moveCursor" << std::endl;
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
     viewport()->update();
     return current;
 }
//**************************************************************************************
//Repaints the view.  Called automatically
//**************************************************************************************
 void ScatterView::paintEvent(QPaintEvent *event)
 {
	 refreshPixmap();
	 QStylePainter painter(viewport());
	 painter.drawPixmap(0, 0, pixmap);
 }

 //*********************************************************************************************************
 //	Scrolls the view if necessary to ensure that the item at index is visible. 
 //	The view will try to position the item according to the given hint
 //*********************************************************************************************************
 void ScatterView::scrollTo(const QModelIndex &index, ScrollHint)
 {
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

     update();
 }

//***********************************************************************************************************************
//   Find the indices corresponding to the extent of the selection.
//	 Applies the selection flags to the items in or touched by the rectangle, rect. 
//   When implementing your own itemview setSelection should call selectionModel()->select(selection, flags) 
//		where selection is either an empty QModelIndex or a QItemSelection that contains all items that are contained in rect. 
//***********************************************************************************************************************
void ScatterView::setSelection(const QRect &rect, QItemSelectionModel::SelectionFlags command)
{	
	//Over-ride command for select only (use clear button to deselect)
	//command = QItemSelectionModel::Select;
	//This is the selection region
    QRect contentsRect = rect.translated( horizontalScrollBar()->value(), verticalScrollBar()->value()).normalized();

	//std::cerr << "left = " << contentsRect.left() << ", right = " << contentsRect.right() << ", top = " << contentsRect.top() << ", bottom = " << contentsRect.bottom() << std::endl;

    int rows = model()->rowCount(rootIndex());

    //QModelIndexList indexes;
	QItemSelection selection;
	selection.clear();
	for (int row = 0; row < rows; ++row) 
	{
		//get the region of the item
		QModelIndex index1 = model()->index(row, 0, rootIndex());
		QModelIndex index2 = model()->index(row, model()->columnCount()-1, rootIndex());
		//QModelIndex index = model()->index(row, 0, rootIndex());
        QRegion region = itemRegion(index1);
		//if it intersects with the selection region, save the index
		if (!region.intersect(contentsRect).isEmpty())
		{
			//indexes.append(index);
			selection.merge(QItemSelection(index1,index2),command);
		}
    }
	selectionModel()->select(selection, command);

	/*
	if (indexes.size() > 0) 
	{	
		QItemSelection selection;
		selection.clear();
        for (int i = 0; i < indexes.size(); ++i) 
		{
			//Add each item's region to the selection
			selection.merge(QItemSelection(indexes[i],indexes[i]),command);
        }
        selectionModel()->select(selection, command);
    } 
	else 
	{
        QModelIndex noIndex;
        QItemSelection selection(noIndex, noIndex);
        selectionModel()->select(selection, command);
    }
	*/
    update();
}

//********************************************************************************************************************
//	 Returns the rectangle on the viewport (in viewport coordinates) occupied by the item at index. 
//   If your item is displayed in several areas then visualRect should return the primary area that contains index 
//		and not the complete area that index might encompasses, touch or cause drawing. 
//   In the base class this is a pure virtual function. 
//********************************************************************************************************************
QRect ScatterView::visualRect(const QModelIndex &index) const
{
	/*  Old Code Saved in case I use scrollbars later
	QRect rect = itemRect(index);
    if (rect.isValid())
        return QRect(rect.left() - horizontalScrollBar()->value(), rect.top() - verticalScrollBar()->value(), rect.width(), rect.height());
    else
        return rect;
	*/
	//std::cerr << "In visualrect" << std::endl;
	return itemRect(index);
	//return viewport()->rect();
}

//********************************************************************************************************************
//     Returns a region corresponding to the selection in viewport coordinates.
//	 Returns the region from the viewport of the items in the given selection.
//********************************************************************************************************************
QRegion ScatterView::visualRegionForSelection(const QItemSelection &selection) const
{
	//std::cerr << "In regionforselection" << std::endl;
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
void ScatterView::refreshPixmap()
{
	pixmap = QPixmap(size());
	//pixmap.fill(this, 0, 0);
	//pixmap.fill(Qt::white);
	pixmap.fill(QColor(220,220,220));
	QPainter painter(&pixmap);
	painter.initFrom(this);
	drawGrid(&painter);
	drawCurves(&painter);
	drawSelection(&painter);
	//painter.drawImage(0,0,bitmap.toImage());
	//update();   //Removed because this calls paintEvent which calls this function
}

//***************************************************************************************************
//The drawGrid() function draws the grid behind the curves and the axes. The area on which we draw
//the grid is specified by rect. If the widget isn't large enough to accommodate the graph, we return
//immediately.
//The first for loop draws the grid's vertical lines and the ticks along the x axis. The second for loop
//draws the grid's horizontal lines and the ticks along the y axis. At the end, we draw a rectangle
//along the margins. The drawText() function is used to draw the numbers corresponding to the tick
//marks on both axes.
//The calls to drawText() have the following syntax:
//painter->drawText(x, y, width, height, alignment, text);
//where (x, y, width, height) define a rectangle, alignment the position of the text within that
//rectangle, and text the text to draw.
//***************************************************************************************************
void ScatterView::drawGrid(QPainter *painter)
{	
	//This rectangle is the graph area
	QRect rect(LMargin, TMargin, viewport()->width() - (LMargin+RMargin), viewport()->height() - (BMargin+TMargin));
	if (!rect.isValid())
		return;
	//Retrieve current plotSettings
	//PlotSettings settings = PlotSettings();
	PlotSettings settings = *mySettings;
	QPen quiteDark = palette().dark().color().light();
	QPen light = palette().light().color();
	//Draw vertical lines
	for (int i = 0; i <= settings.numXTicks; ++i) 
	{
		int x = rect.left() + (i * (rect.width() - 1) / settings.numXTicks);
		double label = settings.minX + (i * settings.spanX() / settings.numXTicks);
		//painter->setPen(quiteDark);
		painter->setPen(QPen(QBrush(Qt::black),1,Qt::DotLine));
		painter->drawLine(x, rect.top(), x, rect.bottom());
		//painter->setPen(light);
		painter->drawLine(x, rect.bottom(), x, rect.bottom() + 5);
		painter->setPen(QPen(QBrush(Qt::black),1));
		painter->drawText(x - 50, rect.bottom() + 5, 100, 15, Qt::AlignHCenter | Qt::AlignTop, QString::number(label));
	}
	for (int j = 0; j <= settings.numYTicks; ++j) 
	{
		int y = rect.bottom() - (j * (rect.height() - 1) / settings.numYTicks);
		double label = settings.minY + (j * settings.spanY() / settings.numYTicks);
		//painter->setPen(quiteDark);
		painter->setPen(QPen(QBrush(Qt::black),1,Qt::DotLine));
		painter->drawLine(rect.left(), y, rect.right(), y);
		//painter->setPen(light);
		painter->drawLine(rect.left() - 5, y, rect.left(), y);
		painter->setPen(QPen(QBrush(Qt::black),1));
		painter->drawText(rect.left() - LMargin, y - 10, LMargin - 5, 20,
			Qt::AlignRight | Qt::AlignVCenter,
			QString::number(label));
	}
	painter->drawRect(rect.adjusted(0, 0, -1, -1));
	
	//Now draw the labels for the x and y axis:
	painter->drawText(rect.left(),rect.bottom() + 20, rect.width(), 20, Qt::AlignHCenter, QString("XXXXXXXXXXXXX"));
	painter->save();
	painter->rotate(-90);
	painter->drawText(-1*rect.bottom(), rect.top() - TMargin, rect.height(), 20, Qt::AlignHCenter, QString("YYYYYYYYYYYYYYYYYYY"));
	painter->restore();
}

//***************************************************************************************************
//The drawCurves() function draws the curves on top of the grid. We start by calling setClipRect() to
//set the QPainter's clip region to the rectangle that contains the curves (excluding the margins and
//the frame around the graph). QPainter will then ignore drawing operations on pixels outside the
//area.
//***************************************************************************************************
void ScatterView::drawCurves(QPainter *painter)
{
	//QItemSelectionModel *selections = selectionModel();

	PlotSettings settings = *mySettings;
	//QRect rect(Margin, Margin, viewport()->width() - 2*Margin, viewport()->height() - 2*Margin);
	QRect rect(LMargin, TMargin, viewport()->width() - (LMargin+RMargin), viewport()->height() - (BMargin+TMargin));
	if (!rect.isValid())
		return;

	painter->setClipRect(rect.adjusted(+1, +1, -1, -1));

	QVector<QRect> rectangles(model()->rowCount(rootIndex()));
	for (int row = 0; row < model()->rowCount(rootIndex()); ++row) 
	{
		QModelIndex indexX = model()->index(row, columnNumForX, rootIndex());
        QModelIndex indexY = model()->index(row, columnNumForY, rootIndex());
		QModelIndex indexC = model()->index(row, columnNumForColoring, rootIndex());
		int c = model()->data(indexC).toInt();
        double valueX = model()->data(indexX).toDouble();
		double valueY = model()->data(indexY).toDouble();

        double dx = valueX - settings.minX;
		double dy = valueY - settings.minY;
		double x = rect.left() + (dx * (rect.width() - 1) / settings.spanX());
		double y = rect.bottom() - (dy * (rect.height() - 1) / settings.spanY());
		rectangles[row] = QRect(x-2,y-2,5,5);

		QColor myColor;
		if (colorMap.size() == 0)
			myColor = Qt::cyan;
		else
			myColor = colorMap.value(c);

		if ( itemInRowIsSelected(row) )
		{
			painter->setPen(myColor);
			painter->setBrush(QBrush(Qt::black,Qt::SolidPattern));
		}
		else
		{
			painter->setBrush(QBrush(myColor,Qt::SolidPattern));
			painter->setPen(Qt::black);
        }
		painter->drawRect(rectangles[row]);
	}
}

//**************************************************************************************
//Create a default colorMap;
//**************************************************************************************
QMap<int, QColor> ScatterView::GetDefaultColors()
{
	QColor defaultColor[6]  = {Qt::yellow, Qt::magenta, Qt::red, Qt::cyan, Qt::blue, Qt::green};

	QMap<int, QColor> classColors;

	//Check to be sure this row contains integers:
	QVariant val = model()->data(model()->index(0,columnNumForColoring));
	std::string type = val.typeName();
	if(type != "int")
	{
		return classColors;
	}
	
	//Now find possible classes(limit: 6 for now)
	std::vector<int> classId(0);
	for (int r=0; r<model()->rowCount(); ++r)
	{
		int clss = model()->data(model()->index(r,columnNumForColoring)).toInt();
		bool found = false;
		for(unsigned int c=0; c<classId.size(); ++c)
		{
			if(classId[c] == clss)
			{
				found = true;
				break;
			}
		}
		if(!found)
		{
			classId.push_back(clss);
			if(classId.size() >= 6)
				break;
		}
	}
	
	//Now assign colors to the classes
	classColors.clear();
	for(unsigned int c=0; c<classId.size(); ++c)
	{
		classColors.insert(classId[c],defaultColor[c]);
	}

	return classColors;
}

//**************************************************************************************
//Checks to see if any of the items in row are selected
//**************************************************************************************
bool ScatterView::itemInRowIsSelected(int row)
{
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

//**************************************************************************************
// This function draws the lines connecting the points in selectionRegion
//**************************************************************************************
void ScatterView::drawSelection(QPainter *painter)
{
	if (selMode == 1)
	{
		painter->setPen(Qt::black);
		for (int i=1; i < selectionRegion.size(); i++)
		{
			painter->drawLine(selectionRegion[i-1],selectionRegion[i]);
		}
	}
	/*for (int i=0; i < selectionRegions.size(); i++)
	{
		for (int j=1; j < selectionRegions[i].size(); j++)
		{
			painter->drawLine(selectionRegions[i][j-1],selectionRegions[i][j]);
		}
	}*/
}


//***************************************************************************************************
//The PlotSettings constructor initializes both axes to the range 0 to 10 with 5 tick marks.
//***************************************************************************************************
PlotSettings::PlotSettings()
{
	minX = 0.0;
	maxX = 10.0;
	numXTicks = 5;
	minY = 0.0;
	maxY = 10.0;
	numYTicks = 5;
}

//***************************************************************************************************
// Allows for the adjustment of minX, maxX, minY, maxY
//***************************************************************************************************
void PlotSettings::setRange(double x1, double x2, double y1, double y2)
{
	minX = x1;
	maxX = x2;
	minY = y1;
	maxY = y2;
}

//***************************************************************************************************
//The adjust() function is called from mouseReleaseEvent() to round the minX, maxX, minY, and maxY
//values to "nice" values and to determine the number of ticks appropriate for each axis. The private
//function adjustAxis() does its work one axis at a time.
//***************************************************************************************************
void PlotSettings::adjust()
{
	adjustAxis(minX, maxX, numXTicks);
	adjustAxis(minY, maxY, numYTicks);
}

//***************************************************************************************************
//The adjustAxis() function converts its min and max parameters into "nice" numbers and sets its
//numTicks parameter to the number of ticks it calculates to be appropriate for the given [min, max]
//range. Because adjustAxis() needs to modify the actual variables (minX, maxX, numXTicks, etc.) and
//not just copies, its parameters are non-const references.
//Most of the code in adjustAxis() simply attempts to determine an appropriate value for the interval
//between two ticks (the "step"). To obtain nice numbers along the axis, we must select the step with
//care. For example, a step value of 3.8 would lead to an axis with multiples of 3.8, which is difficult
//for people to relate to. For axes labeled in decimal notation, "nice" step values are numbers of the
//form 10n, 2·10n, or 5·10n.
//We start by computing the "gross step", a kind of maximum for the step value. Then we find the
//corresponding number of the form 10n that is smaller than or equal to the gross step. We do this by
//taking the decimal logarithm of the gross step, rounding that value down to a whole number, then
//raising 10 to the power of this rounded number. For example, if the gross step is 236, we compute
//log 236 = 2.37291…; then we round it down to 2 and obtain 102 = 100 as the candidate step value
//of the form 10n.
//Once we have the first candidate step value, we can use it to calculate the other two candidates:
//2·10n and 5·10n. For the example above, the two other candidates are 200 and 500. The 500
//candidate is larger than the gross step, so we can't use it. But 200 is smaller than 236, so we use
//200 for the step size in this example.
//It's fairly easy to derive numTicks, min, and max from the step value. The new min value is obtained
//by rounding the original min down to the nearest multiple of the step, and the new max value is
//obtained by rounding up to the nearest multiple of the step. The new numTicks is the number of
//intervals between the rounded min and max values. For example, if min is 240 and max is 1184 upon
//entering the function, the new range becomes [200, 1200], with 5 tick marks.
//This algorithm will give suboptimal resultsin some cases. A more sophisticated algorithm is
//described in Paul S. Heckbert's article "Nice Numbers for Graph Labels" published in Graphics
//Gems (ISBN 0-12-286166-3).
//***************************************************************************************************
void PlotSettings::adjustAxis(double &min, double &max,int &numTicks)
{
	const int MinTicks = 4;
	double grossStep = (max - min) / MinTicks;
	double step = pow(10.0, floor(log10(grossStep)));
	if (5 * step < grossStep) {
		step *= 5;
	} 
	else if (2 * step < grossStep) {
		step *= 2;
	}
	numTicks = int(ceil(max / step) - floor(min / step));
	if (numTicks < MinTicks)
		numTicks = MinTicks;
	min = floor(min / step) * step;
	max = ceil(max / step) * step;
}
