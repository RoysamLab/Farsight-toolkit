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

#ifndef SEGMENTATIONVIEW_H
#define SEGMENTATIONVIEW_H

#include <QtGui/QAbstractItemView>
#include <QtGui/QFont>
#include <QtGui/QScrollBar>
#include <QtGui/QItemSelection>
#include <QtGui/QItemSelectionModel>
#include <QtGui/QWidget>
#include <QtGui/QPainter>
#include <QtGui/QStylePainter>
#include <QtGui/QMouseEvent>
#include <QtGui/QMessageBox>
#include <QtCore/QModelIndex>
#include <QtCore/QRect>
#include <QtCore/QSize>
#include <QtCore/QPoint>
#include <QtCore/QMap>

#include <ftkImage/ftkImage.h>
#include "SegmentationModel.h"

#include <iostream>

class SegmentationView : public QAbstractItemView
{
	Q_OBJECT

public:
	SegmentationView(QWidget *parent = 0);
	void setChannelImage(ftk::Image::Pointer img);
	ftk::Image::Pointer getChannelImage(){return channelImg;};
	void setLabelImage(ftk::Image::Pointer img);
	ftk::Image::Pointer getLabelImage(){return labelImg;}; 
	void setModels(SegmentationModel *sModel);
    QRect visualRect(const QModelIndex &index) const;
    void scrollTo(const QModelIndex &index, ScrollHint hint = EnsureVisible);
    QModelIndex indexAt(const QPoint &point) const;

public slots:
	void setZ(int z);
	void setT(int t);
	void setChannelFlags(std::vector<bool> flags){channelFlags = flags;refreshDisplayImage();};
	void setBoundsVisible(bool val);
	//void setZoom(double zf);

signals:
	void goToZ(int newZ);
	void mouseAt(int x, int y, int z);

protected slots:
    void dataChanged(const QModelIndex &topLeft, const QModelIndex &bottomRight);
	void refreshDisplayImage(void);
    //void rowsInserted(const QModelIndex &parent, int start, int end);
    //void rowsAboutToBeRemoved(const QModelIndex &parent, int start, int end);

protected:
    bool edit(const QModelIndex &index, EditTrigger trigger, QEvent *event);
    QModelIndex moveCursor(QAbstractItemView::CursorAction cursorAction, Qt::KeyboardModifiers modifiers);

    int horizontalOffset() const;
    int verticalOffset() const;

    bool isIndexHidden(const QModelIndex &index) const;

    void setSelection(const QRect&, QItemSelectionModel::SelectionFlags command);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent ( QWheelEvent *e );
	void keyPressEvent( QKeyEvent *event );

    void paintEvent(QPaintEvent *event);
    void resizeEvent(QResizeEvent *event);
    void scrollContentsBy(int dx, int dy);

    QRegion visualRegionForSelection(const QItemSelection &selection) const;

private:
    QRect itemRect(const QModelIndex &item) const;
    QRegion itemRegion(const QModelIndex &index) const;
    void updateGeometries();
	bool itemInRowIsSelected(int row);
	void updateMapping();
	void zoom(double zf);

	void drawImage(QPainter *painter);
	void drawBoundaries(QPainter *painter);
	void drawSelectionMarkers(QPainter *painter);
	void drawObjects(QPainter *painter);

	QImage displayImage;	//Everything that is being displayed in the viewport (includes selections)
	SegmentationModel *resultModel;
	ftk::Image::Pointer labelImg;
	ftk::Image::Pointer channelImg;
	std::vector<bool> channelFlags;
	int currentZ;
	int currentT;
	double currentScale;							//Current scaling of the image 
	double ZoomInFactor;							//Constant zoom-in factor
	double ZoomOutFactor;							//Constant zoom-out factor

	int totalWidth;									//These are the sizes of the original image (scale = 1)
	int totalHeight;

	QPoint origin;

	QColor colorForSelections;
	bool showBounds;
};

#endif 

