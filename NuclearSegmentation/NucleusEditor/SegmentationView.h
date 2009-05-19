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
	//void setZoom(double zf);

signals:
	void goToZ(int newZ);

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
};

#endif 

