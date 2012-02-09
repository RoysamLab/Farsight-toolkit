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

#ifndef _GNRCIMAGEDIALOG_H
#define _GNRCIMAGEDIALOG_H

#include <QMainWindow>
#include <QtGui>
#include <QDialog>
#include <QPaintEvent>
#include <QMouseEvent>
#include <QString> 
#include <QVector>
#include <QPixmap>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsPixmapItem>
#include <QString>
#include <QColor>
#include <QPainter>
#include <QRectF>
#include <QGraphicsScene>
#include <assert.h>
#include "ui_ImageViewArea.h"

#include <string>

// basic header files for reading & writing images
#include "deluxescene.h"

#include <fregl/fregl_space_transformer.h>

class GNRCImageDialog : public QMainWindow, public Ui::IViewForm
{
	Q_OBJECT;

public:
	GNRCImageDialog( QMainWindow *parent = 0);

	~GNRCImageDialog(){};

	void     loadImage( const QString &fileName);
	void     loadXML( const QString &qfileName);
	void     saveROIAs( const QString &in_fileName, const QString &out_fileName );
	void     deleteROI();

	QAction *windowMenuAction() const { return action; }

	QString  strippedName(const QString &fullFileName);
	QString  strippedPath(const QString &fullFileName);

	QImage* currentImage(void);

private:
	QAction *action;
	QString  fileFilters;
	QImage origImg;
	QSize origSize;
	qreal scaling;

	DeluxeScene scene;
	QGraphicsPixmapItem displayItem;
	fregl_space_transformer< InputPixelType > space_transformer_;
	std::string montage_3d_default_;                                          
	std::string montage_path_;

	private slots:
		//void     updateImg(int j);
		void     updateScene(int zoom_power);
		void	   doLastClick();
		void	   doLastDrag();
};

#endif
