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

//************************************************************************
// LabelImageViewQT
//
// This is a QWidget used for displaying an image with corresponding label images.
// Both images may be up to 5 dimensions (up to 6 channels).  When drawing,
// the label images will be used to create boundaries around regions with the
// same label.  When the user clicks in the image, the label value in that
// region will be added to the selection.
//
// The view also has a couple of conveniences like adjusting the intensity
// of the images, show/hide boundaries, retrieve box or retrieve points.
//************************************************************************
#ifndef LABELIMAGEVIEWQT
#define LABELIMAGEVIEWQT

#include <QtGui/QFont>
#include <QtGui/QScrollBar>
#include <QtGui/QScrollArea>
#include <QtGui/QWidget>
#include <QtGui/QPainter>
#include <QtGui/QStylePainter>
#include <QtGui/QMouseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtGui/QToolTip>
#include <QtGui/QSpinBox>
#include <QtGui/QLabel>
#include <QtGui/QGridLayout>
#include <QtGui/QCheckBox>
#include <QtGui/QPushButton>
#include <QtGui/QSlider>
#include <QtCore/QRect>
#include <QtCore/QSize>
#include <QtCore/QPoint>
#include <QtCore/QMap>

#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkObject.h>
#include "ObjectSelection.h"

#include <iostream>
#include <map>

class MyRubberBand;
class IntensityDialog;
class LabelGeometry;

class LabelImageViewQT : public QWidget
{
	Q_OBJECT

public:
	LabelImageViewQT(QWidget *parent = 0);
	void SetChannelImage(ftk::Image::Pointer img);
	ftk::Image::Pointer GetChannelImage(){return channelImg;};
	void SetLabelImage(ftk::Image::Pointer img, ObjectSelection * sels = NULL);
	ftk::Image::Pointer GetLabelImage(){return labelImg;};
	void SetCenterMapPointer(std::map<int, ftk::Object::Point> * cMap = NULL);

public slots:
	void SaveDiplayImageToFile();
	void AdjustImageIntensity();
	void SetBoundsVisible(bool val);
	void SetIDsVisible(bool val);
	void ClearGets(void);
	void GetBox(void);
	void Get2Points(void);
	void update(void);
	void goToSelection(void);
	int GetCurrentZ(void){ return vSpin->value(); };
	int GetCurrentT(void){ return hSpin->value(); };

signals:
	void mouseAt(int x, int y, int z);
	void boxDrawn(int x1, int y1, int x2, int y2, int z);
	void pointsClicked(int x1, int y1, int z1, int x2, int y2, int z2);

protected slots:
	void refreshBaseImage(void);
	void refreshBoundsImage(void);
	void drawObjectIDs(QPainter *painter);
	void selectionChange(void);
	void sliderChange(int v);
	void spinChange(int v);
	void adjustImageIntensity(int threshold, int offset);
	//void refreshFeatures(void);
	void updateVSlider(void);
	void updateHSlider(void);
	void updateChFlags(bool b);

protected:
	void moveEvent (QMoveEvent * event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void keyPressEvent( QKeyEvent *event );
	void paintEvent(QPaintEvent *event);

//private: make protected so we can access from inherited class
	void setupUI(void);
	void zoom(double zf);

	void createChannelWidget(void);
	void removeChannelWidget(void);

	void scaleIntensity(QImage *img, int threshold, int offset);

	void initGrayscaleColorTable();
	QVector<QRgb> grayscaleColorTable;
	QColor colorForSelections;
	QColor colorForNormal;
	QColor colorForIDs;

	//UI Widgets:
	QScrollArea *scrollArea;	//Where the image is displayed
	QLabel *imageLabel;			//Contains the displayed image
    QSlider *vSlider;
	QSpinBox *vSpin;
	QLabel *vLabel;
    QSlider *hSlider;
	QSpinBox *hSpin;
	QLabel *hLabel;

	QWidget *channelWidget;			//pointer to a channelWidget for controlling painted channels
	QCheckBox ** chBoxes;			//pointers to all my checkboxes
	int numChBoxes;					//number of checkboxes (number of channels)
	std::vector<bool> channelFlags;	//is channel is visible or not

	QImage displayImage;				//Currently displayed image
	QImage baseImage;					//The intensity image (2D)
	QImage boundsImage;					//Image containing boundaries (2D)

	ftk::Image::Pointer labelImg;
	std::map<int, ftk::Object::Point> *	centerMap;
	//std::map<int, LabelGeometry> labelGeometries;
	ftk::Image::Pointer channelImg;
	ObjectSelection * selection;

	double currentScale;							//Current scaling of the image 
	double ZoomInFactor;							//Constant zoom-in factor
	double ZoomOutFactor;							//Constant zoom-out factor
	int backgroundThreshold;						//When adjusting intensities, only change values bigger than this
	int foregroundOffset;							//Offset to ADD to intensity values.

	bool showBounds;
	bool showIDs;

	//For collecting two points:
	bool pointsMode;
	std::vector<int> origin3;	//a 3D origin for points mode!!

	//For Getting a Box:
	QPoint origin;
	MyRubberBand *rubberBand;
};

class LabelGeometry
{
public:
	float Centroid[3];			// X, Y, Z
	float BoundingBox[6];		// min(X), max(X), min(Y), max(Y), min(Z), max(Z),...]
};

class MyRubberBand : public QWidget
{
	Q_OBJECT;
public:
	MyRubberBand(QWidget * p = 0);
protected:
     void paintEvent(QPaintEvent *event);
	 void mouseMoveEvent(QMouseEvent *event);
};

class IntensityDialog : public QDialog
{
	Q_OBJECT
public:
	IntensityDialog(int threshold, int offset, QWidget *parent = 0);
signals:
	void valuesChanged(int threshold,int offset);
private:
	QSpinBox *thresholdSpin;
	QSpinBox *offsetSpin;
	QPushButton *hideButton;
private slots:
	void changeThreshold(int v);
	void changeOffset(int v);
};

#endif 

