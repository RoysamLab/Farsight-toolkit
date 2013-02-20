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
#include <QtCore/QSettings>

#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkUtils.h>
#include <ftkFeatures/ftkObject.h>
#include "ObjectSelection.h"

#include "vtkTable.h"
#include "vtkSmartPointer.h"
#include <iostream>
#include <list>
#include <set>
#include <map>
#include "math.h"
#include "float.h"

class MyRubberBand;
class IntensityDialog;
class LabelGeometry;

class LabelImageViewQT : public QWidget
{
	Q_OBJECT

public:
	LabelImageViewQT(QMap<QString, QColor> * new_colorItemsMap = NULL, QWidget *parent = 0);
	~LabelImageViewQT();
	void SetChannelImage(ftk::Image::Pointer img);
	ftk::Image::Pointer GetChannelImage(){return channelImg;};
	void SetLabelImage(ftk::Image::Pointer img, ObjectSelection * sels = NULL);
	ftk::Image::Pointer GetLabelImage(){return labelImg;};
	void SetCenterMapPointer(std::map<int, ftk::Object::Point> * cMap = NULL);
	void SetBoundingBoxMapPointer(std::map<int, ftk::Object::Box> * bMap = NULL);

	// 5D-4D image functions:
	void SetCenterMapVectorPointer(std::vector<std::map<int, ftk::Object::Point> > vectorcenterMap);
	void SetBoundingBoxMapVectorPointer(std::vector<std::map<int, ftk::Object::Box> >  vectorboxMap);

	std::vector<std::map<int, ftk::Object::Point> > GetCenterMapVector(void){return centerMapVector;};
	std::vector<std::map  <int, ftk::Object::Box> > GetBoxMapVector(void){return boxMapVector;};


	void SetCenterMapfromVectorPointer(int time =0);
	void SetBoundingBoxMapfromVectorPointer(int time =0);

	void SetClassMap(vtkSmartPointer<vtkTable> table, std::vector<std::string> columns);
	QString GetColorNameFromTable( int class_num );
	void ClearClassMap(void){ classMap1.clear(); classMap2.clear(); classMap3.clear(); refreshBoundsImage();};

	void SetColorItemsMap(QMap<QString, QColor> * new_colorItemsMap){ colorItemsMap = new_colorItemsMap; refreshBoundsImage(); };
	void SetColorMapForCentroids(QVector<QColor> table){ centroidColorTable = table; refreshBoundsImage(); };
	QVector<QColor> CreateColorTable(void);

	std::vector<bool> GetChannelFlags(void){ return channelFlags; };
	void SetChannelFlags(std::vector<bool> ch_fg);

	QImage * GetDisplayImage(){ return &displayImage; };
	QImage * GetROIMaskImage();
	void SetROIMaskImage( QImage img );

	std::map<int, ftk::Object::Point> * GetCenterMapPointer(){ return centerMap; };
	bool GetCrosshairsVisible(){ return showCrosshairs; };
	bool GetBoundsVisible(){ return showBounds; };
	bool GetIDsVisible(){ return showIDs; };
	bool GetCentroidsVisible(){ return showCentroids; };
	bool GetROIVisible(){ return showROI; };
	bool GetKNeighborsVisible(){ return showKNeighbors; };
	bool GetNucAdjVisible(){ return showNucAdj; };
	bool GetCellAdjVisible(){ return showCellAdj; };
	void SetTable(vtkSmartPointer<vtkTable> table){ Table = table; refreshBoundsImage();};
	void SetNucAdjTable(vtkSmartPointer<vtkTable> NucAdjTable){ NucTable = NucAdjTable; refreshBoundsImage();};
	void SetCellAdjTable(vtkSmartPointer<vtkTable> CellAdjTable){ CellTable = CellAdjTable; refreshBoundsImage();};
	void SetKNeighborTable(vtkSmartPointer<vtkTable> kNeighborsTable){ kNeighborTable = kNeighborsTable;};
	void SetRadNeighborTable(vtkSmartPointer<vtkTable> radNeighborsTable){ radNeighborTable = radNeighborsTable;};
	void DoubleClicksOff(void){ enableDoubleClicks = false;};
	void DoubleClicksOn(void){ enableDoubleClicks = true;};
	//5D Image;
	int GetCurrentTimeVal(void);
	void SetCurrentTimeVal(double time);

	void clearSettings();
		
public slots:
	void SaveDisplayImageToFile(QString fileName);
	void SaveCompositeImageToFile(QString fileName);
	void AdjustImageIntensity();
	void SetBoundsVisible(bool val);
	void SetIDsVisible(bool val);
	void SetCentroidsVisible(bool val);
	void SetCrosshairsVisible(bool val);
	void SetKNeighborsVisibleOn(bool flag);
	void SetRadNeighborsVisibleOn(void);
	void SetQueryViewsOff(void);
	void SetNucAdjVisible(bool val);
	void SetCellAdjVisible(bool val);
	void SetROIVisible(bool val);
	void ClearGets(void);
	void GetBox(void);
	void Get2Points(void);
	void GetROI(void);
	void update();
	void goToSelection(void);
	int GetCurrentZ(void){ return vSpin->value(); };
	int GetCurrentT(void){ return hSpin->value(); };
	void SetColorsToDefaults(void);
	void zoomIn(){ zoom( ZoomInFactor ); };
	void zoomOut(){ zoom( ZoomOutFactor ); };
	void getSnapshots();
	QImage getSnapshotforID(int id);
	QImage getSnapshotforID_1(int id);



signals:
	void mouseAt(int x, int y, int z, int t, std::list<int> v);
	void boxDrawn(int x1, int y1, int x2, int y2, int z);
	void pointsClicked(int x1, int y1, int z1, int x2, int y2, int z2);
	void roiDrawn(void);
	void autoMerge(void);
	void emitTimeChanged(void);

protected slots:
	void refreshBaseImage(void);
	void refreshBoundsImage(void);
	void drawObjectIDs(QPainter *painter);
	void drawObjectBoundaries(QPainter *painter);
	void drawObjectCentroids(QPainter *painter);
	void drawSelectionCrosshairs(QPainter *painter);
	void drawNucAdjacency(QPainter *painter);
	void drawCellAdjacency(QPainter *painter);
	void drawKNeighbors(QPainter *painter);
	void drawRadNeighbors(QPainter *painter);
	void drawROI(QPainter *painter);
	void selectionChange(void);
	void selectionTimeChange(void);
	void sliderChange(int v);
	void hspinChange(int v);
	void vspinChange(int v);
	void adjustImageIntensity(int threshold, int offset);
	//void refreshFeatures(void);
	void updateVSlider(void);
	void updateHSlider(void);
	void initChannelFlags(void);
	void createROIMask(void);


protected:
	void moveEvent (QMoveEvent * event);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseDoubleClickEvent(QMouseEvent *event);
	void keyPressEvent( QKeyEvent *event );
	void paintEvent(QPaintEvent *event);

	void writeSettings();
	void readSettings();

	void setupUI(void);
	void zoom(double zf);

	void createChannelWidget(void);
	void removeChannelWidget(void);

	void scaleIntensity(QImage *img, int threshold, int offset);

	float Distance(int x1, int y1, int x2, int y2);
	float perpDist(int x1, int y1, int x2, int y2, int x3, int y3);

	void initGrayscaleColorTable(void);
	QVector<QRgb> grayscaleColorTable;

	QVector<QColor> centroidColorTable;		//Table of colors for centroids
	QMap<QString, QColor> * colorItemsMap;
	
	//UI Widgets:
	QScrollArea *scrollArea;	//Where the image is displayed
	QLabel *imageLabel;			//Contains the displayed image
    QSlider *vSlider;
	QSpinBox *vSpin;
	QLabel *vLabel;
    QSlider *hSlider;
	QSpinBox *hSpin;
	QLabel *hLabel;

	std::vector<bool> channelFlags;	//is channel is visible or not

	QImage displayImage;				//Currently displayed image
	QImage baseImage;					//The intensity image (2D)
	QImage boundsImage;					//Image containing boundaries and other display items(2D)

	ftk::Image::Pointer labelImg;
	std::map<int, ftk::Object::Point> *	centerMap;
	std::map<int, ftk::Object::Point>::iterator it;
	std::map<int, ftk::Object::Box> * bBoxMap;
	// Amin: 4D images
	std::vector<std::map  <int, ftk::Object::Point> >  centerMapVector;
	std::vector<std::map  <int, ftk::Object::Box> > boxMapVector;

	std::map<int, int> classMap1;
	std::map<int, int> classMap2;
	std::map<int, int> classMap3;
	std::map<int, int> classMap4;
	vtkSmartPointer<vtkTable> Table;
	vtkSmartPointer<vtkTable> NucTable;
	vtkSmartPointer<vtkTable> CellTable;
	vtkSmartPointer<vtkTable> kNeighborTable;
	vtkSmartPointer<vtkTable> radNeighborTable;

	ftk::Image::Pointer channelImg;
	ObjectSelection * selection;

	double currentScale;							//Current scaling of the image 
	double ZoomInFactor;							//Constant zoom-in factor
	double ZoomOutFactor;							//Constant zoom-out factor
	int backgroundThreshold;						//When adjusting intensities, only change values bigger than this
	int foregroundOffset;							//Offset to ADD to intensity values.

	//Saved in settings:
	bool showBounds;
	bool showIDs;
	bool showCentroids;
	bool showCrosshairs;
	bool showROI;		//always comes up false
	bool showKNeighbors;
	bool k_mutual;
	bool showRadNeighbors;
	bool knnDone;
	bool showNucAdj;
	bool showCellAdj;
	bool enableDoubleClicks;
	bool saveSettingsOnExit;

	//For collecting two points:
	bool pointsMode;
	bool roiMode;
	std::vector<int> origin3;	//a 3D origin for points mode!!
	std::vector< ftk::Object::Point > roiPoints;
	QImage roiImage;					//Image containing ROI (2D)

	//For Getting a Box:
	QPoint origin;
	MyRubberBand *rubberBand;
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

