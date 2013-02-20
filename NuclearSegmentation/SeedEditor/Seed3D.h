/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */


#ifndef SEED_EDITOR_H
#define SEED_EDITOR_H

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//QT Includes:
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtCore/QFileInfo>
#include "itkObjectFactory.h"
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLightObject.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"
#include "itkImage.h"

#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkProperty.h"

#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkCommand.h"
#include "vtkRendererCollection.h"

#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkConeSource.h"
#include "vtkCallbackCommand.h"
#include <stdio.h>

#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkImageData.h"
#include "vtkLODActor.h"

#include "vtkImageToStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPoints.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolume.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkPlaybackWidget.h"
#include "vtkPlaybackRepresentation.h"
#include "vtkGlyph3D.h"
#include "vtkFloatArray.h"
#include "vtkCellPicker.h"
#include "vtkCallbackCommand.h"
#include "vtkPointData.h"
#include "vtkSphereWidget.h"
#include "vtkImageReslice.h"
#include "vtkVolume.h"
#include "vtkDataSetMapper.h"
#include "math.h"
#include "vtkCoordinate.h"
#include "vtkImplicitPlaneRepresentation.h"
#include "vtkImplicitPlaneWidget2.h"
#include "vtkPlane.h"
#include "vtkCutter.h"
#include "vtkImageToPolyDataFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkGeometryFilter.h"
#include "vtkSeedRepresentation.h"
#include "vtkOutlineFilter.h"
#include "vtkSeedWidget.h"
#include "vtkImageReslice.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkImageActor.h"
#include "vtkMatrix4x4.h"
#include "vtkImageMapper.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkSphereRepresentation.h"
#include "vtkHandleWidget.h"
#include "vtkInteractorStyleRubberBand2D.h"
#include "vtkPointHandleRepresentation2D.h"
#include "vtkPointHandleRepresentation3D.h"
#include "vtkPlane.h"
#include "vtkProp3DCollection.h"
#include <QObject>
#include <QtGui>
#include <QVTKWidget.h>
#include <iostream>
#include <fstream>


typedef struct
{
float x;
float y;
float z;
} point;

//enum modeStates {noState=0, addState, delState, cpState, fpState};
// when creating a new checkbox, add its index (before BOXMAX_IDX) 
enum StateIndex { idxAdd, idxDel, idxCP, idxFP, idxReset, idxUndoDel, STATE_MAX_IDX};
const short numofStates     = STATE_MAX_IDX;
const float MAX_SEARCH_DIST = 20.0;

// when creating a new glyph instance, add its index (before GLYPHMAX_IDX) and color below
enum GlyphIndex { idxdelAddedglyph, idxdelOrigglyph, idxAddglyph,  idxOrigglyph,  idxCPglyph,    idxFPglyph, idxGTglyph, GLYPHMAX_IDX};
const float colorArray[GLYPHMAX_IDX][3] =
				{{0.5,0.5,0.0},     {0.5,0.2,0.2},  {0.0,0.0,0.9}, {0.9,0.7,0.5}, {0.7,0.0,0.6}, {0.8,0.4,0.0}, {0.35,0.9, 0.9} };
const int NAME_STR_MAX = 64;
const char GLYPH_NAMES[GLYPHMAX_IDX][NAME_STR_MAX] = 
{
	"Deleted From Added list",
		"Deleted From Original list",
		"Added Seeds (False Negatives)",
		"Original remaining Seeds (True Positives)",
		"Cluster Positives",
		"False Positives",
		"Ground Truth"
};

class Seed3D;

 class myGlyph
 {
 public:
	 
	 myGlyph(Seed3D *Owner, const float *color, vtkSmartPointer<vtkSphereSource> sphere);
	 ~myGlyph();
	 vtkSmartPointer<vtkGlyph3D>    Glyph;
	 //void AddPoint(vtkIdType idx, float *p);
	 void AddPoint(float *p);
	 void RemovePoint(vtkIdType idx);
	 vtkIdType Search(double *pickPos, float &mindist, float *finpt);
	 bool isEmpty() {return ( points->GetNumberOfPoints()==0 );}
	 // this is used only for Origglyph, and future GroundTruthglyph
	 void ReadFilePts(std::string sfilename);
	 void SaveToOutfile(std::ofstream *outfstrm );
	 int  GetNumberOfPoints() {return points->GetNumberOfPoints();}
	 void GetPoint( vtkIdType i, float* p) {pcoords->GetTupleValue(i, p);}
	 void Clear();
 private:
	 
	 vtkSmartPointer<vtkFloatArray> pcoords;
	 vtkSmartPointer<vtkPoints> points;
	 vtkSmartPointer<vtkPolyData> polydata;
	 vtkSmartPointer<vtkSphereSource> Sphere;
	 vtkSmartPointer<vtkPolyDataMapper> SphereMapper;
	 vtkSmartPointer<vtkLODActor> SphereActor;
	Seed3D *OwnerClass;
};


class Seed3D : public QMainWindow
{
    Q_OBJECT;
public:
	QVTKWidget *QVTK;
	vtkSmartPointer<vtkPointPicker> PointPicker;
	vtkSmartPointer<vtkRenderer> Renderer;
	Seed3D(QWidget * parent = 0, Qt::WindowFlags flags = 0);
	private slots:
	void loadImage(void);
	void PlaceSeed();
    void ApplyAdd();
	void ApplyDel();

	void AddSeedCheck();
	
	//void Check();
	void loadSeedsGT();
	void saveResult();
	std::string SeedFileDialog(const char *seedtype);

	void PointTransfer(double *pickPos, 
						myGlyph *src1, myGlyph *src2, 
						myGlyph *target1, myGlyph *target2);
    	
private:

	//TaggedGlyph *CPglyph, *FPglyph;
	//OrigGlyph   *Origglyph;
	//DelGlyph    *delAddedglyph, *delOrigglyph;
	//AddedGlyph  *Addglyph;
	
	//this is an array of the myGlyphs needed. They will be initialized 
	// in a loop.
	myGlyph *myGlyphArray[GLYPHMAX_IDX]; 
	// These are the same as above; used as aliases for better readability.
	myGlyph *CPglyph, *FPglyph, *Origglyph, *delAddedglyph, *delOrigglyph, *Addglyph, *GTglyph;

	void createMenus();
	void createStatusBar();
	void CleanUp(); //was DeleteObjects();
	
	//void StateChange(QCheckBox *currBox);
	
	void SetButtonColor(QRadioButton *b, GlyphIndex i);
	
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor1;
    vtkSmartPointer<vtkRenderWindowInteractor> Interactor2;
	vtkSliderRepresentation2D *sliderRep;
	vtkSliderRepresentation2D *sliderRep2;
	vtkSliderRepresentation2D *sliderRep3;
	vtkSliderWidget *sliderWidget;
	vtkSliderWidget *sliderWidget2;
	vtkSliderWidget *sliderWidget3;
    //vtkSmartPointer<vtkCellPicker> CellPicker;
	//vtkSmartPointer<vtkPointPicker> PointPicker;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	//vtkSmartPointer<vtkCallbackCommand> keyPress;

    
	QMenu *fileMenu;
	QAction *loadAction;
	QAction *loadActionGT;
	QAction *saveAction;
	QAction *exitAction;
	QLabel *statusLabel;
	QString lastPath;
	QString fileName;
	
	QWidget *browse;
    
    QVTKWidget *QVTK1;	
    QVTKWidget *QVTK2;

  
	//Qt widgets on the main window
    /*QCheckBox *AddBox;
    QCheckBox *DeleteBox;
    QCheckBox *UndoDelBox;	
    QCheckBox *ResetBox;	
	QCheckBox *CPBox;	
	QCheckBox *FPBox;	*/
	//QCheckBox *FNBox;	
	QRadioButton *AddBtn;
    QRadioButton *DelBtn;
    QRadioButton *UndoDelBtn;	
    QRadioButton *ResetBtn;	
	QRadioButton *CPBtn;	
	QRadioButton *FPBtn;	
	//bool States[numofStates];
	//QRadioButton *BoxArray[numofCheckBoxes];
	QPushButton *PlaceButton;
    QPushButton *ApplyAddButton;
	QPushButton *ApplyDelButton;

	
    vtkSmartPointer<vtkRenderer> Renderer1;
    vtkSmartPointer<vtkRenderer> Renderer2;

	//modeStates mode;
	//int flag; //rename 
	int iRender; 

	//int newseedcounter; 


	vtkSmartPointer<vtkImageData> VTKim;
	//vtkSmartPointer<vtkImageReslice> Reslice;
	//vtkSmartPointer<vtkImageReslice> Reslice1;
	vtkSmartPointer<vtkDataSetMapper>im_mapper;
	vtkSmartPointer<vtkDataSetMapper>im_mapper1;
	vtkSmartPointer<vtkActor> imActor;
	vtkSmartPointer<vtkActor> imActor1;

	//vtkSmartPointer<vtkHandleWidget> widget;
	vtkSmartPointer<vtkPointHandleRepresentation3D> handle;
	//vtkSmartPointer<vtkHandleWidget> widget1;
	//vtkSmartPointer<vtkPointHandleRepresentation2D> handle1;
	//vtkSmartPointer<vtkHandleWidget> widget2;
	//vtkSmartPointer<vtkPointHandleRepresentation2D> handle2;
	vtkSmartPointer<vtkAppendPolyData> apd;
	//vtkSmartPointer<vtkPolyDataMapper> VolumeMapper;
	vtkSmartPointer<vtkActor> VolumeActor;
	vtkSmartPointer<vtkVolume> Volume;
	vtkSmartPointer<vtkPolyData> poly;
	vtkSmartPointer<vtkPolyDataMapper> polymap;

    static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
    void rayCast(char*,char*);
    //std::vector<point> ReadPoints(char* );
	double* y;
	double* wp;
 };


 //TaggedGlyph::Populate()
 //{

 //}

 //class TaggedGlyph : public ParentGlyph
 //{
 //public:
	// TaggedGlyph();
	// virtual void Populate();
	// virtual void Undo();
 //}


 // class OrigGlyph : public ParentGlyph
 //{
 //public:
	// OrigGlyph(char *);
	// virtual void Populate();
	// virtual void Undo();

 //}

 // class AddedGlyph : public ParentGlyph
 //{
 //public:
	// AddedGlyph();
	// virtual void Populate();
	// virtual void Undo();
 //}

 // class DelGlyph : public ParentGlyph
 //{
 //public:
	// DelGlyph();
	// virtual void Populate();
	// virtual void Undo();
 //}


#endif
