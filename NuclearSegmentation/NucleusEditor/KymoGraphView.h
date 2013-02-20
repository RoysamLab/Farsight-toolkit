#ifndef _TRACKING_KYMO_VIEW_H
#define _TRACKING_KYMO_VIEW_H



//QT Includes
#include <QtGui/QMainWindow>
#include <QtGui/QLabel>
#include <QtGui/QSlider>
#include <QtGui/QStatusBar>
#include <QtGui/QGridLayout>
#include <QtGui/QDockWidget>
#include <QtGui>
#include <QObject>
#include <QtGui/QDialog>

#include "vtkImageViewer2.h"
#include "vtkImageActor.h"
#include "vtkSmartPointer.h"
#include "QVTKWidget.h"
#include "vtkImageData.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageMapToColors.h"
#include "vtkCallbackCommand.h"
#include "vtkPointPicker.h"
#include "vtkCellPicker.h"
#include "vtkRegularPolygonSource.h"
#include "vtkCubeSource.h"
#include "vtkInteractorStyleImage.h"
#include "vtkAppendPolyData.h"
#include "vtkImageShiftScale.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation2D.h"
#include "Trace.h"
#include "vtkCamera.h"
#include "vtkStringArray.h"
#include "vtkLabeledDataMapper.h"
#include "vtkSelectVisiblePoints.h"
#include "vtkActor2D.h"
#include "vtkIntArray.h"
#include "vtkLabelPlacementMapper.h"
#include "vtkPointSetToLabelHierarchy.h"
#include "vtkCubeAxesActor.h"
#include "vtkCubeAxesActor2D.h"
#include "vtkAxisActor2D.h"
#include <vtkTextProperty.h>
#include <vtkPointWidget.h>
#include <vtkImagePlaneWidget.h>
#include <vtkTable.h>
#include <vtkCallbackCommand.h>
#include <vtkActor2DCollection.h>
#include "vtkSphereSource.h"
#include <vtkTIFFWriter.h>
#include <vtkImageBlend.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPNGWriter.h>
#include <vtkSmartPointer.h>
//#include <vtkFreeTypeUtilities.h>
#include <vtkTextProperty.h>
#include <vtkImageCanvasSource2D.h>

//#include "vtkObjectFactory.h"

//stl includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

//standard c++ includes
#include <stdio.h>
#include <stdlib.h>

//ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include "itkImageToVTKImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "ftkLabelImageToFeatures.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageToTextureFeaturesFilter.h"
#include <itkScalarToRGBColormapImageFilter.h>

//VTK includes
#include "vtkPolyData.h"
#include "vtkContourFilter.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include <vtkOpenGLVolumeTextureMapper2D.h>
#include <vtkOpenGLVolumeTextureMapper3D.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolume.h>
#include <vtkImageImport.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include "vtkRendererCollection.h"
#include "vtkAppendPolyData.h"

// Farsight Includes:
#include "ftkGUI/LabelImageViewQT.h"
#include "ftkGUI/ObjectSelection.h"



//Macros
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define CLAMP(a,min,max) (((a)<(min))?(min):(((a)>(max))?(max):(a)))
#define DEBUG1 
#define DEBUG2 printf
#define DEBUG3
#define PROGRESS printf


typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;

typedef itk::Image<short,3> LabelImageType;
typedef itk::Image<short,2> Label2DImageType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;


typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;;

typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageToVTKImageFilter<InputImageType> ConnectorType;
typedef itk::ImageToVTKImageFilter<Input2DImageType> Connector2DType;

typedef itk::RGBPixel<unsigned char> RGBPixelType;
typedef itk::Image<RGBPixelType,2> RGB2DImageType;

typedef itk::ScalarToRGBColormapImageFilter<Input2DImageType, RGB2DImageType> RGB2DFilterType;

typedef itk::ImageToVTKImageFilter<RGB2DImageType> RGBConnector2DType;

typedef itk::RescaleIntensityImageFilter< Input2DImageType, Input2DImageType > RescaleFilterType;

class vtkSlider2DKymoCallbackBrightness : public vtkCommand
{
public:
  static vtkSlider2DKymoCallbackBrightness *New() 
    { return new vtkSlider2DKymoCallbackBrightness; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

	  	vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		opacityTransferFunction->AddPoint(2,0.0);
		//opacityTransferFunction->AddPoint(2+256*(1-value),0.2);
		opacityTransferFunction->AddPoint(50,value);
		volume->GetProperty()->SetScalarOpacity(opacityTransferFunction);
    }
  vtkSlider2DKymoCallbackBrightness() {

  }

vtkVolume* volume;
};

class vtkSlider2DKymoCallbackContrast : public vtkCommand
{
public:
  static vtkSlider2DKymoCallbackContrast *New() 
    { return new vtkSlider2DKymoCallbackContrast; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkSliderWidget *sliderWidget = 
        reinterpret_cast<vtkSliderWidget*>(caller);
      double value = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();

		vtkColorTransferFunction *colorTransferFunction = vtkColorTransferFunction::New();
		colorTransferFunction->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
		colorTransferFunction->AddRGBPoint(255*(1-value),1,1,1);
		volume->GetProperty()->SetColor(colorTransferFunction);
    }
  vtkSlider2DKymoCallbackContrast() {

  }

  vtkVolume* volume;
};


class TrackingKymoView: public QObject
{
	Q_OBJECT
public:
	TrackingKymoView(ftk::Image::Pointer myImg,std::vector< std::vector<ftk::IntrinsicFeatures> > vvfeatures,LabelImageViewQT * imview = NULL, ObjectSelection * sels = NULL )
	{
		myfeatures = vvfeatures;
		my4DImg = myImg;
		ImageView = new LabelImageViewQT();
		ImageView = imview;
		Selection = new ObjectSelection();
		Selection = sels;
		this->singleTracksVisible = false;


		m_mainwindow = new QMainWindow();
		m_mainwindow->setWindowTitle(tr("Tracking: KymoGraph"));
		m_imageview = new QVTKWidget(m_mainwindow);
		m_mainwindow->setCentralWidget(m_imageview);
		
		statusLabel = new QLabel(QObject::tr(" Ready"));
		m_mainwindow->statusBar()->addWidget(statusLabel,1);
		
		m_vtkrenderer = vtkSmartPointer<vtkRenderer>::New();
		m_vtkrenderer->BackingStoreOn();
		m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
		m_imageview->GetRenderWindow()->Render();
	
		GenerateImages();
		AddSliders();
		GenerateTracks();

		m_trackpoly = m_tobj->GetVTKPolyData();
		m_trackmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		m_trackmapper->SetInput(m_trackpoly);
		m_trackactor = vtkSmartPointer<vtkActor>::New();
		m_trackactor->SetMapper(m_trackmapper);
		m_trackactor->GetProperty()->SetLineWidth(3.0);
		m_vtkrenderer->AddActor(m_trackactor);

		TraceBitVector = m_tobj->CollectTraceBits();
		CubeActor = getTrackPoints(TraceBitVector);


		m_vtkrenderer->AddActor(CubeActor);
		currentTrackActor = vtkSmartPointer<vtkActor>::New();

		// create widgets for interaction:
		//this->CreatePointer3D();
		this->CreateSphereActors();
		this->CreateInteractorStyle();
	
		m_mainwindow->show();
		m_imageview->show();
		m_imageview->setMinimumSize(my4DImg->GetImageInfo()->numColumns*1,my4DImg->GetImageInfo()->numRows*1);
		m_vtkrenderer->Render();

		SaveAnimation();
		//this->SaveMovie();
		//connect(ImageView, SIGNAL(emitTimeChanged()), this, SLOT(refreshImages()));
		//connect(Selection, SIGNAL(changed()), this, SLOT(refreshSelection()));


	}

	// Public Member Functions:
	void GetTable4DImage(std::vector< vtkSmartPointer<vtkTable> >  table){tableVector = table;};
	void SetTraceBitClasses(int n);// needs to be called after the previous function
	void UpdateTraceBitClassesView(void);// needs to be called after the previous function

	void GetTFTable(vtkSmartPointer<vtkTable> table){TFTable = table;};
	void SetTrackClasses(int n);// needs to be called after the previous function
	//void UpdateTrackClassesView(void);
	
	void CreateInteractorStyle();

	void GenerateImages();
	void GenerateTracks();

	// Edit Functions:
	static void HandleKeyPress(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);
	static void PickCell(vtkObject* caller, unsigned long event, void* clientdata, void* callerdata);

	~TrackingKymoView()
	{
		delete m_mainwindow;
		delete m_imageview;
		delete statusLabel;
	}

public slots:
		void refreshSelection(void);
		void ToggleLabelVisibility(void);
		void ToggleSelectionMode(void);
		void refreshImages(void);

private:

	vtkSmartPointer<vtkActor> getTrackPoints(std::vector<TraceBit> vec);
	void CreatePointer3D(void);
	void CreatePlaneWidget3D(void);
	void ResetOriginalColors(void);
	void UpdateViewColors(vtkSmartPointer<vtkUnsignedCharArray> ColorArray);
	vtkSmartPointer<vtkUnsignedCharArray> CreateClassColorArray(void);
	void AddLabelToVTKImage(vtkSmartPointer<vtkImageData> labelImage, vtkSmartPointer<vtkImageData> newlabelImage,float bbox[]);


	void SaveAnimation();
	void SaveMovie(void);
	void AddSliders();
	std::vector<std::vector<unsigned char> > GetClassColors();
	vtkSmartPointer<vtkVolume> getOneVTKVolume(vtkSmartPointer<vtkImageData> vtkim, float colors[3]);

	QMainWindow * m_mainwindow;
	QVTKWidget * m_imageview;
	QLabel * statusLabel;

	TraceObject * m_tobj;

	LabelImageViewQT * ImageView;
	ObjectSelection * Selection;

	int currentTime;
	int numClasses;
	int numTClasses;
	bool labelsVisible;
	bool singleTracksVisible;
	ftk::Image::Pointer my4DImg;
	std::vector< std::vector<ftk::IntrinsicFeatures> > myfeatures;
	std::vector< vtkSmartPointer<vtkTable> > tableVector;
	vtkSmartPointer<vtkTable> TFTable;		// table of track features

	vtkSmartPointer<vtkRenderer> m_vtkrenderer;
	vtkSmartPointer<vtkActor> m_selectionactor;
	vtkSmartPointer<vtkImageData> m_vtkim;
	vtkSmartPointer<vtkVolume> m_vtkvolume;
	vtkSmartPointer<vtkActor> m_trackactor;
	vtkSmartPointer<vtkPolyDataMapper> m_trackmapper;
	vtkSmartPointer<vtkPolyData> m_trackpoly;
	vtkSmartPointer<vtkPolyData> point_poly;

	vtkSlider2DKymoCallbackContrast *m_callback_contrast;
	vtkSlider2DKymoCallbackBrightness *m_callback_brightness;



	std::map<int,std::vector<TraceBit> > TraceBitMap; // key value is time, mapped values are the trace bits 
	std::vector<TraceBit> TraceBitVector;
	vtkSmartPointer<vtkActor> cubeact;
	vtkSmartPointer<vtkActor2D> labelActor;
	vtkSmartPointer<vtkActor> CubeActor;
	vtkSmartPointer<vtkActor> currentTrackActor;
	vtkSmartPointer<vtkCubeAxesActor2D> cubeAxesActor;
	vtkSmartPointer<vtkPointWidget> pointerWidget3d;
	vtkSmartPointer<vtkRenderWindowInteractor> Interactor;
	vtkSmartPointer<vtkCallbackCommand> keyPress;

	// Editing Stuff:
	void CreateSphereActors(void);
	void RemoveSphereActors(void);
	void DrawSphere(double pos[]);
	void SetupInteractorStyle(void);
	void Delete(void);
	void ConnectNodes(void);
	bool IsValidNodeConnection(void);
	void RelabelDelete(void);
	void TogglePickMode(void);
	void UpdateLabels(void);
	void DeleteAndRelabelData(void);
	void DeleteDiffTData(std::set<unsigned int> difft_selection);
	void DeleteSameTData(std::set<unsigned int> samet_selection);
	TraceLine * DeleteTlineRecursive(TraceLine * tline, unsigned int vtk_cell_id, int new_id,int maxtime);
	void UpdateTracePolyData(void);

	

	int GetMaxId(void);
	vtkSmartPointer<vtkCellPicker> myCellPicker;
	vtkSmartPointer<vtkCallbackCommand> isPicked;
	std::set<unsigned int> TSelection;			// set of selected tracks
	std::set<unsigned int> NSelection;			// set of selected nodes

	vtkSmartPointer<vtkActor> sphereActor1;
	vtkSmartPointer<vtkActor> sphereActor2;

	std::vector<ObjectSelection::Point> points_from_delete;


};

#endif
