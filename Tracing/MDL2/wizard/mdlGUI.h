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

#ifndef MDLGUI_H
#define MDLGUI_H

#include <QtGui/QMainWindow>
#include <QtGui/QLabel>
#include <QtGui/QHBoxLayout>
#include <QtGui/QVBoxLayout>
#include <QtGui/QPushButton>
#include <QtGui/QPixmap>
#include <QtGui/QLineEdit>
#include <QtGui/QFileDialog>
#include <QtGui/QCloseEvent>
#include <QtGui/QGroupBox>

#include "mdlTypes.h"
#include "mdlIntegratedSkeleton.h"
#include "mdlMST.h"
#include "mdlUtils.h"
#include "ftkGUI/PreprocessDialog.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToVTKImageFilter.h"

#include <QVTKWidget.h>
#include "vtkColorTransferFunction.h"
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"
#include "vtkImageData.h"
#include "vtkImageActor.h"
#include "vtkOpenGLVolumeTextureMapper3D.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRendererCollection.h"
#include "vtkSmartPointer.h"
#include "vtkVolumeProperty.h"

class mdlGUI : public QMainWindow
{
    Q_OBJECT;

public:
	mdlGUI(QWidget * parent = 0);
	~mdlGUI();

	void RenderImage(mdl::ImageType::Pointer image);
	vtkSmartPointer<vtkVolumeProperty> NewRGBVolumeProperty(const double range[]);

	void RenderPoints(std::vector<mdl::fPoint3D> points, double r=0.0, double g=0.0, double b=0.0);
	vtkSmartPointer<vtkActor> CreateActorFromPoints(std::vector<mdl::fPoint3D> points);

	void RenderPolyData(std::string filename, double r=0.0, double g=0.0, double b=0.0);
	vtkSmartPointer<vtkActor> CreateActorFromPolyDataFile(const char *filename);

signals:

protected:
	void closeEvent(QCloseEvent *event);

public slots:

private slots:
	void loadImage();
	void preprocess();
	void integratedSkeleton();
	void mstBB();
	void mstSpine();
	void saveAll();
	void showRenderWindow();
    
private:
	QString lastPath;

	QPixmap bannerPixmap;
	QLabel * bannerLabel;

	QHBoxLayout * loadLayout;
	QLabel * loadLabel;
	QLabel * imageLabel;
	QPushButton * loadButton;

	QHBoxLayout * prepLayout;
	QPushButton * prepButton;

	QGroupBox * skelBox;
	QHBoxLayout * skelLayout;
	QLabel * vectMagLabel;
	QLineEdit * vectMagEdit;
	QCheckBox * skelCheckbox;
	QPushButton * skelButton;

	QGroupBox * bbBox;
	QHBoxLayout * bbLayout;
	QLabel * edgeLabel;
	QLineEdit * edgeRangeEdit;
	QLabel * morphLabel;
	QLineEdit * morphStrengthEdit;
	QPushButton * bbButton;

	QGroupBox * spineBox;
	QHBoxLayout * spineLayout;
	QLabel * pruneLabel;
	QLineEdit * pruneEdit;
	QLabel * alphaLabel;
	QLineEdit * alphaEdit;
	QPushButton * spineButton;

	QHBoxLayout * saveLayout;
	QPushButton * saveButton;

	QHBoxLayout * showRenderLayout;
	QPushButton * showRenderButton;

	QVBoxLayout * stepsLayout;
	
	QHBoxLayout * masterLayout;
	QWidget * masterWidget;

	//*****

	mdl::ImageType::Pointer OrigImage;

	QVTKWidget * RenderWidget;
	vtkSmartPointer<vtkRenderer> Renderer;

	typedef itk::ImageToVTKImageFilter< mdl::ImageType > ITKtoVTKType;
	ITKtoVTKType::Pointer ITKtoVTK;

	mdl::ImageType::Pointer PrepImage;

	mdl::IntegratedSkeleton *Skel;
	std::vector<mdl::fPoint3D> SkeletonPoints;

	mdl::MST *MinSpanTree;
	std::vector<mdl::fPoint3D> Nodes;
	std::vector<mdl::pairE> BackbonePairs;
	std::vector<mdl::pairE> SpinePairs;

};


#endif
