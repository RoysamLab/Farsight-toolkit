#ifndef SLICEVIEW5D_H
#define SLICEVIEW5D_H

//qt includes:
#include <QObject>
#include <QtGui>

//vtk includes:
#include <QVTKWidget.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkCamera.h>

//std includes:

//ftk includes:
#include <FTKImage/ftkImage.h>


class SliceView5D: public QWidget
{
	Q_OBJECT;

public:
	SliceView5D(QString filename);
	~SliceView5D();

private slots:
	void setZ(int);

private:
	QVTKWidget * m_imageview;
	QSlider *vSlider;

	typedef vtkSmartPointer<vtkRenderer> RendererPointerType;
	typedef vtkSmartPointer<vtkImageActor> ImageActorPointerType;

	RendererPointerType m_vtkrenderer;
	std::vector<ImageActorPointerType> m_vtkimageactor;

	ftk::Image *img;

	//Helper function
	std::vector<double> RGBtoHSV(std::vector<unsigned char> rgb);

};

#endif 