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

#ifndef IMAGEBROWSER5D_H
#define IMAGEBROWSER5D_H

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
#include <vtkInteractorStyle.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkVolumeTextureMapper2D.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolume.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>

//std includes:

//ftk includes:
#include <ftkImage/ftkImage.h>


class ImageBrowser5D: public QWidget
{
	Q_OBJECT;

public:
	typedef enum { SLICE, VOLUME } RenderMode; 
	ImageBrowser5D(QString filename, RenderMode mode = SLICE);
	ImageBrowser5D(ftk::Image::Pointer img, RenderMode mode = SLICE);
	~ImageBrowser5D();

	void ToggleMode();
	void SetMode(RenderMode mode);
	void ToggleShowChannel(int channel);
	void SetShowChannel(int channel, bool show);

protected:
	//void keyPressEvent( QKeyEvent *event );

private slots:
	void SetZ(int);
	void SetT(int);

private:
	void Setup();
	void CreateObjects();
	void CreateLayout();
	void CreateLookupTable();
	void CreateVolumeProperties();
	void CreateInteractorStyle();
	void UpdateRenderer();
	void Rerender();
	void UpdateImageActors();
	void UpdateImageVolumes();
	void UpdateVSlider();
	void UpdateHSlider();
	void CleanProps();
	std::vector<double> RGBtoHSV(std::vector<unsigned char> rgb);
	static void keyPress(vtkObject * object, unsigned long eid, void *clientdata, void * callerdata);

	ftk::Image::Pointer img;		//smart pointer cleans itself

	QVTKWidget * m_imageview;			//delete this on close
	QSlider * vSlider;					//all children, delete themselves:
	QSpinBox * vSpin;
	QLabel * vLabel;
    QSlider * hSlider;
	QSpinBox * hSpin;
	QLabel * hLabel;

	typedef vtkSmartPointer<vtkRenderer> RendererPointerType;
	RendererPointerType m_vtkrenderer;

	//2D containers:
	typedef vtkSmartPointer<vtkLookupTable> LookupTablePointerType;
	std::vector<LookupTablePointerType> m_lookuptable;		//One table for each channel in image
	typedef vtkSmartPointer<vtkImageActor> ImageActorPointerType;
	std::vector<ImageActorPointerType> m_channelActors;		//One actor for each channel shown

	//3D containers:
	typedef vtkSmartPointer<vtkVolumeProperty> VolumePropertyPointerType;
	std::vector<VolumePropertyPointerType> m_volumeproperty;//One property for each channel in image
	typedef vtkSmartPointer<vtkVolume> VolumePointerType;
	std::vector<VolumePointerType> m_volumes;				//One volume for each channel shown

	vtkSmartPointer<vtkCallbackCommand> m_kycallback;

	RenderMode m_mode;
	int m_T;	//current time slice (0...)
	std::vector<bool> m_chflag;
};

#endif 
