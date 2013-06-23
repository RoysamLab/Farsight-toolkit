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

#include "ImageBrowser5D.h"

//#include "itkImageRegionConstIterator.h"

//Destructor
ImageBrowser5D::~ImageBrowser5D()
{
	m_channelActors.clear();
	m_lookuptable.clear();
	m_volumeproperty.clear();
	m_volumes.clear();
	m_chflag.clear();
	if(m_imageview) delete m_imageview;
}

//Constructors
ImageBrowser5D::ImageBrowser5D(QString filename, RenderMode mode)
{
	img = ftk::Image::New();
	if( !img->LoadFile( filename.toStdString() ) )	return;
	m_mode = mode;
	//if( !img->LoadFileSeries("C:/TestImages/Ying_image - 5D - tiff/021805m5bwt_t%02d.tif",1,10,1) ) return;
	//m_mode = VOLUME;
	this->Setup();
}

ImageBrowser5D::ImageBrowser5D(ftk::Image::Pointer img, RenderMode mode)
{
	this->img = img;
	m_mode = mode;
	this->Setup();
}

void ImageBrowser5D::Setup()
{
	if(!img) return;

	for(int c=0; c<img->GetImageInfo()->numChannels; ++c)
		m_chflag.push_back(true);
	
	CreateObjects();
	CreateLayout();
	CreateLookupTable();
	CreateVolumeProperties();
	CreateInteractorStyle();

	UpdateVSlider();
	UpdateHSlider();
	UpdateRenderer();

	this->resize(img->GetImageInfo()->numColumns, img->GetImageInfo()->numRows);
	this->setAttribute ( Qt::WA_DeleteOnClose );
	this->setWindowTitle(tr("Image Browser"));
}

void ImageBrowser5D::ToggleMode()
{
	if(m_mode == VOLUME)
		SetMode(SLICE);
	else
		SetMode(VOLUME);
}

void ImageBrowser5D::SetMode(RenderMode mode)
{
	if(mode != m_mode)
	{
		m_mode = mode;
		UpdateVSlider();
		UpdateRenderer();
	}
}

void ImageBrowser5D::ToggleShowChannel(int channel)
{
	if(channel >= img->GetImageInfo()->numChannels)
		return;

	if(m_chflag.at(channel) == false)
		SetShowChannel(channel, true);
	else
		SetShowChannel(channel, false);

}
void ImageBrowser5D::SetShowChannel(int channel, bool show)
{
	if(channel >= img->GetImageInfo()->numChannels)
		return;

	if(m_chflag.at(channel) != show)
	{
		m_chflag.at(channel) = show;
		if(m_mode == VOLUME)
			m_volumes.at(channel)->SetVisibility(show);
		else
			m_channelActors.at(channel)->SetVisibility(show);

		Rerender();
	}
}

//Instantiate the objects:
void ImageBrowser5D::CreateObjects()
{
	vSlider = new QSlider(this);
    vSlider->setOrientation(Qt::Vertical);
	vSpin = new QSpinBox(this);
	vSpin->resize( vSpin->minimumSizeHint() );
	vLabel = new QLabel("z",this);

	hSlider = new QSlider(this);
    hSlider->setOrientation(Qt::Horizontal);
	hSpin = new QSpinBox(this);
	hSpin->resize( hSpin->minimumSizeHint() );
	hLabel = new QLabel("t",this);

	//CREATE A QVTK Widget for putting a QVTK renderer in QT window
	m_imageview = new QVTKWidget(this);
	m_vtkrenderer = RendererPointerType::New();
	m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
	//m_vtkrenderer->SetBackground(0.0,0.0,0.0);

	//setup connections
	connect(vSlider, SIGNAL(valueChanged(int)), this, SLOT(SetZ(int)));
	connect(hSlider, SIGNAL(valueChanged(int)), this, SLOT(SetT(int)));
	connect(vSlider, SIGNAL(valueChanged(int)), vSpin, SLOT(setValue(int)));
	connect(hSlider, SIGNAL(valueChanged(int)), hSpin, SLOT(setValue(int)));
	connect(vSpin, SIGNAL(valueChanged(int)), vSlider, SLOT(setValue(int)));
	connect(hSpin, SIGNAL(valueChanged(int)), hSlider, SLOT(setValue(int)));

}
//Create the Layout of all items in this widget.
//Also creates many of the objects!!
void ImageBrowser5D::CreateLayout()
{
	QGridLayout *vsliderLayout = new QGridLayout;
	vsliderLayout->addWidget(vSpin,0,0,1,2);
	vsliderLayout->addWidget(vSlider,1,0,1,1);
	vsliderLayout->addWidget(vLabel,1,1,1,1);

	QHBoxLayout *hsliderLayout = new QHBoxLayout;
	hsliderLayout->addWidget(hSpin);
	hsliderLayout->addWidget(hLabel);
	hsliderLayout->addWidget(hSlider);

	QGridLayout *viewerLayout = new QGridLayout();
    viewerLayout->addWidget(m_imageview, 0, 0);
    viewerLayout->addLayout(vsliderLayout, 0, 1);
	viewerLayout->addLayout(hsliderLayout, 1, 0);

	this->setLayout(viewerLayout);
}

//Creates lookup tables for each channel
void ImageBrowser5D::CreateLookupTable()
{
	if(!img) return;

	m_lookuptable.clear();
	int chs = img->GetImageInfo()->numChannels;
	for(int i=0; i<chs; ++i)
	{
		LookupTablePointerType table = vtkSmartPointer<vtkLookupTable>::New();
		table->SetRange(0, 255);				// image intensity range
		
		std::vector<double> hsv = RGBtoHSV( img->GetImageInfo()->channelColors.at(i) );
		//S and V should always be 1, it is H that I am most concerned about

		table->SetHueRange(hsv[0],hsv[0]);			// Split Range
		table->SetSaturationRange(hsv[1],hsv[1]);	// Full Color Saturation
		table->SetValueRange(hsv[2],hsv[2]);		// from black to white
		table->SetAlphaRange(0,1);
		table->Build();

		m_lookuptable.push_back(table);
	}
}

void ImageBrowser5D::CreateVolumeProperties()
{
	if(!img) return;

	m_volumeproperty.clear();
	int chs = img->GetImageInfo()->numChannels;
	for(int i=0; i<chs; ++i)
	{
		vtkSmartPointer<vtkPiecewiseFunction> opacityTransferFunction = vtkSmartPointer<vtkPiecewiseFunction>::New();
		opacityTransferFunction->AddPoint(2,0.0);
		opacityTransferFunction->AddPoint(255,1.0);

		vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
		colorTransferFunction->AddRGBPoint(0.0,0.0,0.0,0.0);
		std::vector<unsigned char> rgb = img->GetImageInfo()->channelColors.at(i);
		colorTransferFunction->AddRGBPoint(255.0, double(rgb[0])/255.0 , double(rgb[1])/255.0, double(rgb[2])/255.0);

		VolumePropertyPointerType volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
		volumeProperty->SetScalarOpacity(opacityTransferFunction);
		volumeProperty->SetColor(colorTransferFunction);
		volumeProperty->SetInterpolationTypeToLinear();
		volumeProperty->DisableGradientOpacityOn();
		volumeProperty->ShadeOff();

		m_volumeproperty.push_back(volumeProperty);
	}
}

//Sets up the renderer and adds actors to it for each channel.
//This should only be called when T or CH selections change!!
void ImageBrowser5D::UpdateRenderer()
{
	if( !img ) return;

	if(m_mode == SLICE)
		UpdateImageActors();
	else
		UpdateImageVolumes();

	Rerender();
}

void ImageBrowser5D::Rerender()
{
	m_vtkrenderer->ResetCameraClippingRange();	//This fixes the "disappearing image" problem
	m_imageview->GetRenderWindow()->Render();
}

//This functions removes all props from the renderer (2D and 3D) and clears the vector
void ImageBrowser5D::CleanProps(void)
{
	//Clean up existing actors
	for(int i=0; i< (int)m_channelActors.size(); ++i)
	{
		m_vtkrenderer->RemoveActor( m_channelActors.at(i) );
	}
	m_channelActors.clear();

	//Clean up existing volumes
	for(int i=0; i< (int)m_volumes.size(); ++i)
	{
		m_vtkrenderer->RemoveVolume( m_volumes.at(i) );
	}
	m_volumes.clear();
}
//vtkImageActors will display slices of a 3D image.
void ImageBrowser5D::UpdateImageActors(void)
{
	CleanProps();

	//Create new actors:
	int chs = img->GetImageInfo()->numChannels;
	for(int i=0; i<chs; ++i)
	{
		vtkSmartPointer<vtkImageMapToColors> color = vtkSmartPointer<vtkImageMapToColors>::New();
		vtkSmartPointer<vtkImageData> channel = img->GetVtkPtr(m_T,i);
		color->SetInputData( channel );
		color->SetLookupTable( m_lookuptable.at(i) );

		//Create the actor and set its input to the image
		ImageActorPointerType actor = ImageActorPointerType::New();
		actor->SetInputData( color->GetOutput() );
		actor->SetDisplayExtent(0, img->GetImageInfo()->numColumns - 1, 0, img->GetImageInfo()->numRows - 1, 0, 0 );
		actor->SetZSlice( this->vSlider->value() );
		actor->SetVisibility( m_chflag.at(i) );
		actor->RotateWXYZ(180,0,0,1);
		actor->RotateWXYZ(180,0,1,0);
		m_channelActors.push_back(actor);
		m_vtkrenderer->AddActor(actor);
	}
}

//This creates volumes of the 3D image
void ImageBrowser5D::UpdateImageVolumes(void)
{
	CleanProps();

	int chs = img->GetImageInfo()->numChannels;
	for(int i=0; i<chs; ++i)
	{
		/*
		vtkSmartPointer<vtkVolumeTextureMapper3D> vMapper = vtkSmartPointer<vtkVolumeTextureMapper3D>::New();
		vMapper->SetSampleDistance(1);
		vMapper->SetInput(img->GetVtkPtr(m_T,i));
		*/
		/*
		vtkSmartPointer<vtkVolumeRayCastMapper> vMapper = vtkSmartPointer<vtkVolumeRayCastMapper>::New();
		vtkSmartPointer<vtkVolumeRayCastCompositeFunction> compFunction = vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();
		vMapper->SetVolumeRayCastFunction(compFunction);
		vMapper->SetInput(img->GetVtkPtr(m_T,i));
		*/
		vtkSmartPointer<vtkFixedPointVolumeRayCastMapper> vMapper = vtkSmartPointer<vtkFixedPointVolumeRayCastMapper>::New();
		vMapper->SetInputData(img->GetVtkPtr(m_T,i));
		/*
		vtkSmartPointer<vtkVolumeTextureMapper2D> vMapper = vtkSmartPointer<vtkVolumeTextureMapper2D>::New();
		vMapper->SetMaximumNumberOfPlanes(50);
		vMapper->SetInput(img->GetVtkPtr(m_T,i));
		vMapper->Update();
		*/

		vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
		volume->SetMapper(vMapper);
		volume->SetProperty( m_volumeproperty.at(i) );
		volume->SetVisibility( m_chflag.at(i) );
		volume->Update();

		m_volumes.push_back(volume);
		m_vtkrenderer->AddVolume(volume);
	}
}

void ImageBrowser5D::SetZ(int z)
{
	for (int i=0; i<(int)m_channelActors.size(); ++i)
	{
		m_channelActors.at(i)->SetZSlice(z);
	}
	Rerender();
}


void ImageBrowser5D::SetT(int t)
{
	m_T = t;
	UpdateRenderer();
}

//**************************************************************************************************
// 2 functions to update the z and t sliders and their corresponding spin boxes and labels
// according to the current image size.
//**************************************************************************************************
void ImageBrowser5D::UpdateVSlider(void)	//This slider is for z
{
	if(!img || !vSlider || !vSpin || !vLabel) return;
	int z = img->GetImageInfo()->numZSlices;
	vSlider->setRange(0,z-1);
	vSpin->setRange(0,z-1);
	vSlider->setValue(0);
	vSpin->setValue(0);
	if (z > 1 && m_mode == SLICE)
	{
		vSlider->setEnabled(true);
		vSpin->setEnabled(true);
		vLabel->setEnabled(true);
	}
	else
	{
		vSlider->setEnabled(false);
		vSpin->setEnabled(false);
		vLabel->setEnabled(false);
	}
}

void ImageBrowser5D::UpdateHSlider(void)
{
	if(!img || !hSlider || !hSpin || !hLabel) return;
	int t = img->GetImageInfo()->numTSlices;
	hSlider->setRange(0,t-1);
	hSlider->setValue(0);
	hSpin->setRange(0,t-1);
	hSpin->setValue(0);
	if (t > 1)
	{
		hSlider->setEnabled(true);
		hSpin->setEnabled(true);
		hLabel->setEnabled(true);
	}
	else
	{
		hSlider->setEnabled(false);
		hSpin->setEnabled(false);
		hLabel->setEnabled(false);
	}
	m_T = hSlider->value();
}


//THIS FUNCTION COMPUTES THE HSV VALUES FROM RGB VALUES:
std::vector<double> ImageBrowser5D::RGBtoHSV(std::vector<unsigned char> rgb)
{
	std::vector<double> rVal;
	rVal.assign(3,0);

	if( rgb.size() != 3)
		return rVal;

	double h=0.0,s=0.0,v = 0.0;

	//Get RGB values and normalize:
	double r = (double)rgb.at(0)/255.0;
	double g = (double)rgb.at(1)/255.0;
	double b = (double)rgb.at(2)/255.0;

	double vmax = std::max(std::max(r,g),b);
	double vmin = std::min(std::min(r,g),b);
	double del_max = vmax-vmin;

	v = vmax;		//VALUE

	if( del_max == 0 )
	{
		h = 0;		//HUE
		s = 0;		//SATURATION
	}
	else
	{
		s = del_max / vmax;

		double del_R = ( ( ( vmax - r ) / 6 ) + ( del_max / 2 ) ) / del_max;
		double del_G = ( ( ( vmax - g ) / 6 ) + ( del_max / 2 ) ) / del_max;
		double del_B = ( ( ( vmax - b ) / 6 ) + ( del_max / 2 ) ) / del_max;

		if      ( r == vmax ) h = del_B - del_G;
		else if ( g == vmax ) h = ( 1.0 / 3.0 ) + del_R - del_B;
		else if ( b == vmax ) h = ( 2.0 / 3.0 ) + del_G - del_R;

		if ( h < 0 ) h += 1;
		if ( h > 1 ) h -= 1;
	}
	
	rVal.at(0) = h;
	rVal.at(1) = s;
	rVal.at(2) = v;

	return rVal;
}

void ImageBrowser5D::CreateInteractorStyle(void)
{
	m_kycallback = vtkSmartPointer<vtkCallbackCommand>::New();
	m_kycallback->SetCallback(keyPress);
	m_kycallback->SetClientData(this);
	//m_imageview->GetRenderWindow()->GetInteractor()->SetInteractorStyle(NULL);
	//I want to keep mouse command observers, but change the key ones:
	m_imageview->GetRenderWindow()->GetInteractor()->RemoveObservers(vtkCommand::KeyPressEvent);
	m_imageview->GetRenderWindow()->GetInteractor()->RemoveObservers(vtkCommand::KeyReleaseEvent);
	m_imageview->GetRenderWindow()->GetInteractor()->RemoveObservers(vtkCommand::CharEvent);
	m_imageview->GetRenderWindow()->GetInteractor()->AddObserver(vtkCommand::KeyPressEvent, m_kycallback);
}

void ImageBrowser5D::keyPress(vtkObject * object, unsigned long eid, void *clientdata, void * callerdata)
{
	ImageBrowser5D * client = (ImageBrowser5D*) clientdata;
	int keypressed = client->m_imageview->GetInteractor()->GetKeyCode();
	switch(keypressed)
	{
	case '1': case '2': case '3': case '4': case '5': 
	case '6': case '7': case '8': case '9':
		client->ToggleShowChannel( int(keypressed) - int(48) - 1 );	//Convert key to integer
		break;
	case 'v': case 'V':
		client->SetMode(VOLUME);
		break;
	case 's': case 'S':
		client->SetMode(SLICE);
		break;
	}
}

