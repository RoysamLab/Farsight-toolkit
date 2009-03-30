#include "SliceView5D.h"

//#include "itkImageRegionConstIterator.h"

//Destructor
SliceView5D::~SliceView5D()
{
	m_vtkimageactor.clear();
	if(img) delete img;
	if(m_imageview) delete m_imageview;
	if(vSlider) delete vSlider;
}

//Constructor
SliceView5D::SliceView5D(QString filename)
{
	//CREATE A QVTK Widget for putting a QVTK renderer in QT window
	m_imageview = new QVTKWidget();

	//Load up an image:
	img = new ftk::Image();
	bool forDisplay = true;
	if( !img->LoadFile(filename.toStdString(), forDisplay) )	return;

	/***************************************************************
	// TEST FOR CHARLENE
	typedef itk::Image<unsigned char, 3> ImageType;
	ImageType::Pointer i = img->GetItkPtr<unsigned char>(0,0);

	//Try iterating through image
	typedef itk::ImageRegionConstIterator<ImageType> IteratorType;
	IteratorType inputIt(i, i->GetLargestPossibleRegion() );

	for(inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
	{
		unsigned char val = inputIt.Get();
	}
	*******************************************************************/

	//Create the renderer and add the image to it.
	m_vtkrenderer = RendererPointerType::New();

	int chs = img->GetImageInfo()->numChannels;
	for(int i=0; i<chs; ++i)
	{
		vtkSmartPointer<vtkLookupTable> table = vtkSmartPointer<vtkLookupTable>::New();
		table->SetRange(0, 255);				// image intensity range
		
		std::vector<double> hsv = RGBtoHSV( img->GetImageInfo()->channelColors.at(i) );
		//S and V should always be 1, it is H that I am most concerned about

		table->SetHueRange(hsv[0],hsv[0]);			// Split Range
		table->SetSaturationRange(hsv[1],hsv[1]);	// Full Color Saturation
		table->SetValueRange(hsv[2],hsv[2]);		// from black to white
		table->SetAlphaRange(0,1);
		table->Build();

		vtkSmartPointer<vtkImageMapToColors> color = vtkSmartPointer<vtkImageMapToColors>::New();
		color->SetInput( img->GetVtkPtr(0,i) );
		color->SetLookupTable(table);

		//Create the actor and set its input to the image
		ImageActorPointerType actor = ImageActorPointerType::New();
		actor->SetInput( color->GetOutput() );
		actor->SetDisplayExtent(0, img->GetImageInfo()->numColumns - 1, 0, img->GetImageInfo()->numRows - 1, 0, 0 );
		actor->SetZSlice(0);

		m_vtkimageactor.push_back(actor);
		m_vtkrenderer->AddActor(actor);
	}
	
	m_vtkrenderer->SetBackground(0.0,0.0,0.0);
	
	//Create a new Interactor Style to use for the image:
	vtkInteractorStyleImage *style = vtkInteractorStyleImage::New();
	m_imageview->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);

	//Set the renderer & render
	m_imageview->GetRenderWindow()->AddRenderer(m_vtkrenderer);
	m_imageview->GetRenderWindow()->Render();

	//Setup Vertical slider widget
    vSlider = new QSlider();
    vSlider->setOrientation(Qt::Vertical);
	vSlider->setDisabled(false);
	vSlider->setRange(0,img->GetImageInfo()->numZSlices-1);
	vSlider->setValue(0);
	connect(vSlider, SIGNAL(valueChanged(int)), this, SLOT(setZ(int)));

	QGridLayout *viewerLayout = new QGridLayout();
    viewerLayout->addWidget(m_imageview, 0, 0);
    viewerLayout->addWidget(vSlider, 0, 1);
	viewerLayout->setColumnMinimumWidth(0, img->GetImageInfo()->numColumns);
	viewerLayout->setRowMinimumHeight(0, img->GetImageInfo()->numRows);

	this->setLayout(viewerLayout);
	//this->resize(img->GetImageInfo()->numColumns, img->GetImageInfo()->numRows);
	this->setAttribute ( Qt::WA_DeleteOnClose );
	this->setWindowTitle(tr("Image Browser"));
}

void SliceView5D::setZ(int z)
{
	for (int i=0; i<(int)m_vtkimageactor.size(); ++i)
	{
		m_vtkimageactor.at(i)->SetZSlice(z);
	}

	m_vtkrenderer->ResetCameraClippingRange();	//This fixes the "disappearing image" problem
	//m_vtkrenderer->ResetCamera();					//This will reset the "zoom" as well - not good
	m_imageview->GetRenderWindow()->Render();

}



//THIS FUNCTION COMPUTES THE HSV VALUES FROM RGB VALUES:
std::vector<double> SliceView5D::RGBtoHSV(std::vector<unsigned char> rgb)
{
	std::vector<double> rVal;
	rVal.assign(3,0);

	if( rgb.size() != 3)
		return rVal;

	double h,s,v;

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