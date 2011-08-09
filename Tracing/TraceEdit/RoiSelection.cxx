#include "ROISelection.h"

ROISelectionDialog::ROISelectionDialog(double *tedistances, std::string &ImageName, std::vector<double> &cellposxy, QWidget *parent, Qt::WindowFlags flags): QMainWindow(parent,flags){
	distances = tedistances;
	MyImageName = ImageName;
	Mycellpos = cellposxy;
	num_cells = (int)(cellposxy.size()/2);
	QSettings settings("Farsight", "ROI Selection GUI");
	setWindowTitle(tr("ROI Selection GUI"));
	resize(1024, 1024);
	segView = new LabelImageViewQT(NULL);
	CreateMenus();
	CreateMaxIntProj();
	ftk::Image::Pointer myImg = ftk::Image::New();
	myImg->LoadFile("max_pr.tif");
	segView->SetChannelImage(myImg);
	setCentralWidget(segView);
	segView->update();
	raise();
	activateWindow();
}

ROISelectionDialog::~ROISelectionDialog(){
}


void ROISelectionDialog::CreateMaxIntProj(){
	typedef  unsigned char  PixelType;
	typedef itk::Image<PixelType, 3> InputImageType;
	typedef itk::Image<PixelType, 2> OutputImageType;
	typedef itk::ImageFileReader< InputImageType >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	
	reader->SetFileName( MyImageName.c_str() );
	reader->Update();

	InputImageType::Pointer inp_im = reader->GetOutput();

	int size1 = inp_im->GetLargestPossibleRegion().GetSize()[0],
		size2 = inp_im->GetLargestPossibleRegion().GetSize()[1],
		size3 = inp_im->GetLargestPossibleRegion().GetSize()[2],
		small_one = 3;

	if( size1<size3 && size1<size2 )
		small_one = 1;
	if( size2<size3 && size2<=size1 )
		small_one = 2;
	
	typedef itk::ImageRegionIteratorWithIndex< InputImageType > IteratorType;
	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > IteratorType1;

	OutputImageType::Pointer im_out = OutputImageType::New();
	OutputImageType::PointType origin;
	origin[0] = 0;
	origin[1] = 0;
	im_out->SetOrigin( origin );

	OutputImageType::IndexType start;
	start[0] = 0;  // first index on X
	start[1] = 0;  // first index on Y
	OutputImageType::SizeType  sizes;
	sizes[0] = size1;  // size along X
	sizes[1] = size2;  // size along X

	OutputImageType::RegionType region;
	region.SetSize( sizes );
	region.SetIndex( start );

	im_out->SetRegions( region );
	im_out->Allocate();
	im_out->FillBuffer(0);
	im_out->Update();

	InputImageType::RegionType::IndexType outputi;

	IteratorType pix_buf( inp_im, inp_im->GetRequestedRegion() );
	IteratorType1 pix_buf1( im_out, im_out->GetRequestedRegion() );

	for( int k=0; k<size3; ++k ){
		pix_buf1.GoToBegin();
		for( int i=0; i<size1; ++i ){
			for( int j=0; j<size2; ++j ){
				outputi[0] = i;	outputi[1] = j;	outputi[2] = k;
				pix_buf.SetIndex( outputi );
				if(k == 0){
					pix_buf1.Set( pix_buf.Get() );
				}
				else if( pix_buf1.Get() < pix_buf.Get() )
					pix_buf1.Set( pix_buf.Get() );
				++pix_buf1;
			}
		}
	}
	writer->SetInput( im_out );
	writer->SetFileName( "max_pr.tif" );
	writer->Update();
}

void ROISelectionDialog::CreateMenus(void){	
	QMenu *ROIMenu = menuBar()->addMenu("ROI");
	QAction *StROISelAction = new QAction("Start ROI Selection", this);
	connect(StROISelAction, SIGNAL(triggered()), this, SLOT(StartROISel()));
	ROIMenu->addAction(StROISelAction);
	QAction *ClrROISelAction = new QAction("Clear ROI Selection", this);
	connect(ClrROISelAction, SIGNAL(triggered()), this, SLOT(ClearROISel()));
	ROIMenu->addAction(ClrROISelAction);
	QAction *DoneROISelAction = new QAction("Done ROI Selection", this);
	connect(DoneROISelAction, SIGNAL(triggered()), this, SLOT(EndROISel()));
	ROIMenu->addAction(DoneROISelAction);
}

void ROISelectionDialog::StartROISel(void){
	segView->SetROIVisible( true );
	segView->GetROI();
	connect(segView, SIGNAL(roiDrawn()), this, SLOT(EndROISel()));
}

void ROISelectionDialog::ClearROISel(void){
	segView->GetROIMaskImage()->fill(Qt::white);
	segView->SetROIVisible(false);
	segView->update();
}

void ROISelectionDialog::EndROISel(void){
	segView->ClearGets();
	disconnect(segView, SIGNAL(roiDrawn()), this, SLOT(EndROISel()));
	UpdatedTable();

}

void ROISelectionDialog::UpdatedTable(void){
	QImage * t_img = segView->GetROIMaskImage();
	if(t_img == NULL)
		return;

	typedef itk::Image<unsigned char, 2> ImageType;

	//Change QImage into ITK IMAGE:
	ImageType::Pointer roiImg = ImageType::New();
	ImageType::PointType origin;
   	origin[0] = 0; 
	origin[1] = 0;
    roiImg->SetOrigin( origin );
	ImageType::IndexType start = {{ 0,0 }};    
	ImageType::SizeType size = {{ t_img->width(), t_img->height() }};
    ImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( start );
    roiImg->SetRegions( region );
	roiImg->Allocate();
	roiImg->FillBuffer(0);
	for(int i=0; i<t_img->width(); ++i){
		for(int j=0; j<t_img->height(); ++j){
			ImageType::IndexType ind = {{ i,j }};
			roiImg->SetPixel(ind, t_img->pixelIndex( i, j ) );
		}
	}

	//Testing
	/*typedef itk::ImageFileWriter< ImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( roiImg );
	writer->SetFileName( "test.tif" );
	writer->Update();*/

	typedef itk::Image< int, 2 > IntImageType;
	typedef itk::Image< float, 2 > FltImageType;
	typedef itk::BinaryThresholdImageFilter< ImageType, IntImageType > BinaryThresholdFilterType;
	BinaryThresholdFilterType::Pointer binarythreshfilter = BinaryThresholdFilterType::New();
	binarythreshfilter->SetInsideValue( INT_MAX );
	binarythreshfilter->SetOutsideValue( 0 );
	binarythreshfilter->SetLowerThreshold( 1 );
	binarythreshfilter->SetUpperThreshold( INT_MAX );
	binarythreshfilter->SetInput( roiImg );
	binarythreshfilter->Update();

	typedef itk::SignedMaurerDistanceMapImageFilter< IntImageType, FltImageType > SignedMaurerDistanceMapFilterType;
	SignedMaurerDistanceMapFilterType::Pointer distancemapfilter = SignedMaurerDistanceMapFilterType::New();
	distancemapfilter->SetInput(binarythreshfilter->GetOutput());
	distancemapfilter->SquaredDistanceOff();
	distancemapfilter->Update();

/*	//Testing
	typedef itk::RescaleIntensityImageFilter< FltImageType, ImageType > RescaleFltIType;
	RescaleFltIType::Pointer RescaleIntIO1 = RescaleFltIType::New();
	RescaleIntIO1->SetOutputMaximum( 255 );
	RescaleIntIO1->SetOutputMinimum( 0 );
	RescaleIntIO1->SetInput( distancemapfilter->GetOutput() ); //watershedfilter->GetOutput() image1
	RescaleIntIO1->Update();

	typedef itk::ImageFileWriter< ImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetInput( RescaleIntIO1->GetOutput() );
	writer->SetFileName( "test.jpg" );
	writer->Update();
*/
	std::string ColName = "DistanceToDevice";
	FltImageType::Pointer dImg = distancemapfilter->GetOutput();

	typedef itk::ImageRegionIteratorWithIndex< FltImageType > IteratorType;
	IteratorType pix_buf( dImg, dImg->GetRequestedRegion() );
	FltImageType::RegionType::IndexType outputi;
	for( int i=0; i<num_cells; ++i ){
		outputi[0] = (int)Mycellpos.at(i*2); outputi[1] = (int)Mycellpos.at(i*2+1);
		pix_buf.SetIndex( outputi );
		if( pix_buf.Get()>0 )
			distances[i] = pix_buf.Get();
		else
			distances[i] = 0;
	}
	remove("max_pr.tif");
	QMainWindow::close();
}

void ROISelectionDialog::closeEvent(QCloseEvent *event){
	emit dialogClosed();
	delete segView;
	MyImageName.clear();
	event->accept();
} 