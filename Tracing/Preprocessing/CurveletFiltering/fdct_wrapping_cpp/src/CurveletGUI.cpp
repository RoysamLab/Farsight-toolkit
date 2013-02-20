#include "CurveletGUI.h"
CurveletGUI::CurveletGUI(QWidget *parent)
: QMainWindow(parent)
{
	this->resize(800,600);
	this->FileListView = new QListWidget;
	this->FileListView->isSortingEnabled();
	this->setCentralWidget(this->FileListView);
	this->InputFileList.clear();
	QToolBar * mainToolBar = new QToolBar(this);

	this->exitAction = new QAction(tr("Exit"), this);
	connect(this->exitAction, SIGNAL(triggered()), this, SLOT(close()));
	this->exitAction->setShortcut(QKeySequence::Close);
	this->menuBar()->addAction(this->exitAction);

	this->loadImages = new QAction("Load Images", this);
	connect(this->loadImages, SIGNAL(triggered()), this, SLOT(BrowseFiles()));
	this->menuBar()->addAction(this->loadImages);

	this->ProcessImages = new QAction("Run Curvelets on Images", this);
	connect(this->ProcessImages, SIGNAL(triggered()), this, SLOT(ProcessFiles()));
	this->ProcessImages->setEnabled(false);
	this->menuBar()->addAction(this->ProcessImages);

	this->SigmaValue = new QDoubleSpinBox(this);
	this->SigmaValue->setRange(0,1);
	this->SigmaValue->setValue(.03);
	this->SigmaValue->setSingleStep(.01);
	mainToolBar->addWidget(new QLabel("Sigma Value: "));
	mainToolBar->addWidget(this->SigmaValue);
	this->addToolBar(mainToolBar);
}

void CurveletGUI::BrowseFiles()
{
	this->imageDir = this->CurveletGUISettings.value("imageDir", ".").toString();
	this->InputFileList = QFileDialog::getOpenFileNames(this , "Choose Images", this->imageDir, 
	  tr("Image File ( *.tiff *.tif *.pic *.PIC ) "));
	if(!this->InputFileList.isEmpty())
	{
		this->FileListView->addItems(this->InputFileList);
		this->ProcessImages->setEnabled(true);
		this->imageDir = QFileInfo(this->InputFileList[0]).absolutePath();
		this->CurveletGUISettings.setValue("imageDir", imageDir);
		this->CurveletGUISettings.sync();
	}
}
void CurveletGUI::ProcessFiles()
{
	if (!this->InputFileList.isEmpty())
	{
		float sigma = (float) this->SigmaValue->value();
		QProgressDialog progress("Converting Images", "Abort", 0, this->InputFileList.size(), this);
		progress.setWindowModality(Qt::WindowModal);
		for (int i = 0; i < this->InputFileList.size(); i++)
		{
			progress.setValue(i);
			if (progress.wasCanceled())
			{
				break;
			}

			QString outputname = this->InputFileList[i];
			outputname.insert(outputname.size() -4, QString("_cv"));

			InputReaderType::Pointer reader = InputReaderType::New();
			OutputWriterType::Pointer writer = OutputWriterType::New();
			reader->SetFileName(this->InputFileList[i].toStdString().c_str());
			reader->Update();

			Curvelet curvletfilter = Curvelet();
			curvletfilter.SetSigma(sigma);
			TiffImageType::Pointer outputim = curvletfilter.RunOnInputImage(reader->GetOutput());
			std::cout<< "writing "<<outputname.toStdString().c_str() << std::endl;

			writer->SetInput(outputim);
			writer->SetFileName(outputname.toStdString().c_str());
			writer->Update();
		}
		progress.setValue(this->InputFileList.size());
	}//end file list empty
}
void CurveletGUI::closeEvent(QCloseEvent *event)
{
	event->accept();
}
