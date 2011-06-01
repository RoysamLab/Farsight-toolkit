#include "ImageManager.h"
ImageFileManger::ImageFileManger(QWidget *parent)
: QMainWindow(parent)
{
	this->FileListView = new QListWidget;
	this->FileListView->isSortingEnabled();
	this->setCentralWidget(this->FileListView);
	this->InputFileList.clear();

	this->exitAction = new QAction(tr("Exit"), this);
	connect(this->exitAction, SIGNAL(triggered()), this, SLOT(close()));
	this->exitAction->setShortcut(QKeySequence::Close);
	this->menuBar()->addAction(this->exitAction);

	this->loadImages = new QAction("Load Images", this);
	connect(this->loadImages, SIGNAL(triggered()), this, SLOT(BrowseFiles()));
	this->menuBar()->addAction(this->loadImages);

	this->ConvertImages = new QAction("Convert Images", this);
	connect(this->ConvertImages, SIGNAL(triggered()), this, SLOT(ConvertFiles()));
	this->ConvertImages->setEnabled(false);
	this->menuBar()->addAction(this->ConvertImages);

	this->outputDirectories;
	this->outputDirectories << "DAPI" << "Cy5"<< "TRITC"<<"GFP";
}
void ImageFileManger::BrowseFiles()
{
	this->imageDir = this->ImageManagerSettings.value("imageDir", ".").toString();
	this->InputFileList = QFileDialog::getOpenFileNames(this , "Choose Images", this->imageDir, 
	  tr("Image File ( *.tiff *.tif *.pic *.PIC ) "));
	if(!this->InputFileList.isEmpty())
	{
		this->ConvertImages->setEnabled(true);
		this->imageDir = QFileInfo(this->InputFileList[0]).absolutePath();
		this->ImageManagerSettings.setValue("imageDir", imageDir);
		this->FileListView->addItems(this->InputFileList);
		this->ImageManagerSettings.sync();
	}
}
void ImageFileManger::ConvertFiles()
{
	if (!this->InputFileList.isEmpty())
	{
		//QString newDir = QString("DAPI");
		QDir directory(this->imageDir);
		for (int k = 0; k < this->outputDirectories.size(); k++)
		{
			directory.mkdir(this->outputDirectories[k]);
		}
		for (int i = 0; i < this->InputFileList.size(); i++)
		{
			QString currentDir;
			currentDir.clear();
			QString tempname = QFileInfo(this->InputFileList[i]).fileName();
			for (int j = 0; j < this->outputDirectories.size(); j++)
			{
				if (tempname.contains(this->outputDirectories.at(j), Qt::CaseInsensitive))
				{
					currentDir = this->outputDirectories.at(j);
				}//end search
			}
			if (currentDir.isEmpty())
			{
				continue;
			}
			tempname.replace(" ", "_");
			tempname.prepend("8Bit");
			QString outputFilename = QString(this->imageDir +"/"+ currentDir +"/"+tempname);
			std::cout<< "Output name\t" <<outputFilename.toStdString().c_str()<< "\n";


			InputReaderType::Pointer reader = InputReaderType::New();
			OutputWriterType::Pointer writer = OutputWriterType::New();
			reader->SetFileName(this->InputFileList[i].toStdString().c_str());
			reader->Update();

			RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
			rescaleFilter->SetInput(reader->GetOutput());
			rescaleFilter->SetOutputMinimum(0);
			rescaleFilter->SetOutputMaximum(255);
			rescaleFilter->Update();

			writer->SetInput(rescaleFilter->GetOutput());
			writer->SetFileName(outputFilename.toStdString().c_str());
			writer->Update();
		}
	}
	
	//16 to 8 bit conversion
}
void ImageFileManger::closeEvent(QCloseEvent *event)
{
	event->accept();
}
