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

	this->append = new QAction("Append L Measure txt", this);
	connect(this->append, SIGNAL(triggered()), this, SLOT(appendLists()));
	this->menuBar()->addAction(this->append);

	this->StartPreprocessing = new QAction("Preprocess..", this);
	connect(this->StartPreprocessing, SIGNAL(triggered()), this, SLOT(Preprocess()));
	this->StartPreprocessing->setEnabled(false);
	this->menuBar()->addAction(this->StartPreprocessing);

	preprocessdialog = new PreprocessDialog("",this);

	this->ProcessImages = new QAction("Process Images", this);
	connect(this->ProcessImages, SIGNAL(triggered()), this, SLOT(ProcessFiles()));
	this->ProcessImages->setEnabled(false);
	this->menuBar()->addAction(this->ProcessImages);

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
		this->ProcessImages->setEnabled(true);
		this->StartPreprocessing->setEnabled(true);
		this->imageDir = QFileInfo(this->InputFileList[0]).absolutePath();
		this->ImageManagerSettings.setValue("imageDir", imageDir);
		this->FileListView->addItems(this->InputFileList);
		this->ImageManagerSettings.sync();
	}
}

void ImageFileManger::Preprocess(void)
{

	if(!preprocessdialog->exec())
	{
		std::cerr << "Preprocess dialog canceled." << std::endl;

	}
}

void ImageFileManger::ProcessFiles()
{
	if (!this->InputFileList.isEmpty())
	{
		QDir directory(this->imageDir);
		for (int k = 0; k < this->outputDirectories.size(); k++)
		{
			directory.mkdir(this->outputDirectories[k]);
		}
		QProgressDialog progress("Converting Images", "Abort", 0, this->InputFileList.size(), this);
		progress.setWindowModality(Qt::WindowModal);
		for (int i = 0; i < this->InputFileList.size(); i++)
		{
			progress.setValue(i);
			if (progress.wasCanceled())
			{
				break;
			}
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

			InputReaderType::Pointer reader = InputReaderType::New();
			OutputWriterType::Pointer writer = OutputWriterType::New();
			reader->SetFileName(this->InputFileList[i].toStdString().c_str());

			//16 to 8 bit conversion
			RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
			rescaleFilter->SetInput(reader->GetOutput());
			rescaleFilter->SetOutputMinimum(0);
			rescaleFilter->SetOutputMaximum(255);
			rescaleFilter->Update();

			preprocessdialog->SetImage( rescaleFilter->GetOutput() );
			
			writer->SetFileName(outputFilename.toStdString().c_str());
			try
			{
				preprocessdialog->Process();
				writer->SetInput(preprocessdialog->GetImage());
				writer->Update();
			}
			catch(itk::ExceptionObject & e) 
			{
				std::cerr << "Exception in ITK Pipeline: " << e << std::endl;
				continue;
			}
		}
		progress.setValue(this->InputFileList.size());
	}
	
	
}
std::vector<QString> ImageFileManger::readDataFile(QString FileName)
{
	QFile inFile(FileName);
	if (!inFile.open(QIODevice::ReadOnly | QIODevice::Text))
	{
	}
	std::vector<QString> Lines;
	QTextStream in(&inFile);
	QString lineIn = in.readLine();
	while (!lineIn.isNull()) 
	{
		QString newList = lineIn;//.split("\t");
		Lines.push_back(newList);
		lineIn = in.readLine();
	}
	inFile.close();
	return Lines;
}
void ImageFileManger::appendLists()
{
	QString MainFileName = QFileDialog::getOpenFileName(this , "Choose file", "", tr("Text ( *.txt ) "));
	QString appendFileName = QFileDialog::getOpenFileName(this , "Choose file", "", tr("Text ( *.txt ) "));
	std::vector<QString> MainList = this->readDataFile(MainFileName);
	std::vector<QString> AppendList = this->readDataFile(appendFileName);
	for (int i = 0; i < MainList.size() ; i++)
	{
		for (int j = 0; j < AppendList.size(); j++)
		{
			QString appendString = AppendList.at(j);
			QString key = QString(AppendList.at(j).section('\t', 0, 0) + "\t");
			//std::cout<< "Key\t"<<j << "\t" << key.toStdString().c_str() << "\n";
			if (MainList.at(i).contains(key))
			{
				std::cout<< "matched " << key.toStdString().c_str() << "\n";
				appendString.remove(key);
				MainList.at(i).append(appendString);
				continue;
			}
		}//end j loop
	}//end for loop
	QFile outFile(MainFileName);
	if(outFile.open(QIODevice::WriteOnly))
	{
		QTextStream out(&outFile);
		for (int i = 0; i < MainList.size() ; i++)
		{
			out << MainList.at(i) << "\n";
		}
	}
	outFile.close();
}
void ImageFileManger::closeEvent(QCloseEvent *event)
{
	event->accept();
}
