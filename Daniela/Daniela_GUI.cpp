#include "Daniela_GUI.h"

Daniela_GUI::Daniela_GUI(QWidget * parent, Qt::WindowFlags flags) : QMainWindow(parent,flags)
{
	QSettings settings("Farsight", "Daniela GUI");
	setWindowTitle(tr("Daniela GUI"));
	resize(1024, 1024);

	selection = new ObjectSelection();
	//QMap<QString, QColor> colorItemsMap;
	segView = new LabelImageViewQT(NULL);
	createMenus();
	segView->SetChannelImage(NULL);
	setCentralWidget(segView);
	segView->update();

	raise();
	activateWindow();
}


void Daniela_GUI::loadImage(void)
{
	QString lastPath = settings.value("lastPath", ".").toString();
	//std::cout << lastPath.toStdString() << std::endl;
	QString fileName = QFileDialog::getOpenFileName(this, "Open Image", lastPath, "*.*");
	ftk::Image::Pointer myImg = ftk::Image::New();
	myImg->LoadFile(fileName.toStdString());
	segView->SetChannelImage(myImg);
	setCentralWidget(segView);
	segView->update();
}


void Daniela_GUI::createMenus(void)
{
	
	QMenu *fileMenu = menuBar()->addMenu("File");
	QAction *loadImageAction = new QAction("Load Image", this);
	connect(loadImageAction, SIGNAL(triggered()), this, SLOT(loadImage()));
	fileMenu->addAction(loadImageAction);
	QMenu *fileMenu2 = menuBar()->addMenu("Load Result");
	QAction *loadImageAction2 = new QAction("Load label Image", this);
	connect(loadImageAction2, SIGNAL(triggered()), this, SLOT(loadResult()));
	fileMenu2->addAction(loadImageAction2);


}

void Daniela_GUI::loadResult(void)
{
	QSettings settings;
	QString lastPath = settings.value("lastPath", ".").toString();
	QString originalImageFileName = QFileDialog::getOpenFileName(this, "Open Original Image", lastPath, "*.*");
	QString labelImageFileName = QFileDialog::getOpenFileName(this, "Open Label Image", lastPath, "*.*");
	//std::cout << fileName.toStdString() << std::endl;
	ftk::Image::Pointer labelImg = ftk::Image::New();
	ftk::Image::Pointer origImg = ftk::Image::New();
	labelImg->LoadFile(labelImageFileName.toStdString());
	origImg->LoadFile(originalImageFileName.toStdString());
	/*QMap<QString, QColor> colorItemsMap;
	segView = new LabelImageViewQT(&colorItemsMap);*/
	selection->clear();
	segView->SetChannelImage(origImg);
	segView->SetLabelImage(labelImg,selection);
	segView->update();
}
    