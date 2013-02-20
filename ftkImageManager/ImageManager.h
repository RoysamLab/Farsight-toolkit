#include <QAction>
#include <QtGui>
#include <QFile>
#include <QTextStream>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"
#include "ftkPreprocess2.h"

#include "ftkGUI/PreprocessDialog.h"

typedef unsigned short InputPixelType;
typedef unsigned char OutputPixelType;
const unsigned Dimension = 3;
typedef itk::Image<InputPixelType, Dimension> InputImageType;
typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

typedef itk::ImageFileReader<InputImageType> InputReaderType;
typedef itk::ImageFileWriter<OutputImageType> OutputWriterType;
typedef itk::RescaleIntensityImageFilter<InputImageType, OutputImageType> RescaleFilterType;


class ImageFileManger : public QMainWindow 
{
Q_OBJECT;
public:
	ImageFileManger(QWidget * parent = 0);
	QListWidget * FileListView;
	std::vector<QString> readDataFile(QString FileName);
public slots:
	void BrowseFiles();
	void ProcessFiles();
	void appendLists();
	void Preprocess();
protected:
	void closeEvent(QCloseEvent *event);
	PreprocessDialog *preprocessdialog;
private:
	QSettings ImageManagerSettings;
	QString imageDir;
	QAction *exitAction;
	QAction *loadImages;
	QAction *append;
	QAction *ProcessImages;
	QAction *StartPreprocessing;
	QStringList InputFileList;
	QStringList outputDirectories;
};