#include <QAction>
#include <QtGui>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"

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
public slots:
	void BrowseFiles();
	void ConvertFiles();
protected:
	void closeEvent(QCloseEvent *event);
private:
	QSettings ImageManagerSettings;
	QString imageDir;
	QAction *exitAction;
	QAction *loadImages;
	QAction *ConvertImages;
	QStringList InputFileList;
	QStringList outputDirectories;
};