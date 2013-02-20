#include <QAction>
#include <QtGui>
#include "Curvelet.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkRescaleIntensityImageFilter.h"

typedef unsigned char TiffPixelType;
const unsigned Dimension = 3;
typedef itk::Image<InputPixelType, Dimension> TiffImageType;

typedef itk::ImageFileReader<TiffImageType> InputReaderType;
typedef itk::ImageFileWriter<TiffImageType> OutputWriterType;

class CurveletGUI : public QMainWindow 
{
Q_OBJECT;
public:
	CurveletGUI(QWidget * parent = 0);
	QListWidget * FileListView;
public slots:
	void BrowseFiles();
	void ProcessFiles();
protected:
	void closeEvent(QCloseEvent *event);
private:
	QSettings CurveletGUISettings;
	QString imageDir;
	QAction *exitAction;
	QAction *loadImages;
	QAction *ProcessImages;
	QDoubleSpinBox * SigmaValue;
	QStringList InputFileList;
};