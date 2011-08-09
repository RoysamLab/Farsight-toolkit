#ifndef ROISELECTION_H
#define ROISELECTION_H

#include "ftkGUI/LabelImageViewQT.h"

#include <QtGui/QMainWindow>
#include <QtGui/QMenubar>

#include <ftkImage/ftkImage.h>
#include <ftkCommon/ftkUtils.h>

#include "itkImageRegionIteratorWithIndex.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

class ROISelectionDialog : public QMainWindow
{
	Q_OBJECT
public:
	ROISelectionDialog(std::vector<double> &distances, std::string &ImageName, std::vector<double> &cellposxy, QWidget *parent, Qt::WindowFlags flags = 0);
	~ROISelectionDialog();

private slots:
	void StartROISel(void);
	void ClearROISel(void);
	void EndROISel(void);
	void closeEvent(QCloseEvent *event);

public:
	LabelImageViewQT *segView;

private:
	void CreateMenus(void);
	void LoadImage(std::string file_name);
	void UpdatedTable(void);
	void CreateMaxIntProj(void);

	QSettings settings;

	std::vector<double> distances;
	std::string MyImageName;
	ftk::Image::Pointer MyImg;
	std::vector<double> Mycellpos;
	int num_cells;

signals:
	void dialogClosed();

};
#endif