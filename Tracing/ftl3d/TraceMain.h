/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
=========================================================================*/

#ifndef TRACEMAIN_H_
#define TRACEMAIN_H_

#include <iostream>
#include <fstream>

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkMedianImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkArray.h"
#include "itkCovariantVector.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkSymmetricEigenAnalysis.h"
#include "vnl/vnl_math.h"
#include <vnl/vnl_vector_fixed.h>

//#include "TraceConfig.h" //no longer used
#include "SeedContainer3D.h"
#include "Seed2Seg.h"
#include "TraceContainer3D.h"
#include "TraceNode.h"

#include <ftkCommon/ftkProjectManager.h>
#include <ftkCommon/ftkParameters.h>
#include "tinyxml/tinyxml.h"

#include <QAction>
#include <QtGui>
#include <QMainWindow>
#include <QWidget>
#include <QString>
#include <QObject>


class TraceSEMain : public QMainWindow 
{
	Q_OBJECT;

public:
	TraceSEMain(QWidget * parent = 0);
	~TraceSEMain(){};
	void ImageDenoise(ImageType3D::Pointer&, int );
	void UpdateMultiscale( ImageType3D::Pointer& , ImageType3D::Pointer& );
	//void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer& , const float);
	void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer&);
	void ImageStatistics(ImageType3D::Pointer & );
	void WriteSeedImage(ImageType2D::Pointer, SeedContainer3D::Pointer,std::string );
	void SetProjectName(QString projName);
	void LoadTraceSettingsFile(std::string filename);
public slots:
	bool runSETracing();
	void addFileToTrace();
	void GetInputFileName();
	void GetOutputFileName();
	void LoadFromTraceProject();
	void FileSuffixChanged(QString s);

private:
	QSettings TraceSESettings;
	bool ReadNodeXMLFile(std::string xmlfname, std::vector<TraceNode*>& NodeContainer);
	void WriteSWCFile(std::string SWCFilename, const std::vector<TraceNode*>& NodeContainer);
	void CreateSettingsLayout();
	void createFileActions();
//define  structures
	typedef	 float PixelType ;
	typedef itk::Image< PixelType, 2 >   ImageType2D;
	typedef itk::Image< PixelType, 3 >   ImageType3D;
	typedef vnl_vector_fixed<double,3> Vect3;
	typedef vnl_matrix_fixed <double,3,3> Mat33;
	ftk::ProjectManager * Project;
//GUI objects
	QGroupBox * FileActions;
	QGroupBox * settingsBox;

	QComboBox * fileSuffixBox;
	QLabel * InputFileNameLine;
	QLineEdit * OutputFileNameLine;
	QTextEdit * FileListWindow;
	unsigned int numDataFiles;
//file names
	QString newInput, newOutput, FileSuffix, ProjectName;
	std::vector<std::string> InputFileNames;
    std::vector<std::string> OutputFileNames;
    std::vector<std::string> OutputSWCFileNames;


//tracing parameters and gui object
	int GridSpacing;
	QSpinBox * GetGridSpacing;
	double StepRatio;
	QDoubleSpinBox * GetStepRatio;
	double AspectRatio;
	QDoubleSpinBox * GetAspectRatio;
	double THRESHOLD;
	QDoubleSpinBox * GetTHRESHOLD;
	double minContrast;
	QDoubleSpinBox * GetminContrast;
	double MaximumVesselWidth;
	QDoubleSpinBox * GetMaximumVesselWidth;
	double MinimumVesselLength;
	QDoubleSpinBox * GetMinimumVesselLength;
	int FitIterations;
	QSpinBox * GetFitIterations;
	double MinimumVesselWidth;
	QDoubleSpinBox * GetMinimumVesselWidth;
	double StartTHRESHOLD;
	QDoubleSpinBox * GetStartTHRESHOLD;
	itk::FixedArray<double, 3> Spacing;
	int UseHessian;
	QCheckBox * GetUseHessian;
	bool ConvertToSWC;
	QCheckBox * GetConvertToSWC;
};

#endif
