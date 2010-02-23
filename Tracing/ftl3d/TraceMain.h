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

#include "TraceConfig.h"
#include "SeedContainer3D.h"
#include "Seed2Seg.h"
#include "TraceContainer3D.h"
#include <QAction>
#include <QtGui>
#include <QMainWindow>


typedef	 float PixelType ;
typedef itk::Image< PixelType, 2 >   ImageType2D;
typedef itk::Image< PixelType, 3 >   ImageType3D;
typedef vnl_vector_fixed<double,3> Vect3;
typedef vnl_matrix_fixed <double,3,3> Mat33;
class TraceSEMain: public QMainWindow 
{
Q_OBJECT;
public:
	TraceSEMain();
	void ImageDenoise(ImageType3D::Pointer&, int );
	void UpdateMultiscale( ImageType3D::Pointer& , ImageType3D::Pointer& );
	//void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer& , const float);
	void GetFeature( ImageType3D::Pointer&, ImageType3D::Pointer&);
	void ImageStatistics(ImageType3D::Pointer & );
	void WriteSeedImage(ImageType2D::Pointer, SeedContainer3D::Pointer,std::string );
public slots:
	bool runSETracing();
	void addFileToTrace();
	void GetInputFileName();
	void GetOutputFileName();

private:
	QGroupBox * FileActions;
	void createFileActions();
	QString newInput, newOutput;
	QLabel * InputFileNameLine;
	QLineEdit * OutputFileNameLine;
	QTextEdit * FileListWindow;

	QGroupBox * settingsBox;
	void CreateSettingsLayout();
	std::vector<std::string> InputFileNames;
    std::vector<std::string> OutputFileNames;
	unsigned int numDataFiles;

    //tracing parameters
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
};
#endif
