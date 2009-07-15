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

#ifndef __HELPERS_H
#define __HELPERS_H

#pragma warning(disable: 4996)
#pragma warning(disable: 4018)

//stl includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>



//standard c++ includes
#include <stdio.h>
#include <stdlib.h>

//ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkConstNeighborhoodIterator.h>
#include "itkImageToVTKImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "ftkLabelImageToFeatures.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkScalarImageTextureCalculator.h"

//VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkContourFilter.h"
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include <vtkOpenGLVolumeTextureMapper2D.h>
#include <vtkOpenGLVolumeTextureMapper3D.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkVolume.h>
#include <vtkImageImport.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkImageData.h>

//Macros
#define MAX(a,b) (((a) > (b))?(a):(b))
#define MIN(a,b) (((a) < (b))?(a):(b))
#define CLAMP(a,min,max) (((a)<(min))?(min):(((a)>(max))?(max):(a)))
#define DEBUG1 
#define DEBUG2 printf
#define DEBUG3
#define PROGRESS printf

//typdefs
typedef unsigned char InputPixelType;
typedef unsigned char OutputPixelType;

typedef itk::Image<InputPixelType,3> InputImageType;
typedef itk::Image<InputPixelType,2> Input2DImageType;
typedef itk::Image<short,3> LabelImageType;
typedef itk::Image<short,2> Label2DImageType;

typedef itk::Vector<unsigned char, 3> VectorPixelType;
typedef itk::Image<VectorPixelType, 3> ColorImageType;
typedef itk::Image<VectorPixelType, 2> Color2DImageType;

typedef itk::ImageRegionConstIterator<LabelImageType> ConstLabelIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelIteratorType;
typedef itk::ImageRegionIterator<InputImageType> IteratorType;


typedef itk::ImageLinearIteratorWithIndex< Color2DImageType > LinearColorIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< ColorImageType > SliceColorIteratorType;

typedef itk::ImageLinearIteratorWithIndex< Input2DImageType > LinearIteratorType;
typedef itk::ImageSliceConstIteratorWithIndex< InputImageType > SliceIteratorType;;

typedef itk::ImageRegionConstIterator<Input2DImageType> Const2DIteratorType;
typedef itk::ImageRegionIterator<Input2DImageType> twoDIteratorType;

typedef itk::ImageToVTKImageFilter<InputImageType> ConnectorType;
typedef itk::ImageToVTKImageFilter<Input2DImageType> Connector2DType;

//structure declarations
struct TrackPoint{
	double x,y,z;
	int t,id;
};

struct cubecoord{
	unsigned short sx,sy,sz;
	unsigned short ex,ey,ez;
};

struct Feature{
	int x,y,z,t;
	int volume,tag;
	int avg_intensity,min_intensity,max_intensity,sd_intensity;
	int num;
	cubecoord bbox;
	float xaxis_length,yaxis_length,zaxis_length;
	float texture_energy,texture_entropy,texture_inv_diff_moment,texture_inertia, texture_cluster_shade, texture_cluster_prominence;

	void Print(FILE *fp = stdout)
	{
		fprintf(fp,"\n\n");
		fprintf(fp,"num = %d x = %d y = %d z = %d volume = %d\navg_intensity = %d min_intensity = %d max_intensity = %d\n",num, x,y,z,volume,avg_intensity,min_intensity,max_intensity);
		fprintf(fp,"axis_lengths = %0.3f %0.3f %0.3f\n",xaxis_length, yaxis_length, zaxis_length);
		fprintf(fp,"texture: energy %0.3f entropy %0.3f inv_diff_moment %0.3f inertia %0.3f cluster_shade %0.3f cluster_prominence %0.3f\n",texture_energy, texture_entropy, texture_inv_diff_moment, texture_inertia, texture_cluster_shade, texture_cluster_prominence);
		fprintf(fp,"\n\n");
	}
	//double diff(Feature &other)
	//{
	//	double sum;
	//	sum = sqrt(double((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z)));
	//	if(sum > 35)
	//	{
	//		return 1e10;
	//	}

	//	sum = 1-exp(-double((x-other.x)*(x-other.x)/1000.0+(y-other.y)*(y-other.y)/1000.0+(z-other.z)*(z-other.z)/300.0+(intensity-other.intensity)*(intensity-other.intensity)/400.0+(volume-other.volume)*(volume-other.volume)/200000.0));
	//	if(sum<0)
	//	{
	//		Print();
	//		other.Print();
	//	}
	//	printf("Returned SUM = %0.3lf\n",sum);
	//	return sum;
	//}
	//double diff_withoverlap(Feature &other)
	//{
	//	double sum;
	//	sum = sqrt(double((x-other.x)*(x-other.x)+(y-other.y)*(y-other.y)+(z-other.z)*(z-other.z)));
	//	if(sum > 35)
	//	{
	//		return 1e10;
	//	}
	//	int overlap = box_overlap(other);
	//	if(overlap==0)
	//		return 1e10;
	//	double val = overlap*1.0/(0.5*(box_volume()+other.box_volume()));

	//	sum = (1.0-val)*(1-exp(-double((x-other.x)*(x-other.x)/1000.0+(y-other.y)*(y-other.y)/1000.0+(z-other.z)*(z-other.z)/300.0+(intensity-other.intensity)*(intensity-other.intensity)/2/100+(volume-other.volume)*(volume-other.volume)/200000.0)));
	//	
	//	return sum;
	//}
	int box_volume()
	{
		return (static_cast<int>(bbox.ex-bbox.sx)*(bbox.ey-bbox.sy)*(bbox.ez-bbox.sz));
	}
	int box_overlap(Feature &other)
	{
		cubecoord o;
		o.sx = MAX(bbox.sx,other.bbox.sx);
		o.sy = MAX(bbox.sy,other.bbox.sy);
		o.sz = MAX(bbox.sz,other.bbox.sz);
		o.ex = MIN(bbox.ex,other.bbox.ex);
		o.ey = MIN(bbox.ey,other.bbox.ey);
		o.ez = MIN(bbox.ez,other.bbox.ez);
		if((o.ex>o.sx) && (o.ey>o.sy) && (o.ez > o.sz))
			return ((int)(o.ex-o.sx))*(o.ey-o.sy)*(o.ez-o.sz);
		else
			return 0;
	}
};

//function declarations

template <typename T> typename T::Pointer readImage(const char*);
template <typename T> typename int writeImage(typename T::Pointer,const char*);
//Color2DImageType::Pointer getColorBoundaryImage(LabelImageType::Pointer labelled, InputImageType::Pointer im, int channel);
vtkSmartPointer<vtkPolyData> get2DBoundary(LabelImageType::Pointer label);
Input2DImageType::Pointer get2DBoundaryImage(LabelImageType::Pointer label);
InputImageType::Pointer getEmpty(int s1,int s2, int s3);
Input2DImageType::Pointer get2DEmpty(int s1, int s2);
Input2DImageType::Pointer getProjection(InputImageType::Pointer);
vtkSmartPointer<vtkActor> getActorForPolyData(vtkSmartPointer<vtkPolyData>);

template InputImageType::Pointer readImage<InputImageType>(const char*);
template LabelImageType::Pointer readImage<LabelImageType>(const char*);
template int writeImage<LabelImageType>(LabelImageType::Pointer, const char *);

void getFeatureVectors(LabelImageType::Pointer im,InputImageType::Pointer in_image,std::vector<Feature> &feature_vector,int time,int tag);
void getFeatureVectorsNew(LabelImageType::Pointer im,InputImageType::Pointer in_image,std::vector<Feature> &feature_vector,int time,int tag);
vtkSmartPointer<vtkTextActor> getActorForFeature(ftk::LabelImageFeatures f);
vtkSmartPointer<vtkPolyData> getRectangle(double x1, double y1, double x2, double y2);
Input2DImageType::Pointer get2DMaskedImage(InputImageType::Pointer im, LabelImageType::Pointer l);
vtkSmartPointer<vtkVolume> getOneVTKVolume(vtkSmartPointer<vtkImageData> vtkim, float colors[3]);
void getFeatureVectorsFarsight(LabelImageType::Pointer im, InputImageType::Pointer in_image, std::vector<ftk::LabelImageFeatures> & feature_vector, int time, int tag);
// for efficiency they've been implemented as different functions. Function overloading could be done.

bool deleteLabels1(LabelImageType::Pointer im, int n1);
bool deleteLabels2(LabelImageType::Pointer im, int n1, int n2);
bool deleteLabels3(LabelImageType::Pointer im, int n1, int n2, int n3);
bool deleteLabels4(LabelImageType::Pointer im, int n1, int n2, int n3, int n4);
bool deleteLabels5(LabelImageType::Pointer im, int n1, int n2, int n3, int n4, int n5);

bool mergeLabels1(LabelImageType::Pointer im, int merge_to,int n1);
bool mergeLabels2(LabelImageType::Pointer im, int merge_to,int n1, int n2);
bool mergeLabels3(LabelImageType::Pointer im, int merge_to,int n1, int n2, int n3);
bool mergeLabels4(LabelImageType::Pointer im, int merge_to,int n1, int n2, int n3, int n4);
bool mergeLabels5(LabelImageType::Pointer im, int merge_to,int n1, int n2, int n3, int n4, int n5);

#endif
