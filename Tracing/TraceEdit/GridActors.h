#ifndef GRIDACTORS_H_
#define GRIDACTORS_H_

#include "ImageActors.h"

#include "vtkActor.h"
#include "vtkFloatArray.h"
#include "vtkImageActor.h"
#include "vtkImageCast.h"
#include "vtkImageMapToColors.h"
#include "vtkImageGridSource.h"
#include "vtkLineSource.h"
#include "vtkLookupTable.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"

/**
 * @author Audrey Cheong
 */

/*! Gridlines class */
class GridlineActors
{
public:
	GridlineActors();
	~GridlineActors();
	//void createGrid2D(double left_bound, double right_bound, double bottom_bound, double top_bound,int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGridxy(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGridxz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGridyz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGrid3D(double bounds[],int height_spacing, int width_spacing, int depth_spacing, int line_width, int r, int g, int b, int opacity);
	vtkSmartPointer<vtkActor> GetHorizontalGridlines(int i);
	vtkSmartPointer<vtkActor> GetVerticalGridlines(int i);
	vtkSmartPointer<vtkActor> GetDepthGridlines(int i);
	int NumberOfLines();
	int NumberOfHorizontalLines();
	int NumberOfVerticalLines();
	int NumberOfDepthLines();

private:
	vtkSmartPointer<vtkActor>* GridlineActorVectorHorizontal;	/*!< Contains horizontal gridline actors */
	vtkSmartPointer<vtkActor>* GridlineActorVectorVertical;		/*!< Contains vertical gridline actors */
	vtkSmartPointer<vtkActor>* GridlineActorVectorDepth;		/*!< Contains depth gridline actors */

	int num_horizontal_lines;	/*!< Number of horizontal lines */
	int num_vertical_lines;		/*!< Number of vertical lines */
	int num_depth_lines;		/*!< Number of depth lines */
};
#endif
