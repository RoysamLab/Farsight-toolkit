/*=========================================================================

=========================================================================*/
#ifndef GRIDACTORS_H_
#define GRIDACTORS_H_

#include "ImageActors.h"

#include "vtkActor.h"
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

class GridlineActors
{

public:
	GridlineActors();
	~GridlineActors();
	void createGrid(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGridxz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	void createGridyz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value);
	vtkSmartPointer<vtkActor> GetHorizontalGridlines(int i);
	vtkSmartPointer<vtkActor> GetVerticalGridlines(int i);
	int NumberOfLines();
	int NumberOfHorizontalLines();
	int NumberOfVerticalLines();

private:
	vtkSmartPointer<vtkActor>* GridlineActorVectorHorizontal;
	vtkSmartPointer<vtkActor>* GridlineActorVectorVertical;

	int num_horizontal_lines;
	int num_vertical_lines;
};
#endif
