// Create a grid

#include "GridActors.h"

GridlineActors::GridlineActors()
{
	/*!
	 * Default constructor
	 */
	num_horizontal_lines = 0;
	num_vertical_lines = 0;
	num_depth_lines = 0;
	GridlineActorVectorHorizontal = 0;
	GridlineActorVectorVertical = 0;
	GridlineActorVectorDepth = 0;
}
GridlineActors::~GridlineActors()
{
	/*!
	 * Destructor
	 */
	if(GridlineActorVectorHorizontal != 0)
	{
		delete [] GridlineActorVectorHorizontal;
	}
	if(GridlineActorVectorVertical != 0)
	{
		delete [] GridlineActorVectorVertical;
	}
	if(GridlineActorVectorDepth != 0)
	{
		delete [] GridlineActorVectorDepth;
	}
}
//void GridlineActors::createGrid2D(double left_bound, double right_bound, double bottom_bound, double top_bound,int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value)
//{
//	/*!
//	 * Create a 2D grid
//	 */
//	num_depth_lines = 0;
//
//	double r_color = r/256.0;
//	double g_color = g/256.0;
//	double b_color = b/256.0;
//	double line_opacity = opacity/100.0;
//
//	vtkSmartPointer<vtkProperty> lineproperty = vtkSmartPointer<vtkProperty>::New();
//	lineproperty->SetColor(r_color,g_color,b_color);
//	lineproperty->SetOpacity(line_opacity);
//	lineproperty->SetLineWidth(line_width);
//
//	//Manually create lines
//	/// horizontal lines
//	num_horizontal_lines = (top_bound - bottom_bound) / height_spacing + 1;
//	GridlineActorVectorHorizontal = new vtkSmartPointer<vtkActor>[num_horizontal_lines];
//	
//	unsigned int horizontal_line_index = 0;
//	for (double increment = bottom_bound; increment < top_bound; increment+=height_spacing, horizontal_line_index++)
//	{	
//		// Create two points, P0 and P1
//		double p0[3] = {left_bound, increment, z_plane_value};
//		double p1[3] = {right_bound, increment, z_plane_value};
//
//		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
//		lineSource->SetPoint1(p0);
//		lineSource->SetPoint2(p1);
//		lineSource->Update();
//
//		// Visualize
//		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//		mapper->SetInputConnection(lineSource->GetOutputPort());
//		//Create GridlineActor
//		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
//		GridlineActor->SetProperty(lineproperty);
//		//Store it to the GridlineActorVector
//		GridlineActorVectorHorizontal[horizontal_line_index] = GridlineActor;
//		
//		//Add the GridlineActor to the renderer
//		GridlineActor->SetMapper(mapper);
//		GridlineActor->SetPickable(0);
//	}
//
//	/// vertical lines
//	num_vertical_lines = (right_bound - left_bound) / width_spacing + 1;
//	GridlineActorVectorVertical = new vtkSmartPointer<vtkActor>[num_vertical_lines];
//
//	unsigned int vertical_line_index = 0;
//	for (double increment = left_bound; increment < right_bound; increment+=width_spacing, vertical_line_index++)
//	{	
//		// Create two points, P0 and P1
//		double p0[3] = {increment, bottom_bound, z_plane_value};
//		double p1[3] = {increment, top_bound, z_plane_value};
//
//		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
//		lineSource->SetPoint1(p0);
//		lineSource->SetPoint2(p1);
//		lineSource->Update();
//
//		// Visualize
//		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//		mapper->SetInputConnection(lineSource->GetOutputPort());
//		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
//		GridlineActor->SetProperty(lineproperty);
//		GridlineActorVectorVertical[vertical_line_index] = GridlineActor;
//		GridlineActor->SetMapper(mapper);
//		GridlineActor->SetPickable(0);
//	}
//}
void GridlineActors::createGridxy(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value)
{
	/*!
	 * Create a 2D xy grid
	 */
	num_depth_lines = 0;

	double r_color = r/256.0;
	double g_color = g/256.0;
	double b_color = b/256.0;
	double line_opacity = opacity/100.0;

	vtkSmartPointer<vtkProperty> lineproperty = vtkSmartPointer<vtkProperty>::New();
	lineproperty->SetColor(r_color,g_color,b_color);
	lineproperty->SetOpacity(line_opacity);
	lineproperty->SetLineWidth(line_width);

	//Manually create lines
	/// horizontal lines
	num_horizontal_lines = (bounds[3] - bounds[2]) / height_spacing + 1;
	if(GridlineActorVectorHorizontal != 0)
	{
		delete [] GridlineActorVectorHorizontal;
	}
	GridlineActorVectorHorizontal = new vtkSmartPointer<vtkActor>[num_horizontal_lines];
	
	unsigned int horizontal_line_index = 0;
	for (double increment = bounds[2]; increment < bounds[3]; increment+=height_spacing, horizontal_line_index++)
	{	
		// Create two points, P0 and P1
		double p0[3] = {bounds[0], increment, z_plane_value};
		double p1[3] = {bounds[1], increment, z_plane_value};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		//Create GridlineActor
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		//Store it to the GridlineActorVector
		GridlineActorVectorHorizontal[horizontal_line_index] = GridlineActor;
		
		//Add the GridlineActor to the renderer
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}

	/// vertical lines
	num_vertical_lines = (bounds[1] - bounds[0]) / width_spacing + 1;
	if(GridlineActorVectorVertical != 0)
	{
		delete [] GridlineActorVectorVertical;
	}
	GridlineActorVectorVertical = new vtkSmartPointer<vtkActor>[num_vertical_lines];

	unsigned int vertical_line_index = 0;
	for (double increment = bounds[0]; increment < bounds[1]; increment+=width_spacing, vertical_line_index++)
	{	
		// Create two points, P0 and P1
		double p0[3] = {increment, bounds[2], z_plane_value};
		double p1[3] = {increment, bounds[3], z_plane_value};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		GridlineActorVectorVertical[vertical_line_index] = GridlineActor;
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}
}
void GridlineActors::createGridxz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value)
{
	/*!
	 * Create a 2D xz grid
	 */
	num_depth_lines = 0;
	double r_color = r/256.0;
	double g_color = g/256.0;
	double b_color = b/256.0;
	double line_opacity = opacity/100.0;

	vtkSmartPointer<vtkProperty> lineproperty = vtkSmartPointer<vtkProperty>::New();
	lineproperty->SetColor(r_color,g_color,b_color);
	lineproperty->SetOpacity(line_opacity);
	lineproperty->SetLineWidth(line_width);

	//Manually create lines
	/// horizontal lines
	num_horizontal_lines = (bounds[5] - bounds[4]) / height_spacing + 1;
	if(GridlineActorVectorHorizontal != 0)
	{
		delete [] GridlineActorVectorHorizontal;
	}
	GridlineActorVectorHorizontal = new vtkSmartPointer<vtkActor>[num_horizontal_lines];
	
	unsigned int horizontal_line_index = 0;
	for (double increment = bounds[4]; increment < bounds[5]; increment+=height_spacing, horizontal_line_index++)
	{	
		// Create two points, P0 and P1
		double p0[3] = {bounds[0], 0.0, increment};
		double p1[3] = {bounds[1], 0.0, increment};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		//Create GridlineActor
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		//Store it to the GridlineActorVector
		GridlineActorVectorHorizontal[horizontal_line_index] = GridlineActor;
		
		//Add the GridlineActor to the renderer
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}

	/// vertical lines
	num_vertical_lines = (bounds[1] - bounds[0]) / width_spacing + 1;
	if(GridlineActorVectorVertical != 0)
	{
		delete [] GridlineActorVectorVertical;
	}
	GridlineActorVectorVertical = new vtkSmartPointer<vtkActor>[num_vertical_lines];

	unsigned int vertical_line_index = 0;
	for (double increment = bounds[0]; increment < bounds[1]; increment+=width_spacing, vertical_line_index++)
	{
		// Create two points, P0 and P1
		double p0[3] = {increment, 0.0, bounds[4]};
		double p1[3] = {increment, 0.0, bounds[5]};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		GridlineActorVectorVertical[vertical_line_index] = GridlineActor;
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}
}
void GridlineActors::createGridyz(double bounds[],int height_spacing, int width_spacing, int line_width, int r, int g, int b, int opacity, int z_plane_value)
{
	/*!
	 * Create a 2D yz grid
	 */
	num_depth_lines = 0;
	double r_color = r/256.0;
	double g_color = g/256.0;
	double b_color = b/256.0;
	double line_opacity = opacity/100.0;

	vtkSmartPointer<vtkProperty> lineproperty = vtkSmartPointer<vtkProperty>::New();
	lineproperty->SetColor(r_color,g_color,b_color);
	lineproperty->SetOpacity(line_opacity);
	lineproperty->SetLineWidth(line_width);

	//Manually create lines
	/// horizontal lines
	num_horizontal_lines = (bounds[5] - bounds[4]) / height_spacing + 1;
	if(GridlineActorVectorHorizontal != 0)
	{
		delete [] GridlineActorVectorHorizontal;
	}
	GridlineActorVectorHorizontal = new vtkSmartPointer<vtkActor>[num_horizontal_lines];
	
	unsigned int horizontal_line_index = 0;
	for (double increment = bounds[4]; increment < bounds[5]; increment+=height_spacing, horizontal_line_index++)
	{	
		// Create two points, P0 and P1
		double p0[3] = {0.0, bounds[2], increment};
		double p1[3] = {0.0, bounds[3], increment};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		//Create GridlineActor
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		//Store it to the GridlineActorVector
		GridlineActorVectorHorizontal[horizontal_line_index] = GridlineActor;
		
		//Add the GridlineActor to the renderer
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}

	/// vertical lines
	num_vertical_lines = (bounds[3] - bounds[2]) / width_spacing + 1;
	if(GridlineActorVectorVertical != 0)
	{
		delete [] GridlineActorVectorVertical;
	}
	GridlineActorVectorVertical = new vtkSmartPointer<vtkActor>[num_vertical_lines];

	unsigned int vertical_line_index = 0;
	for (double increment = bounds[2]; increment < bounds[3]; increment+=width_spacing, vertical_line_index++)
	{	
		// Create two points, P0 and P1
		double p0[3] = {0.0, increment, bounds[4]};
		double p1[3] = {0.0, increment, bounds[5]};

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		lineSource->SetPoint1(p0);
		lineSource->SetPoint2(p1);
		lineSource->Update();

		// Visualize
		vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection(lineSource->GetOutputPort());
		vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
		GridlineActor->SetProperty(lineproperty);
		GridlineActorVectorVertical[vertical_line_index] = GridlineActor;
		GridlineActor->SetMapper(mapper);
		GridlineActor->SetPickable(0);
	}
}
void GridlineActors::createGrid3D(double bounds[],int height_spacing, int width_spacing, int depth_spacing, int line_width, int r, int g, int b, int opacity)
{
	/*!
	 * Create a 3D grid.
	 * @param bounds[] image bounds.
	 */
	double r_color = r/256.0;
	double g_color = g/256.0;
	double b_color = b/256.0;
	double line_opacity = opacity/100.0;

	vtkSmartPointer<vtkProperty> lineproperty = vtkSmartPointer<vtkProperty>::New();
	lineproperty->SetColor(r_color,g_color,b_color);
	lineproperty->SetOpacity(line_opacity);
	lineproperty->SetLineWidth(line_width);

	//Manually create lines

	int temp_depth_lines = (bounds[5] - bounds[4]) / depth_spacing + 1;
	int temp_horizontal_lines = (bounds[3] - bounds[2]) / height_spacing + 1;
	int temp_vertical_lines = (bounds[1] - bounds[0]) / width_spacing + 1;

	/// horizontal lines
	num_horizontal_lines = temp_horizontal_lines*temp_depth_lines;
	if(GridlineActorVectorHorizontal != 0)
	{
		delete [] GridlineActorVectorHorizontal;
	}
	GridlineActorVectorHorizontal = new vtkSmartPointer<vtkActor>[num_horizontal_lines];

	/// vertical lines
	num_vertical_lines = temp_vertical_lines*temp_depth_lines;
	if(GridlineActorVectorVertical != 0)
	{
		delete [] GridlineActorVectorVertical;
	}
	GridlineActorVectorVertical = new vtkSmartPointer<vtkActor>[num_vertical_lines];

	/// depth lines
	num_depth_lines = temp_vertical_lines*temp_horizontal_lines;
	if(GridlineActorVectorDepth != 0)
	{
		delete [] GridlineActorVectorDepth;
	}
	GridlineActorVectorDepth = new vtkSmartPointer<vtkActor>[num_depth_lines];

	unsigned int horizontal_line_index = 0;
	unsigned int vertical_line_index = 0;
	unsigned int depth_line_index = 0;
	for (double depthIncrement = bounds[4]; depthIncrement < bounds[5]; depthIncrement+= depth_spacing)
	{
		// horizontal lines
		for (double increment = bounds[2]; increment < bounds[3]; increment+=height_spacing, horizontal_line_index++)
		{
			// Create two points, P0 and P1
			double p0[3] = {bounds[0], increment, depthIncrement};
			double p1[3] = {bounds[1], increment, depthIncrement};

			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			lineSource->SetPoint1(p0);
			lineSource->SetPoint2(p1);
			lineSource->Update();

			// Visualize
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(lineSource->GetOutputPort());
			//Create GridlineActor
			vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
			GridlineActor->SetProperty(lineproperty);
			//Store it to the GridlineActorVector
			GridlineActorVectorHorizontal[horizontal_line_index] = GridlineActor;
			
			//Add the GridlineActor to the renderer
			GridlineActor->SetMapper(mapper);
			GridlineActor->SetPickable(0);
		}

		// vertical lines
		for (double increment = bounds[0]; increment < bounds[1]; increment+=width_spacing, vertical_line_index++)
		{	
			// Create two points, P0 and P1
			double p0[3] = {increment, bounds[2], depthIncrement};
			double p1[3] = {increment, bounds[3], depthIncrement};

			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			lineSource->SetPoint1(p0);
			lineSource->SetPoint2(p1);
			lineSource->Update();

			// Visualize
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(lineSource->GetOutputPort());
			vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
			GridlineActor->SetProperty(lineproperty);
			GridlineActorVectorVertical[vertical_line_index] = GridlineActor;
			GridlineActor->SetMapper(mapper);
			GridlineActor->SetPickable(0);			
		}
	}
	for (double y_increment = bounds[2]; y_increment < bounds[3]; y_increment+=height_spacing)
	{
		for (double x_increment = bounds[0]; x_increment < bounds[1]; x_increment+=width_spacing, depth_line_index++)
		{
			//depth lines
			// Create two points, P0 and P1
			double p0[3] = {x_increment, y_increment, bounds[4]};
			double p1[3] = {x_increment, y_increment, bounds[5]};

			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			lineSource->SetPoint1(p0);
			lineSource->SetPoint2(p1);
			lineSource->Update();

			// Visualize
			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(lineSource->GetOutputPort());
			vtkSmartPointer<vtkActor> GridlineActor = vtkSmartPointer<vtkActor>::New();
			GridlineActor->SetProperty(lineproperty);
			GridlineActorVectorDepth[depth_line_index] = GridlineActor;
			
			//Add the GridlineActor to the renderer
			GridlineActor->SetMapper(mapper);
			GridlineActor->SetPickable(0);
		}
	}
}
vtkSmartPointer<vtkActor> GridlineActors::GetHorizontalGridlines(int i)
{
	/*!
	 * @return horizontal gridlines
	 */
	if (i == -1)
	{
		return NULL;
	}
	return this->GridlineActorVectorHorizontal[i];
}
vtkSmartPointer<vtkActor> GridlineActors::GetVerticalGridlines(int i)
{
	/*!
	 * @return vertical gridlines
	 */
	if (i == -1)
	{
		return NULL;
	}
	return this->GridlineActorVectorVertical[i];
}
vtkSmartPointer<vtkActor> GridlineActors::GetDepthGridlines(int i)
{
	/*!
	 * @return depth gridlines
	 */
	if (i == -1)
	{
		return NULL;
	}
	return this->GridlineActorVectorDepth[i];
}
int GridlineActors::NumberOfLines()
{
	/*!
	 * @return number of gridlines
	 */
	return (num_horizontal_lines + num_vertical_lines + num_depth_lines);
}
int GridlineActors::NumberOfHorizontalLines()
{
	/*!
	 * @return number of horizontal gridlines
	 */
	return num_horizontal_lines;
}
int GridlineActors::NumberOfVerticalLines()
{
	/*!
	 * @return number of vertical gridlines
	 */
	return num_vertical_lines;
}
int GridlineActors::NumberOfDepthLines()
{
	/*!
	 * @return number of depth gridlines
	 */
	return num_depth_lines;
}
