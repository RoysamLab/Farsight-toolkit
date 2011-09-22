
#include <iostream>

#include "itkTimeProbe.h"

#include "Common.h"

/**
 * This class implements a vessel tracer using a spherical structuring element for tracing the vasuclature 
 * from 3D microscopy images. This is a C++ implementation of Amit's vessel tracing work. 
 *
 * Major steps of the algorithm are:
 * 
 * Pre-processing:
 *  1. Calculate gradient vector field
 *  2. Median filter (and flat-field correction?) of the input data
 *  3. Save images as *.mhd files
 * 
 * Part 1:
 *	1. Generate primary nodes (initialization using a grid) - using energy minimization for the spherical model
 *	2. Generate secondary nodes - using ODFs and additional constraints
 *	Result: Final set of nodes
 *
 * Part 2:
 *	1. Generate affinity graph
 *	2. Create MST from graph
 *	3. Prune the tree - based on heuristics
 * 
 * Post processing:
 * 1. Complete loops and generate a vascular mask (using a cylinder model)
 * Result: Final tracing of vasculature 
 */

class ftkVesselTracer{

public:

	ftkVesselTracer();
	ftkVesselTracer(std::string, bool);
	~ftkVesselTracer();

	/** Preprocessing of data: 1. Reading the data from TIFF files 2. Median filtering 3. Edge enhancing
	 * 4. GVF computations 5. Saving all data and GVF as MHD files
	 * (data path with extension, empty data pointer)
	 */
	int PreprocessData(std::string, ImageType3D::Pointer&);
	
	/** Load gx, gy, gz, data and oribin for further processing.
	 * (data path)
	 */
	void LoadPreprocessedData(std::string);

	/* Returns maximum intensity projection image
	 * (ITK image ptr)
	 */
	void ComputeMaximumIntensityProjectionImage(void);

	/* Render MIP image
	 * (ITK data ptr)
	 */
	void RenderMIPImage(void);

	typedef float PixelType;
	typedef itk::Image<PixelType, 3> ImageType3D;
	typedef unsigned char RenderPixelType;
	typedef itk::Image<RenderPixelType, 3> RenderImageType3D;
	
	typedef itk::ImageToVTKImageFilter<RenderImageType3D> ITKToVTKConnectorType;
	
private:
	
	ImageType3D::Pointer originalData, inputData, gx, gy, gz;
	RenderImageType3D::Pointer originalDataForRendering, inputDataForRendering, MIPImage;

};