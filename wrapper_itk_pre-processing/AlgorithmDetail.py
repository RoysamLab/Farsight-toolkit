# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# AlgorithmDetail.py
# This file contains the details of the algorithm execution and implementation. 
# The main reason was to move the implementation away from the main script.

# Author 	: Adarsh K. Ramasubramonian
# Date 		: 22 May, 2009
#
# Things to do.
# 1. Check if IO model can be added to Input/Output Image Stack 
# 2. Include useCaster  - DONE 6/3/09.

# Last modified on 4 June, 2009 - GetReader in Resample - only had two arguments,
# added one more.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

from InsightToolkit import *	# Insight Toolkit
import sys
import cmath

ITK = "itk"
# Reader
READER = "itkImageFileReader"			# Class Name
READER_PIXELTYPE = "F"					# Input Image pixeltype (UC, US, F)
# READER_DIM							# Input Image dimension (2, 3)
# Writer
WRITER = "itkImageFileWriter"			# Class Name
# WRITER_TYPE							# Output Image type (UC, US, F)
# WRITER_DIM							# Output Image dimension (2, 3)


# Default execution for all algorithms (until 6/3/09) except resampling filter.

def ExecuteCommands(thisalgorithm):
	""" Default execution of commands """	

	dimension 	= thisalgorithm.GetDimension()
	pixelType 	= thisalgorithm.GetOutputImagePixelType()

	reader = GetReader(thisalgorithm, dimension, pixelType)			# Reader Object
	writer = GetWriter(thisalgorithm, dimension, pixelType)			# Writer Object
	filter = GetFilter(thisalgorithm, dimension, pixelType)			# Filter Object

	parameters = thisalgorithm.GetParameters()
	object = {}
	keys = parameters.keys()[:]
	if thisalgorithm.GetHasStructuringElement():
		keys.remove("Radius")
		radius = eval( str( parameters["Radius"]["value"] ) )
		AssignStructuringElement(filter, READER_PIXELTYPE, dimension, radius)

	for key in keys:
		if key == "Order":
			orderVar = eval( parameters[key]["value"] )
			AssignOrderOfFilter(filter, orderVar)

		else:
# Added on 19 May, 2009
			if thisalgorithm.IsTupleParameter(key):
				value = parameters[key]["value"]
				value = eval(str(value))	
				# Create the object				
				tempObject = eval(parameters[key]["type"] + dimension  + "()" )
				for num in range(len(value)):
					tempObject.SetElement(num,  value[num])	
				eval( "filter.Set" + key + "( tempObject )" )
			else:
				eval( "filter.Set" + key + "( " + str(parameters[key]["value"]) + " )" )

	filter.SetInput ( reader.GetOutput() )
	
	if thisalgorithm.GetUseCaster():
		caster = eval ( "itkRescaleIntensityImageFilter" + READER_PIXELTYPE + dimension + pixelType + dimension + "_New()" )
	
		caster.SetOutputMinimum(   0 );
		if pixelType == "UC":
			caster.SetOutputMaximum( 255 );
		elif pixelType == "US":
			caster.SetOutputMaximum( 65535 );

		caster.SetInput ( filter.GetOutput() )
		writer.SetInput ( caster.GetOutput() )
	else:
		writer.SetInput ( filter.GetOutput() )

	writer.Update()

	print "Algorithm successfully implemented:\n"
	print thisalgorithm

def GetReader(thisalgorithm, dimension, pixelType):
	""" Get Reader Object to read the input file """
	# Check dimension
	# Check if input stack or single image
	# If input stack
	if dimension == "3" and thisalgorithm.GetInputImageStack():
		READER = "itkImageSeriesReader"
		READER_PIXELTYPE = "F"

	#	ImageSeriesReader			
		if thisalgorithm.GetUseCaster():
			reader = eval( READER + READER_PIXELTYPE + dimension + "_New()" )
		else:
			reader = eval( READER + pixelType + dimension + "_New()" )
	# 	Name generator
		readerNameGenerator = itkNumericSeriesFileNames_New()
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(readerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1) doesn't seem to be working
	# 	Set file name
		reader.SetFileNames( readerNameGenerator.GetFileNames() )

	elif dimension == "2" or dimension == "3":
	# If single image
	# 	ImageFileReader
		READER = "itkImageFileReader"
		READER_PIXELTYPE = "F"
		if thisalgorithm.GetUseCaster():
			reader = eval( READER + READER_PIXELTYPE + dimension + "_New()" ) 
		else:
			reader = eval( READER + pixelType + dimension + "_New()" )
	#   Set FileName
		reader.SetFileName( thisalgorithm.GetInputFileName() )

	return reader


def GetWriter(thisalgorithm, dimension, pixelType):
	""" Get Writer Object to write the output file """
	# Check dimension
	# Check if input stack or single image
	# If input stack
	if dimension == "3" and thisalgorithm.GetOutputImageStack():
		WRITER = "itkImageSeriesWriter"
		WRITER_PIXELTYPE = pixelType

	#	ImageSeriesReader			
		writer = eval( WRITER + pixelType + dimension + pixelType + "2" + "_New()" )
	# 	Name generator
		writerNameGenerator = itkNumericSeriesFileNames_New()
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(writerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1)  doesn't seem to be working

	# 	Set file name
		writer.SetFileNames( writerNameGenerator.GetFileNames() )

	elif dimension == "2" or dimension == "3":
	# If single image
	# 	ImageFileReader
		WRITER = "itkImageFileWriter"
		WRITER_PIXELTYPE = pixelType
		writer = eval( WRITER + pixelType + dimension + "_New()" ) 
	#   Set FileName
		writer.SetFileName( thisalgorithm.GetOutputFileName() )
#		print pixelType, " is the pxiel type"

	return writer

def GetFilter(thisalgorithm, dimension, pixelType):	
	""" Get the filter object - the main algorithm """
	if thisalgorithm.GetOutputType():
		if thisalgorithm.GetUseCaster():
			filter = eval ( thisalgorithm.GetName() + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )
		else:
			filter = eval ( thisalgorithm.GetName() + pixelType + dimension + pixelType + dimension + "_New()" )
	else:
		filter = eval ( thisalgorithm.GetName() + READER_PIXELTYPE + dimension + "_New()" )
	return filter

def AssignOrderOfFilter(filter, order):
	""" Order of filter - required for Gaussian filters """
	if order == 0:
		filter.SetZeroOrder()
	elif order == 1:
		filter.SetFirstOrder()
	elif order== 2:
		filter.SetSecondOrder()
	else:
		print "Please check order of filter !"
		sys.exit()

# Required for Grayscale Erode and Dilate Filter.
def AssignStructuringElement(filter, READER_PIXELTYPE, dimension, radius):
	element = eval( "itkBinaryBallStructuringElement" + READER_PIXELTYPE + dimension + "()" )
	element.SetRadius( radius )
	element.CreateStructuringElement()
	filter.SetKernel(element)

def SetSeriesAndIndices(Generator, algorithm):
	Generator.SetSeriesFormat( algorithm.GetInputFileName() )
	Generator.SetStartIndex( algorithm.GetStartIndex() )
	Generator.SetEndIndex( algorithm.GetEndIndex() )
	Generator.SetIncrementIndex( algorithm.GetIncrementIndex() )

# Separate Execution strategy for Resample filters.
def AlgorithmDetailResampleFilter(thisalgorithm):

	dimension = thisalgorithm.GetDimension()
	pixelType = thisalgorithm.GetOutputImagePixelType()

	reader = GetReader(thisalgorithm, dimension, pixelType)
	writer = GetWriter(thisalgorithm, dimension, pixelType)

	# Define the smoothing filters.
	smootherX = eval( "itkRecursiveGaussianImageFilter" + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )
	smootherY = eval( "itkRecursiveGaussianImageFilter" + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )

	smootherX.SetInput( reader.GetOutput() )
	smootherY.SetInput( smootherX.GetOutput() )

	# Assign the parameters of the filters
	inputImage = reader.GetOutput()
	inputSpacing = inputImage.GetSpacing()

	print inputSpacing.GetElement(0), " and ", inputSpacing.GetElement(1), " and ", inputSpacing.GetElement(2)
	isoSpacing = float( abs( cmath.sqrt( inputSpacing.GetElement(2) * inputSpacing.GetElement(0) ) ) )
	smootherX.SetSigma( isoSpacing )
	smootherY.SetSigma( isoSpacing )

	smootherX.SetDirection( 0 )
	smootherY.SetDirection( 1 )

	smootherX.SetNormalizeAcrossScale( True )
	smootherY.SetNormalizeAcrossScale( True )

	# Define the filter - ResampleImageFilter
	filter = eval( thisalgorithm.GetName() + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )

	# Set transform of resampler to identity
	transform = itkIdentityTransform3_New()
	transform.SetIdentity()
	filter.SetTransform( transform )

	# Define and set interpolator
	interpolator = eval( itkLinearInterpolateImageFunction + READER_PIXELTYPE + dimension + "D_New()" )
	filter.SetInterpolator( interpolator )

	filter.SetDefaultPixelValue( 255 )		# Default pixel value

	spacing[0] = isoSpacing 
	spacing[1] = isoSpacing 
	spacing[2] = isoSpacing 

	filter.SetOutputSpacing( spacing )
	print filter.GetOutputSpacing( spacing )

	filter.SetOutputOrigin( inputImage.GetOrigin() )
	filter.SetOutputDirection( inputImage.GetDirection() )

	inputSize = inputImage.GetLargestPossibleRegion().GetSize()

	dx = ( inputSize[0] * inputSpacing[0] / isoSpacing )
	dy = ( inputSize[1] * inputSpacing[1] / isoSpacing )
	dz = ( (inputSize[2] - 1) * inputSpacing[2] / isoSpacing )
	size = [dx, dy, dz]
	filter.SetSize(size)
	filter.SetInput ( smootherY.GetOutput() )
	writer.SetInput ( filter.GetOutput() )
	writer.Update()


