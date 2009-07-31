# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# AlgorithmDetail.py
# This file contains the details of the algorithm execution and implementation. 
# The main reason was to move the implementation away from the main script. The use 
# of wrappers is restricted to this module, and if we need to move from CSwig to 
# WrapITK, we only need to change this file. 

# Author 	: Adarsh K. Ramasubramonian
# Date 		: 22 May, 2009
#
# Things to do.
# 1. Check if IO model can be added to Input/Output Image Stack - To be done XXX
# 2. Include useCaster  - DONE 6/3/09.

# Last modified on 4 June, 2009 - GetReader in Resample - only had two arguments,
# added one more.
#
# To be done: Vesselness - wrappers are still not there.
#             
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

from InsightToolkit import *	# Insight Toolkit
import sys
import cmath
import pdb
import FilterObjectXML
import FilterAlgorithm
import warnings

ITK = "itk"
# Reader
READER = "itkImageFileReader"			# Class Name
READER_PIXELTYPE = "F"					# Input Image pixeltype (UC, US, F)
                                        # We choose this by default, and use a resampler to change 
                                        # pixel type back to unsigned char or unsigned int.
# READER_DIM							# Input Image dimension (2, 3)
# Writer    
WRITER = "itkImageFileWriter"			# Class Name
# WRITER_TYPE							# Output Image type (UC, US, F)
# WRITER_DIM							# Output Image dimension (2, 3)


def defaultAlgorithmExecution(XMLFileName):
    """ Details of algorithm in XML file. Read it and execute the algorithm """
    warnings.warn("Some of the algorithms for Python may not be supported in CSwig (without additional wrappers). Please use WrapITK.")
    xmlObject = FilterObjectXML.FilterObjectXML()
    thisalgorithm = FilterAlgorithm.FilterAlgorithm()
    xmlObject.ParseDocument(thisalgorithm, XMLFileName)

    if thisalgorithm.GetKey() == "resample":	# Special steps for Resampling
        ExecuteResampleAlgorithm(thisalgorithm)
    elif thisalgorithm.GetKey() == "anisotropicdiffusionvesselenhancement":
        print "Anisotropic Diffusion Vessel Enhancement algorithm not available for CSwig. Sorry !!!"
        sys.exit(1)
#        ExecuteVesselnessAlgorithm(thisalgorithm)
    else:
        ExecuteCommands(thisalgorithm)			# Default algorithm execution



# Default execution for all algorithms (until 6/3/09) except resampling filter.

def ExecuteCommands(thisalgorithm):
    """ Default execution of commands """   

    dimension   = thisalgorithm.GetDimension()
    pixelType   = thisalgorithm.GetOutputImagePixelType()

    reader = GetReader(thisalgorithm, dimension, pixelType)         # Reader Object
    writer = GetWriter(thisalgorithm, dimension, pixelType)         # Writer Object
    filter = GetFilter(thisalgorithm, dimension, pixelType)         # Filter Object

    parameters = thisalgorithm.GetParameters()
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
                SetElementFromList(tempObject, value, len(value))
#               for num in range(len(value)):
#                   tempObject.SetElement(num,  value[num]) 
                eval( "filter.Set" + key + "( tempObject )" )
            else:
                eval( "filter.Set" + key + "( " + str(parameters[key]["value"]) + " )" )

    
    if not thisalgorithm.GetUseRescaler():
        if pixelType != "F":
            warnings.warn("\nWARNING: The output may lose some details if not stored as float type !\n")

    rescaler = GetRescaler( READER_PIXELTYPE, dimension, pixelType )
    PipelineInputOutput([reader, filter, rescaler, writer])
        
    writer.Update()

    print thisalgorithm
    print "Algorithm successfully executed!\n"

def GetReader(thisalgorithm, dimension, pixelType):
	""" Get Reader Object to read the input file """
	# Check dimension
	# Check if input stack or single image
	# If input stack
	if dimension == "3" and thisalgorithm.GetInputImageStack():
		READER = "itkImageSeriesReader"
		READER_PIXELTYPE = "F"

															#	ImageSeriesReader			
		if thisalgorithm.GetUseRescaler():
			reader = eval( READER + READER_PIXELTYPE + dimension + "_New()" )
		else:
			reader = eval( READER + pixelType + dimension + "_New()" )
		readerNameGenerator = itkNumericSeriesFileNames_New()			# Name generator
		readerNameGenerator.SetSeriesFormat( thisalgorithm.GetInputFileName() )														
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(readerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1) doesn't seem to be working
		reader.SetFileNames( readerNameGenerator.GetFileNames() )		# Set Filename

	elif dimension == "2" or dimension == "3":
	# If single image
		READER = "itkImageFileReader"
		READER_PIXELTYPE = "F"
		if thisalgorithm.GetUseRescaler():
			reader = eval( READER + READER_PIXELTYPE + dimension + "_New()" ) 
		else:
			reader = eval( READER + pixelType + dimension + "_New()" )
		reader.SetFileName( thisalgorithm.GetInputFileName() )			# Set Filename

	return reader


def GetWriter(thisalgorithm, dimension, pixelType):
	""" Get Writer Object to write the output file """
	# If input stack
	if dimension == "3" and thisalgorithm.GetOutputImageStack():

		WRITER = "itkImageSeriesWriter"
		WRITER_PIXELTYPE = pixelType

		writer = eval( WRITER + pixelType + dimension + pixelType + "2" + "_New()" )	# ImageSeriesReader			
		writerNameGenerator = itkNumericSeriesFileNames_New()							# Name generator
		writerNameGenerator.SetSeriesFormat( thisalgorithm.GetOutputFileName() )
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(writerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1)  doesn't seem to be working

		writer.SetFileNames( writerNameGenerator.GetFileNames() ) 		# Set file name

	elif dimension == "2" or dimension == "3":							# If single image
		WRITER = "itkImageFileWriter"
		WRITER_PIXELTYPE = pixelType
		writer = eval( WRITER + pixelType + dimension + "_New()" ) 
		writer.SetFileName( thisalgorithm.GetOutputFileName() )			# Set FileName

	return writer

def GetFilter(thisalgorithm, dimension, pixelType):	
	""" Get the filter object - the main algorithm """
	if thisalgorithm.GetOutputType():    
		if thisalgorithm.GetUseRescaler():
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
	Generator.SetStartIndex( algorithm.GetStartIndex() )
	Generator.SetEndIndex( algorithm.GetEndIndex() )
	Generator.SetIncrementIndex( algorithm.GetIncrementIndex() )

# Separate Execution strategy for Resample filters.
def ExecuteResampleAlgorithm(thisalgorithm):

	dimension = thisalgorithm.GetDimension()
	pixelType = thisalgorithm.GetOutputImagePixelType()

	reader = GetReader(thisalgorithm, dimension, pixelType)
	writer = GetWriter(thisalgorithm, dimension, pixelType)

	# Define the smoothing filters.
	smootherX = eval( "itkRecursiveGaussianImageFilter" + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )
	smootherY = eval( "itkRecursiveGaussianImageFilter" + READER_PIXELTYPE + dimension + READER_PIXELTYPE + dimension + "_New()" )


	# Assign the parameters of the filters
	inputImage = reader.GetOutput()
	inputImage.Update()

	inputSpacing = inputImage.GetSpacing()
	isoSpacing = float( abs( cmath.sqrt( inputSpacing.GetElement(2) * inputSpacing.GetElement(0) ) ) )
	smootherX.SetSigma( isoSpacing )
	smootherY.SetSigma( isoSpacing )

	smootherX.SetDirection( 0 )
	smootherY.SetDirection( 1 )

	smootherX.SetNormalizeAcrossScale( True )
	smootherY.SetNormalizeAcrossScale( True )

	# Define the filter - ResampleImageFilter
	filter = GetFilter(thisalgorithm, dimension, pixelType)

	# Set transform of resampler to identity
	transform = eval( "itkIdentityTransform" + dimension + "_New()" )
	transform.SetIdentity()
	filter.SetTransform( transform.GetPointer() )

	# Define and set interpolator
	interpolator = eval( "itkLinearInterpolateImageFunction" + READER_PIXELTYPE + dimension + "_New()" )
	filter.SetInterpolator( interpolator.GetPointer() )

	filter.SetDefaultPixelValue( 255 )		# Default pixel value

	spacing = filter.GetOutput().GetSpacing()		# Just to get the class type correct.
	SetElementFromList(spacing, [isoSpacing]*3, 3)

	filter.SetOutputSpacing( spacing )
	filter.SetOutputOrigin( inputImage.GetOrigin() )
	filter.SetOutputDirection( inputImage.GetDirection() )

	inputSize = inputImage.GetLargestPossibleRegion().GetSize()

	dx = ( inputSize.GetElement(0) * inputSpacing.GetElement(0) / isoSpacing )
	dy = ( inputSize.GetElement(1) * inputSpacing.GetElement(1) / isoSpacing )
	dz = ( (inputSize.GetElement(2) - 1) * inputSpacing.GetElement(2) / isoSpacing )

	size = eval( "itkSize" + dimension + "()" )
	SetElementFromList(size, [dx, dy, dz], 3)
	filter.SetSize(size)

	caster = GetRescaler(READER_PIXELTYPE, dimension, pixelType)

	PipelineInputOutput([reader, smootherX, smootherY, filter, caster, writer])

	writer.Update()
	print thisalgorithm
	print "Algorithm successfully executed!\n"

def GetRescaler( inputPixelType, dimension, outputPixelType):

	caster = eval ( "itkRescaleIntensityImageFilter" + inputPixelType + dimension + outputPixelType + dimension + "_New()" )
	
	caster.SetOutputMinimum(   0 )
	if outputPixelType == "UC":
		caster.SetOutputMaximum( 255 )
	elif outputPixelType == "US":
		caster.SetOutputMaximum( 65535 )
	return caster

def PipelineInputOutput(listOfObjects):
	for i in range(len(listOfObjects)-1):
		listOfObjects[i+1].SetInput( listOfObjects[i].GetOutput() )

#def ExecuteVesselnessAlgorithm(thisalgorithm):
#
#	pdb.set_trace()
#	dimension 	= thisalgorithm.GetDimension()
#	pixelType 	= thisalgorithm.GetOutputImagePixelType()
#
#	reader = GetReader(thisalgorithm, dimension, pixelType)
#	writer = GetWriter(thisalgorithm, dimension, pixelType)
#
#	parameters = thisalgorithm.GetParameters()
#	object = {}
#	keys = parameters.keys()[:]
#
#	hessianfilter = itkHessianRecursiveGaussianImageFilterF3_New()
#	eval ( "hessianfilter.SetSigma( " + str(parameters["Sigma"]["value"]) + ")" )
#
#	vesselnessfilter = itkHessian3DToVesselnessMeasureImageFilterF_New()
#	eval ( "vesselnessfilter.SetAlpha1( " + str(parameters["Alpha1"]["value"]) + ")" )
#	eval ( "vesselnessfilter.SetAlpha2( " + str(parameters["Alpha2"]["value"]) + ")" )
#
#	PipelineInputOutput([reader, hessianfilter, vesselnessfilter, writer])
#
#	writer.Update()
#
#	print thisalgorithm
#	print "Algorithm successfully executed!\n"

def SetElementFromList(object, list, size):
	for i in range(size):
		object.SetElement(i, list[i])
		
