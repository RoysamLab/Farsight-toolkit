# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ExecuteAlgorithmWrapITK.py - 

# The execution of the algorithm is performed using WrapITK instead of CableSwig.
# The default WrapITK installation does not include unsigned char; however, saving 
# files as 8-bits is possible. The filters are not defined for template parameters
# with unsigned char. Moreover, some of the filters don't need rescale intensities.
# In such cases, we could either have (i) a special case for unsigned char output, or 
# (ii) separate the implementations of Rescale intensities and Image Caster
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import ExternalProjects
import itk
import sys, cmath
import FilterObjectXML
import FilterAlgorithm
import pdb
import warnings

READER = "ImageFileReader"
READER_PIXELTYPE = "F"

WRITER = "ImageFileWriter"


def defaultAlgorithmExecution(XMLFileName):
    """ Details of algorithm in XML file. Read it and execute the algorithm """
    xmlObject = FilterObjectXML.FilterObjectXML()
    thisalgorithm = FilterAlgorithm.FilterAlgorithm()
    xmlObject.ParseDocument(thisalgorithm, XMLFileName)

    if thisalgorithm.GetKey() == "resample":	# Special steps for Resampling
        ExecuteResampleAlgorithm(thisalgorithm)
#    elif thisalgorithm.GetKey() == "anisotropicdiffusionvesselenhancement":
#        print "Vesselness algorithm not tested. Sorry !!!"
#        ExecuteVesselnessAlgorithm(thisalgorithm)
    else:
        ExecuteCommands(thisalgorithm)			# Default algorithm execution

def ExecuteCommands(thisalgorithm):
    """ Default execution of commands """
    
    dimension = thisalgorithm.GetDimension()
    pixelType = thisalgorithm.GetOutputImagePixelType()
    
    reader = GetReader(thisalgorithm, dimension, pixelType)         # Reader Object
    writer = GetWriter(thisalgorithm, dimension, pixelType)         # Writer Object
    
    
    parameters = thisalgorithm.GetParameters()
    keys = parameters.keys()[:]
    
    if thisalgorithm.GetHasStructuringElement():
        keys.remove("Radius")
        radius = eval( str( parameters["Radius"]["value"] ) )
        kernel = GetKernel(dimension, radius)
        filter = GetFilterWithKernel(thisalgorithm, dimension, pixelType, kernel)
    else:
        filter = GetFilter(thisalgorithm, dimension, pixelType)         # Filter Object
    
    
    for key in keys:
        if key == "Order":
            orderVar = eval( parameters[key]["value"] )
            AssignOrderOfFilter(filter, orderVar)
    
        else:
        # Added on 19 May, 2009
            if thisalgorithm.IsTupleParameter(key):
                value = parameters[key]["value"]
                value = eval(str(value))            # Converted to tuple
                eval( "filter.Set" + key + "( value )" )
                    # In wrapITK, we can directly pass the tuple as argument.
            else:
                eval( "filter.Set" + key + "( " + str(parameters[key]["value"]) + " )" )
    
    
    if not thisalgorithm.GetUseRescaler():          # Function name should be changed to rescaler
        if pixelType != "F":
            warnings.warn("\nWARNING: The output may lose some details if not stored as float type !!!\n")
                                  
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
		READER = "itk.ImageSeriesReader"
		READER_PIXELTYPE = "F"

															#	ImageSeriesReader			
		inputPixelType = "itk." + READER_PIXELTYPE
		inputImageType = "itk.Image[" + inputPixelType + "," + dimension +"]"
		reader = eval( READER + "[" + inputImageType + "].New()" )
		readerNameGenerator = itk.NumericSeriesFileNames.New()			# Name generator
		readerNameGenerator.SetSeriesFormat( thisalgorithm.GetInputFileName() )														
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(readerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1) doesn't seem to be working
		reader.SetFileNames( readerNameGenerator.GetFileNames() )		# Set Filename
								# XXX Need to do something.

	elif dimension == "2" or dimension == "3":
	# If single image
		READER = "itk.ImageFileReader"
		READER_PIXELTYPE = "F"
		inputPixelType = "itk." + READER_PIXELTYPE
		inputImageType = "itk.Image[" + inputPixelType + "," + dimension + "]"
		reader = eval( READER + "[" + inputImageType + "].New()" )
		reader.SetFileName( thisalgorithm.GetInputFileName() )			# Set Filename

	return reader


def GetWriter(thisalgorithm, dimension, pixelType):
	""" Get Writer Object to write the output file """
	# If input stack
	if dimension == "3" and thisalgorithm.GetOutputImageStack():

		WRITER = "itk.ImageSeriesWriter"
		WRITER_PIXELTYPE = pixelType
		inputPixelType = "itk." + pixelType
		inputImageType = "itk.Image[" + inputPixelType + "," + dimension + "]"
		outputPixelType = "itk." + pixelType 
		outputImageType = "itk.Image[" + outputPixelType + "," + "2" + "]"
		writer = eval( WRITER + "[" + inputImageType + "," + outputImageType + "].New()" )
		writerNameGenerator = itk.NumericSeriesFileNames.New()							# Name generator
		writerNameGenerator.SetSeriesFormat( thisalgorithm.GetOutputFileName() )
	# 	Set start, end, and increment index, and the format
		SetSeriesAndIndices(writerNameGenerator, thisalgorithm)

	# 	Assign IO model		- XXX (1)  doesn't seem to be working

		writer.SetFileNames( writerNameGenerator.GetFileNames() ) 		# Set file name

	elif dimension == "2" or dimension == "3":							# If single image
		WRITER = "itk.ImageFileWriter"
		WRITER_PIXELTYPE = pixelType
		outputPixelType = "itk." + pixelType 
		outputImageType = "itk.Image[" + outputPixelType + "," + dimension + "]"
		writer = eval( WRITER + "[" + outputImageType + "].New()" )
		writer.SetFileName( thisalgorithm.GetOutputFileName() )			# Set FileName

	return writer


def GetFilter(thisalgorithm, dimension, pixelType):	
	""" Get the filter object - the main algorithm """
	filterName = thisalgorithm.GetName()
	filterName = filterName[0:3] + "." + filterName[3:]
	inputPixelType = "itk." + READER_PIXELTYPE
	inputImageType = "itk.Image[" + inputPixelType + "," + dimension + "]"
	if thisalgorithm.GetOutputType():    
		outputPixelType= "itk." + READER_PIXELTYPE 
		outputImageType= "itk.Image[" + outputPixelType + "," + dimension + "]"
		filter = eval( filterName + "[" + inputImageType + "," + outputImageType + "].New()" )
	else:
		filter = eval( filterName + "[" + inputImageType + "].New()" )
	return filter

def GetFilterWithKernel(thisalgorithm, dimension, pixelType, kernel):	
	""" Get the filter object - the main algorithm """
	filterName = thisalgorithm.GetName()
	filterName = filterName[0:3] + "." + filterName[3:]
	inputPixelType = "itk." + READER_PIXELTYPE
	inputImageType = "itk.Image[" + inputPixelType + "," + dimension + "]"
	if thisalgorithm.GetOutputType():    
		outputPixelType= "itk." + READER_PIXELTYPE 
		outputImageType= "itk.Image[" + outputPixelType + "," + dimension + "]"
		filter = eval( filterName + "[" + inputImageType + "," + outputImageType + ",kernel].New()" )
	else:
		filter = eval( filterName + "[" + inputImageType + ",kernel].New()" )
	filter.SetKernel(kernel)
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


def GetKernel(dimension, radius):
	element = itk.strel( int(dimension), radius)
	return element


def SetSeriesAndIndices(Generator, algorithm):
	Generator.SetStartIndex( algorithm.GetStartIndex() )
	Generator.SetEndIndex( algorithm.GetEndIndex() )
	Generator.SetIncrementIndex( algorithm.GetIncrementIndex() )


def GetRescaler( inputPixelType, dimension, outputPixelType):

	inputPixelType = eval( "itk." + inputPixelType )
	outputPixelType= eval( "itk." + outputPixelType )
	inputImageType = itk.Image[inputPixelType, int(dimension)]
	outputImageType= itk.Image[outputPixelType, int(dimension)]
	rescaler = itk.RescaleIntensityImageFilter[inputImageType,outputImageType].New()
	
	rescaler.SetOutputMinimum(   0 )
	if outputPixelType == "UC":
		rescaler.SetOutputMaximum( 255 )
	elif outputPixelType == "US":
		rescaler.SetOutputMaximum( 65535 )
	return rescaler

def PipelineInputOutput(listOfObjects):
	for i in range(len(listOfObjects)-1):
		listOfObjects[i+1].SetInput( listOfObjects[i].GetOutput() )

def SetElementFromList(object, list, size):
	for i in range(size):
		object.SetElement(i,list[i])


def ExecuteResampleAlgorithm(thisalgorithm):

	dimension = thisalgorithm.GetDimension()
	pixelType = thisalgorithm.GetOutputImagePixelType()

	reader = GetReader(thisalgorithm, dimension, pixelType)
	writer = GetWriter(thisalgorithm, dimension, pixelType)

	itkpixelType = "itk." + READER_PIXELTYPE
	imageType = "itk.Image[" + itkpixelType + "," + dimension + "]"

	# Define the smoothing filters.
	smootherX = eval( "itk.RecursiveGaussianImageFilter[" + imageType + "," + imageType + "].New()" )
	smootherY = eval( "itk.RecursiveGaussianImageFilter[" + imageType + "," + imageType + "].New()" )


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
	transform = eval( "itk.IdentityTransform[itk.D, " + dimension + "].New()" )
	transform.SetIdentity()
	filter.SetTransform( transform.GetPointer() )

	# Define and set interpolator

	interpolator = eval( "itk.LinearInterpolateImageFunction[" + imageType + ",itk.D].New()" )
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

	size = eval( "itk.Size[" + dimension + "]()" )
	SetElementFromList(size, [dx, dy, dz], 3)
	filter.SetSize(size)

	if not pixetType == "F":
		print """\nWARNING: The output may lose some details if not stored as float type !!!\n """
		
	rescaler = GetRescaler(READER_PIXELTYPE, dimension, pixelType)

	PipelineInputOutput([reader, smootherX, smootherY, filter, rescaler, writer])

	writer.Update()
	print thisalgorithm
	print "Algorithm successfully executed!\n"


