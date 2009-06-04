# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# process.py - this python script is written in order to execute different algorithms
#             in the ITK framework. The idea is to  let the user type very simple commands. 
# 			  For further help on the execution of the algorithm, please check 
# 			  www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python.
# 			  The wiki page also includes details of what each algorithm does and how to run
#			  them.
# 
#
# Author 	: Adarsh K. Ramasubramonian
# Date 		: 15 May, 2009
#
# Things to do:
# 1. Export algorithm details to 
#
# Last modified on 3 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import sys
import FilterAlgorithm # Algorithm class
import ListOfItems
import FilterObjectXML
from sys import argv			# Command line arguments
from basic import *
from AlgorithmDetail import *
# import pdb						# For debugging


# List of algorithms included
listOfAlgorithms = ("Mean", "Median", "Laplacian", "Sigmoid", "GradientAnisotropicDiffusion", "CurvatureAnisotropicDiffusion", "CurvatureFlow", "MinMaxCurvatureFlow","RecursiveGaussian", "SmoothingRecursiveGaussian", "GrayscaleErode", "GrayscaleDilate", "Flip", "Normalize", "ShiftScale", "SobelEdgeDetection", "Resample")

XMLFileName = "parameters.xml"
xmlObject = FilterObjectXML.FilterObjectXML()

def defaultAlgorithmExecution(thisalgorithm):

	if thisalgorithm.GetKey() == "resample":					# Special steps for Resampling
		AlgorithmDetailResampleFilter(thisalgorithm)
	else:
		ExecuteCommands(thisalgorithm)							# Default algorithm execution


def usage():
	# Incorrect number of arguments
	print \
	"""Incorrect number of arguments"
	Correct usage: *.py algorithm_key input_filename output_filename (The input and output filenames are optional)
	"""
	sys.exit()



def GetAlgorithmDetails(choiceOfAlgorithm):
	choiceOfAlgorithm = choiceOfAlgorithm.lower()
	dprint("Choice of algorithm :" + choiceOfAlgorithm)
	
	# List of all the main algorithms and their essential code. - this needs to be moved into an XML file, along with an interface that will permit anybody to add a new algorithm with its parameters to the XML file.
	
	algorithm = FilterAlgorithm.FilterAlgorithm()		# Algorithm Object
	if len(argv) > 2:
		print "Input : " + argv[2]
		algorithm.SetInputFileName(argv[2])
	if len(argv) > 3:
		print "Output: " + argv[3]
		algorithm.SetOutputFileName(argv[3])
	
	if choiceOfAlgorithm == "mean":								# Mean filter
		algorithm.SetKey("mean") 
		algorithm.SetName("itkMeanImageFilter")
		# algorithm.AddParameters(["xRadius", "yRadius"])
		algorithm.AddTupleParameter("Radius", 2, "itkSize")
		defaultExecution = True
	
	elif choiceOfAlgorithm == "median":							# Median filter
		algorithm.SetKey("median")
		algorithm.SetName("itkMedianImageFilter")
		algorithm.AddTupleParameter("Radius", 2, "itkSize") 
		defaultExecution = True
	
	elif choiceOfAlgorithm == "laplacian":						# Laplacian filter
		algorithm.SetKey("laplacian")
		algorithm.SetName("itkLaplacianImageFilter")
		defaultExecution = True
	
	elif choiceOfAlgorithm == "sigmoid":						# Sigmoid filter
		algorithm.SetKey("sigmoid")
		algorithm.SetName("itkSigmoidImageFilter")
		algorithm.AddParameters(["OutputMinimum", "OutputMaximum", "Alpha", "Beta"])
		defaultExecution = True
	
	elif choiceOfAlgorithm == "gradientanisotropicdiffusion":	# Gradient Anisotropic Diffusion filter
		algorithm.SetKey("gradientanisotropicdiffusion")
		algorithm.SetName("itkGradientAnisotropicDiffusionImageFilter")
		algorithm.AddParameters(["NumberOfIterations", "TimeStep", "ConductanceParameter"])
		defaultExecution = True
	
	elif choiceOfAlgorithm == "curvatureanisotropicdiffusion":	# Gradient Anisotropic Diffusion filter
		algorithm.SetKey("curvatureanisotropicdiffusion")
		algorithm.SetName("itkCurvatureAnisotropicDiffusionImageFilter")
		algorithm.AddParameters(["NumberOfIterations", "TimeStep", "ConductanceParameter"])
		defaultExecution = True

	elif choiceOfAlgorithm == "curvatureflow":					# Curvature Flow filter
		algorithm.SetKey("curvatureflow")
		algorithm.SetName("itkCurvatureFlowImageFilter")
		algorithm.AddParameters(["NumberOfIterations", "TimeStep"])
		defaultExecution = True
	
	elif choiceOfAlgorithm == "minmaxcurvatureflow":			# Min-max curvature flow
		algorithm.SetKey("minmaxcurvatureflow")
		algorithm.SetName("itkMinMaxCurvatureFlowImageFilter")
		algorithm.AddParameters(["NumberOfIterations", "TimeStep", "StencilRadius"]) # XXX
		defaultExecution = True
	
	
	elif choiceOfAlgorithm == "smoothingrecursivegaussian":		# Smoothing Recursive Gaussian filter
		algorithm.SetKey("smoothingrecursivegaussian")
		algorithm.SetName("itkSmoothingRecursiveGaussianImageFilter")
		algorithm.AddParameters(["Sigma","NumberOfThreads", "NormalizeAcrossScale"])
		defaultExecution = True
	
	elif choiceOfAlgorithm == "grayscaleerode":					# Erode
		algorithm.SetKey("grayscaleerode")
		algorithm.SetName("itkGrayscaleErodeImageFilter")
		algorithm.SetHasStructuringElement(True)
		algorithm.AddParameter("Radius")
		defaultExecution = True
	
	elif choiceOfAlgorithm == "grayscaledilate":				# Dilate
		algorithm.SetKey("grayscaledilate")
		algorithm.SetName("itkGrayscaleDilateImageFilter")
		algorithm.SetHasStructuringElement(True)
		algorithm.AddParameter("Radius")
		defaultExecution = True

	elif choiceOfAlgorithm == "flip":							# Flip filter
		algorithm.SetKey("flip")
		algorithm.SetName("itkFlipImageFilter")
		algorithm.AddTupleParameter("FlipAxes",2,"itkFixedArrayB")
		algorithm.SetOutputType(None)
		defaultExecution = True

	elif choiceOfAlgorithm == "normalize":						# Normalize filter
		algorithm.SetKey("normalize")
		algorithm.SetName("itkNormalizeImageFilter")
		defaultExecution = True

	elif choiceOfAlgorithm == "shiftscale":						# Shiftscale filter
		algorithm.SetKey("shiftscale")
		algorithm.SetName("itkShiftScaleImageFilter")
		algorithm.AddParameters(["Shift","Scale"])
		algorithm.SetUseCaster(False)
		defaultExecution = True

	elif choiceOfAlgorithm == "sobeledgedetection":				# Sobel Edge Detector
		algorithm.SetKey("sobeledgedetection")
		algorithm.SetName("itkSobelEdgeDetectionImageFilter")
		algorithm.AddParameter("NumberOfThreads")
		defaultExecution = True
	
	elif choiceOfAlgorithm == "resample":						# Resample to remove anisotropy
		algorithm.SetKey("resample")
		algorithm.SetName("itkResampleImageFilter")
		defaultExecution = True

	
	dprint("choiceOfAlgorithm was "+choiceOfAlgorithm)
	
	if defaultExecution:
		if interactiveMode:
			algorithm.interactiveInput()
		else:
			algorithm.GUIInputwx()
		xmlObject.UpdateDocument(algorithm, XMLFileName)
	

# This function checks if the argument name contains the extension .xml, in which case it returns a 1.
def CheckXMLFile( name):
	lengthOfName = len(name)
	if name[lengthOfName-4:] == ".xml":							 
		return True
	else:
		return False

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 					
# main
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 					

interactiveMode = False 					# Command line UI
GUIMode 		= not interactiveMode		# GUI

if len(argv) < 1 or len(argv) > 4:
	usage()

if len(argv) == 1:
	choiceOfAlgorithm  = ListOfItems.DisplayListOfAlgorithms(listOfAlgorithms)
	GetAlgorithmDetails(choiceOfAlgorithm)
else:
	# Check if the first argument passed is an XML file.
	if CheckXMLFile(argv[1]):
		XMLFileName = argv[1]
	else:
		choiceOfAlgorithm = argv[1]					# Choice of algorithm/filter
		GetAlgorithmDetails(choiceOfAlgorithm)

newalgorithm = FilterAlgorithm.FilterAlgorithm()			# New algorithm object
xmlObject.ParseDocument( newalgorithm, XMLFileName)			# Parse XML file to read the parameters
defaultAlgorithmExecution(newalgorithm)						# Execute the algorithm

raw_input("Press a key to exit")




#	elif choiceOfAlgorithm == "recursivegaussian":				# Recursive Gaussian filter
#		algorithm.SetKey("recursivegaussian")
#		algorithm.SetName("itkRecursiveGaussianImageFilter")
#		algorithm.AddParameters(["Sigma","Order", "NormalizeAcrossScale"])
#		defaultExecution = True
