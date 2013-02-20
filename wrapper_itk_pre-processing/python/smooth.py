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
# 1. Export algorithm details to XML file - DONE
# 2. Allow user to specify XML for output - TBD XXX
#
# Last modified on 3 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import sys
import os
import FilterAlgorithm # Algorithm class
import ListOfAlgorithmsWithURL 
import FilterAlgorithmGUIwx
import basic 
import ExecuteAlgorithmModule
import pdb						# For debugging
import AlgorithmXML

# List of algorithms included
# listOfAlgorithms = ("Mean", "Median", "Laplacian", "Sigmoid", "GradientAnisotropicDiffusion", "CurvatureAnisotropicDiffusion", "CurvatureFlow", "MinMaxCurvatureFlow", "SmoothingRecursiveGaussian", "GrayscaleErode", "GrayscaleDilate", "Flip", "Normalize", "ShiftScale", "SobelEdgeDetection", "Resample", "Vesselness")

XMLFileName = "parameters.xml"
USINGXMLFILEFORDETAILS = True
algorithmXMLObject = AlgorithmXML.AlgorithmXML()

def usage():
	# Incorrect number of arguments
	print \
	"""Incorrect number of arguments"
	Correct usage: *.py algorithm_key input_filename output_filename (The input and output filenames are optional)
	"""
	sys.exit()


def GetAlgorithmDetails(choiceOfAlgorithm, argv):
	global XMLFileName 
	global algorithmXMLObject
	choiceOfAlgorithm = choiceOfAlgorithm.lower()
	basic.dprint("Choice of algorithm :" + choiceOfAlgorithm)
                                                       

	if not USINGXMLFILEFORDETAILS:
		# List of all the main algorithms and their essential code. - this needs to be moved into an XML file, along with an interface that will permit anybody to add a new algorithm with its parameters to the XML file. - DONE.
		
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
			algorithm.AddTupleParameter("Radius", 2, "itkSize")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/MeanFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.6.1")
		
		elif choiceOfAlgorithm == "median":							# Median filter
			algorithm.SetKey("median")
			algorithm.SetName("itkMedianImageFilter")
			algorithm.AddTupleParameter("Radius", 2, "itkSize") 
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/MedianFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.6.2")
		
		elif choiceOfAlgorithm == "laplacian":						# Laplacian filter
			algorithm.SetKey("laplacian")
			algorithm.SetName("itkLaplacianImageFilter")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/LaplacianFilter")
		
		elif choiceOfAlgorithm == "sigmoid":						# Sigmoid filter
			algorithm.SetKey("sigmoid")
			algorithm.SetName("itkSigmoidImageFilter")
			algorithm.AddParameters(["OutputMinimum", "OutputMaximum", "Alpha", "Beta"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/SigmoidFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.3.2")
		
		elif choiceOfAlgorithm == "gradientanisotropicdiffusion":	# Gradient Anisotropic Diffusion filter
			algorithm.SetKey("gradientanisotropicdiffusion")
			algorithm.SetName("itkGradientAnisotropicDiffusionImageFilter")
			algorithm.AddParameters(["NumberOfIterations", "TimeStep", "ConductanceParameter"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/GradientAnisotropicDiffusionFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.7.3")
		
		elif choiceOfAlgorithm == "curvatureanisotropicdiffusion":	# Curvature Anisotropic Diffusion filter
			algorithm.SetKey("curvatureanisotropicdiffusion")
			algorithm.SetName("itkCurvatureAnisotropicDiffusionImageFilter")
			algorithm.AddParameters(["NumberOfIterations", "TimeStep", "ConductanceParameter"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/CurvatureAnisotropicDiffusionFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.7.3")
	
		elif choiceOfAlgorithm == "curvatureflow":					# Curvature Flow filter
			algorithm.SetKey("curvatureflow")
			algorithm.SetName("itkCurvatureFlowImageFilter")
			algorithm.AddParameters(["NumberOfIterations", "TimeStep"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/CurvatureFlowFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.7.3")
		
		elif choiceOfAlgorithm == "minmaxcurvatureflow":			# Min-max curvature flow
			algorithm.SetKey("minmaxcurvatureflow")
			algorithm.SetName("itkMinMaxCurvatureFlowImageFilter")
			algorithm.AddParameters(["NumberOfIterations", "TimeStep", "StencilRadius"]) # XXX
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/MinMaxCurvatureFlowFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.7.3")
		
		
		elif choiceOfAlgorithm == "smoothingrecursivegaussian":		# Smoothing Recursive Gaussian filter
			algorithm.SetKey("smoothingrecursivegaussian")
			algorithm.SetName("itkSmoothingRecursiveGaussianImageFilter")
			algorithm.AddParameters(["Sigma","NumberOfThreads", "NormalizeAcrossScale"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/SmoothingRecursiveGaussianFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.7.1")
		
		elif choiceOfAlgorithm == "grayscaleerode":					# Erode
			algorithm.SetKey("grayscaleerode")
			algorithm.SetName("itkGrayscaleErodeImageFilter")
			algorithm.SetHasStructuringElement(True)
			algorithm.AddParameter("Radius")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/GrayscaleErodeFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.6.3")
		
		elif choiceOfAlgorithm == "grayscaledilate":				# Dilate
			algorithm.SetKey("grayscaledilate")
			algorithm.SetName("itkGrayscaleDilateImageFilter")
			algorithm.SetHasStructuringElement(True)
			algorithm.AddParameter("Radius")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/GrayscaleDilateFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.6.3")
	
		elif choiceOfAlgorithm == "flip":							# Flip filter
			algorithm.SetKey("flip")
			algorithm.SetName("itkFlipImageFilter")
			algorithm.AddTupleParameter("FlipAxes",2,"itkFixedArrayB")
			algorithm.SetOutputType(None)
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/FlipFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.9.3")
	
		elif choiceOfAlgorithm == "normalize":						# Normalize filter
			algorithm.SetKey("normalize")
			algorithm.SetName("itkNormalizeImageFilter")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/NormalizeFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.3.1")
	
		elif choiceOfAlgorithm == "shiftscale":						# Shiftscale filter
			algorithm.SetKey("shiftscale")
			algorithm.SetName("itkShiftScaleImageFilter")
			algorithm.AddParameters(["Shift","Scale"])
			algorithm.SetUseRescaler(False)
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/ShiftScaleFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.3.1")
	
		elif choiceOfAlgorithm == "sobeledgedetection":				# Sobel Edge Detector
			algorithm.SetKey("sobeledgedetection")
			algorithm.SetName("itkSobelEdgeDetectionImageFilter")
			algorithm.AddParameter("NumberOfThreads")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/SobelEdgeDetectionFilter")
		
		elif choiceOfAlgorithm == "resample":						# Resample to remove anisotropy
			algorithm.SetKey("resample")
			algorithm.SetName("itkResampleImageFilter")
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/ResampleFilter")
			algorithm.SetAdvancedHelpURL("www.itk.org/ItkSoftwareGuide.pdf#6.9.4")
	
		elif choiceOfAlgorithm == "vesselness":
			algorithm.SetKey("vesselness")
			algorithm.SetName("itkHessian3DToVesselnessMeasureImageFilter")
			algorithm.AddParameters(["Sigma","Alpha1","Alpha2"])
			defaultExecution = True
			algorithm.SetHelpURL("www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python/VesselnessFilter")
	else:
		algorithm = algorithmXMLObject.GetAlgorithm(choiceOfAlgorithm)     
		if len(argv) > 2:
			print "Input : " + argv[2]
			algorithm.SetInputFileName(argv[2])
		if len(argv) > 3:
			print "Output: " + argv[3]
			algorithm.SetOutputFileName(argv[3])
	   	defaultExecution = True

	if defaultExecution:
		root = FilterAlgorithmGUIwx.FilterAlgorithmApp(False, algorithm)
		root.MainLoop()
	else:
		# Nothing as of now. May be used later
		print "Should not have come here in smooth.py. Check !"
		sys.exit()
	


def main(argv = None):
	global XMLFileName 
	global algorithmXMLObject
	if argv is None:
		argv = sys.argv

	if len(argv) < 1 or len(argv) > 4:
		usage()
	

	if len(argv) == 1:
		listOfAlgorithms = algorithmXMLObject.GetListOfAlgorithms()
		choiceOfAlgorithm  = ListOfAlgorithmsWithURL.DisplayList(listOfAlgorithms)
		GetAlgorithmDetails(choiceOfAlgorithm, argv)
	else:
		if basic.CheckFileExtension(argv[1], ".xml"):             # If XML file is passed as input.
			XMLFileName = argv[1]
			ExecuteAlgorithmModule.defaultAlgorithmExecution(XMLFileName)
		else:
			choiceOfAlgorithm = argv[1]			
			GetAlgorithmDetails(choiceOfAlgorithm, argv)
	
	raw_input("Press a key to exit")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 					
# main
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 					
if __name__ == "__main__":
	main()
	
