# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This file contains a list of URLs that will be accessed for each filter.

# Author: Adarsh K. Ramasubramonian
# Date  : 3 June, 2009
#
# Comments: The list of algorithms specified is only for initial purposes. Later, when all the webpages are created, we can just directly return the name.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import sys

listOfAlgorithms = ("Mean", "Median", "Laplacian", "Sigmoid", "GradientAnisotropicDiffusion", "CurvatureAnisotropicDiffusion", "CurvatureFlow", "MinMaxCurvatureFlow","RecursiveGaussian", "SmoothingRecursiveGaussian", "GrayscaleErode", "GrayscaleDilate", "Flip", "Normalize", "ShiftScale", "SobelEdgeDetection", "Resample")

class URLOfAlgorithmHelp(object):
	def __init__(self):
		self.MAIN_URL = "http://www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python"

	# This will be used to return the URL of the wiki page for an algorithm.
	def GetURLOfAlgorithm(self, name):
		if name.lower() == "mainpage":
			return self.MAIN_URL
		else:
			key = self.IsNameInList(name)
			if key:
				return self.MAIN_URL + "/" + key + "Filter"
			else:
				print "Please check name of the algorithm."
				sys.exit()
	
	# Check if given name of algorithm is the given list of algorithms.						  
	def IsNameInList(self, name):
		returnValue = ""
		for algorithm in listOfAlgorithms:
			if name.lower() == algorithm.lower():
				returnValue = algorithm
				break
		return returnValue
	


