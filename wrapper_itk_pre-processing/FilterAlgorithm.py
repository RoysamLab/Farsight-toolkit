# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This class contains the definition of the class filter object that will contain all the member functions and variables required for a filter , algorithm 

# Author 	: Adarsh K. Ramasubramonian
# Date		: 21 May, 2009

# Last modified on 21 May, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
from FilterObject import *		# Base class
import FilterAlgorithmGUIwx		# Contains all the widgets for the GUI.

class FilterAlgorithm(FilterObject):
	""" Filter algorithm - derived from Filter object, has additional data structures """
	# Public members
	def __init__(self):
		self.__outputType = True
		super(FilterAlgorithm, self).__init__()	# Super class constructor

		self.__requiresStructuringElement	= False
		self.__structuringElement 			= None
		self.__inputFileName 				= ""
		self.__outputFileName 				= ""
		self.__dimension 					= "2"			# Defualt - 2D image
		self.__inputImageStack 				= False
		self.__outputImageStack 			= False
		self.__outputImagePixelType 		= "UC"
		self.__useCaster 					= True			# Rescale filter output image
		self.__startIndex 					= 0
		self.__endIndex 					= 0
		self.__incrementIndex 				= 0


	def __str__(self):
		returnValue = super(FilterAlgorithm, self).__str__()
		returnValue += "\n" + "Input file name/format      : " + self.GetInputFileName()											
		returnValue += "\n" + "Output file name/format     : " + self.GetOutputFileName() 
		returnValue += "\n" + "Output Image type    : "
		if self.GetOutputImagePixelType() == "UC":
			returnValue += "8-bit"
		elif self.GetOutputImagePixelType() == "US":
			returnValue += "16-bit"
		elif self.GetOutputImagePixelType() == "F":
			returnValue += "float"
		returnValue += "\n" + "Input/output image dimension  : " + self.GetDimension()
		return returnValue

	# Long list of Get() and Set() functions for the member objects. Would be nice 
	# if there is a shortcut or a macro.

	def SetInputFileName(self, name):
		self.__inputFileName = name
	def SetOutputFileName(self, name):
		self.__outputFileName = name

	def GetInputFileName(self):
		return self.__inputFileName
	def GetOutputFileName(self):
		return self.__outputFileName 

	def GetStartIndex(self):
		return self.__startIndex
	def GetEndIndex(self):
		return self.__endIndex
	def GetIncrementIndex(self):
		return self.__incrementIndex

	def SetStartIndex(self, num):
		self.__startIndex = num
	def SetEndIndex(self, num):
		self.__endIndex = num
	def SetIncrementIndex(self, num):
		self.__incrementIndex = num

	def SetInputImageStack(self, state):
		self.__inputImageStack = state
	def SetOutputImageStack(self, state):
		self.__outputImageStack = state

	def GetInputImageStack(self):
		return self.__inputImageStack
	def GetOutputImageStack(self):
		return self.__outputImageStack

	def GetHasStructuringElement(self):
		return self.__requiresStructuringElement
	def SetHasStructuringElement(self, state):
		self.__requiresStructuringElement = state

	def SetStructuringElement(self, element):
		self.__structuringElement = element
	def GetStructuringElement(self):
		return self.__structuringElement

	def SetOutputType(self, toWhat):
		self.__outputType = toWhat
	def GetOutputType(self):
		return self.__outputType 

	def SetOutputImagePixelType(self, imageType):
		self.__outputImagePixelType = imageType
	def GetOutputImagePixelType(self):
		return self.__outputImagePixelType

	def SetDimension(self, dimension):
		self.__dimension = dimension
	def GetDimension(self):
		return self.__dimension

	def SetUseCaster(self, state):
		self.__useCaster = state
	def GetUseCaster(self):
		return self.__useCaster

	# GUI input - used with Tkinter - not any more.
#	def GUIInput(self):
#		root = Tk()
#		root.grid()
#		root.geometry("500x500")
#		app = FilterAlgorithmGUI.FilterAlgorithmGUI(root, self)
#		app.CreateWidgets()											
#		root.mainloop()

	def GUIInputwx(self):
		root = FilterAlgorithmGUIwx.FilterAlgorithmApp(False, self)
		root.MainLoop()

# interactiveInput needs to be activated. 					   
#	def interactiveInput(self):
#		print "Algorithm            : ", self.GetName()
#		print "Key                  : ", self.GetKey()
#		print "Number of parameters : ", self.numberParameters
#		keys = self.GetParameters().keys()
#		for num in range(self.numberParameters):
#			flag = 0
#			while flag == 0:
#				value = raw_input("Enter " + keys[num] + " : ")
#				if value:
#					flag = 1
#			self.AssignParameter(keys[num],value)
#		print "This function needs to be updated. Do not use it !!!"
#		import sys
#		sys.exit

