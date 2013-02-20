# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# FilterAlgorithmGUIwx - 
#
# This class contains the definition of the class filter object that will contain 
# all the member functions and variables required for a FilterAlgorithm 

# Author 	: Adarsh K. Ramasubramonian
# Date		: 21 May, 2009

# Last modified on 3 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import os
import wx
import FilterObjectXML
import ExecuteAlgorithmModule
import FilterObjectGUIwx		# Super class
import HyperLinkModule 

# Application Object.
class FilterAlgorithmApp(wx.App):	# Constructor
	def __init__(self, state, algorithm):
		""" Constructor for FilterAlgorithmApp """
		wx.App.__init__(self, state)
		app = FilterAlgorithmGUIwx(None, -1, "Algorithm Details", algorithm)
		app.CreateWidgets()
		app.Show(True)

# Main GUI object					  
class FilterAlgorithmGUIwx(FilterObjectGUIwx.FilterObjectGUIwx):
	""" GUI for parameters for FilterAlgorithm using wxPython. """
	def __init__(self, parent, id, title, algorithm):
		""" Constructor for FilterAlgorithmGUIwx """
		super(FilterAlgorithmGUIwx, self).__init__(parent, id, title, algorithm)
		self.inputImageStack 	= False	
		self.outputImageStack 	= False	
		self.dirname 			= ''

	def CreateWidgets(self):
		self.DisplayAlgorithmLabel()            # From base class
		self.DisplayHelpURL()
		self.DisplayAlgorithmName()				# From base class
		self.DisplayAdvancedHelpURL()
		self.DisplayAlgorithmKey()				# From base class
		self.ParametersWidgets()				# From base class
		self.OutputImagePixelTypeWidgets()		# From base class
		self.ImageDimension()					# From base class
		self.InputOutputImageStack()
		self.InputFileNameWidgets()
		self.OutputFileNameWidgets()
		self.SeriesIndices()
		self.PressFilterButtonWidget()
		self.DisplayQuitButtonWidget()		

	def DisplayHelpURL(self):
		HyperLinkModule.hyperlink.HyperLinkCtrl(
            self.panel, -1, 
            label = "Help", 
            URL = self.algorithm.GetHelpURL(),
            pos = (400,self.tempRow + 15))	

	def DisplayAdvancedHelpURL(self):
		HyperLinkModule.hyperlink.HyperLinkCtrl(
            self.panel, -1, 
            label = "Advanced help", 
            URL = self.algorithm.GetAdvancedHelpURL(),
            pos = (400,self.tempRow + 15))	

	def DisplayQuitButtonWidget(self):
		self.NextRow()
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Press Quit to exit")
		self.quitButton = wx.Button(self.panel, id = 1, label = "Quit", pos = (200, self.tempRow), size = (75, 40))
		self.quitButton.Bind(wx.EVT_BUTTON, self.QuitNow)
		self.quitButton.SetToolTip(wx.ToolTip("Click 'Quit' to quit."))

	def QuitNow(self, event):
		self.Destroy()

	def InputOutputImageStack(self):
		""" Widget for enabling input/output images to be stacks rather than a single image """
		self.NextRow()
		self.inputImageStackCheckBox = wx.CheckBox(self.panel, -1, "Input Stack Image", pos = (20, self.tempRow))
		self.Bind(wx.EVT_CHECKBOX, self.updateInputImageStack, self.inputImageStackCheckBox)
										
		self.outputImageStackCheckBox = wx.CheckBox(self.panel, -1, "Output Stack Image",  pos = (150, self.tempRow))
		self.Bind(wx.EVT_CHECKBOX, self.updateOutputImageStack, self.outputImageStackCheckBox)
		self.inputImageStackCheckBox.SetValue(wx.CHK_UNCHECKED)
		self.inputImageStackCheckBox.Enable(False)
		self.outputImageStackCheckBox.SetValue(wx.CHK_UNCHECKED)
		self.outputImageStackCheckBox.Enable(False)

	def updateInputImageStack(self, event):
		if event.GetEventObject().GetValue():
			self.inputImageStack = True
			self.inputFileNameLabel.SetLabel("Enter input file format")
			self.EnableFileIndexWidgets()
		else:
			self.inputImageStack = False
			self.inputFileNameLabel.SetLabel("Enter input file name")
			if not self.outputImageStack:		# If neither input nor output is stacked, disable stack indices
				self.DisableFileIndexWidgets()

	def updateOutputImageStack(self, event):
		if event.GetEventObject().GetValue():
			self.outputImageStack = True
			self.outputFileNameLabel.SetLabel("Enter output file format")
			self.EnableFileIndexWidgets()
		else:
			self.outputImageStack = False
			self.outputFileNameLabel.SetLabel("Enter output file name")
			if not self.inputImageStack:
				self.DisableFileIndexWidgets()	# If neither input nor output is stacked, disable stack indices

	# Image dimension - may have to change to input and output dimension later.																	   
	def ImageDimension(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Input/Output Image dimension")
		radio2D = wx.RadioButton(self.panel, -1, '2D', pos = (200, self.tempRow), style = wx.RB_GROUP)
		radio3D = wx.RadioButton(self.panel, -1, '3D', pos = (250, self.tempRow))
		for eachRadio in [radio2D, radio3D]:
			self.Bind(wx.EVT_RADIOBUTTON, self.SetDimension, eachRadio)
	
	def SetDimension(self, event):
		localDimensionObject = event.GetEventObject()
		if localDimensionObject.GetLabel() == "2D":		# For 2D images, no stacks - hence disable them
			self.algorithm.SetDimension("2")
			self.inputImageStackCheckBox.SetValue(wx.CHK_UNCHECKED)
			self.inputImageStackCheckBox.Enable(False)
			self.outputImageStackCheckBox.SetValue(wx.CHK_UNCHECKED)
			self.outputImageStackCheckBox.Enable(False)
		elif localDimensionObject.GetLabel() == "3D":
			self.algorithm.SetDimension("3")
			self.inputImageStackCheckBox.Enable(True)
			self.outputImageStackCheckBox.Enable(True)
	
	def InputFileNameWidgets(self):
		self.NextRow()
		self.inputFileNameLabel = wx.StaticText(self.panel, -1, "Enter input file name", pos = (20, self.tempRow))

		self.inputFileNameEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
		if self.algorithm.GetInputFileName():
			self.inputFileNameEntry.SetValue(self.algorithm.GetInputFileName())

		self.inputFileButton = wx.Button(self.panel, id = 1, label = "Browse", pos = (400, self.tempRow), size = (75, 20))
		self.inputFileButton.Bind(wx.EVT_BUTTON, self.AskInputFileName)
		self.inputFileButton.SetToolTip(wx.ToolTip("Browse for an input image"))

	def OutputFileNameWidgets(self):
		self.NextRow()
		self.outputFileNameLabel = wx.StaticText(self.panel, -1, "Enter output file name", pos = (20, self.tempRow))

		self.outputFileNameEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
		if self.algorithm.GetOutputFileName():
			self.outputFileNameEntry.SetValue(self.algorithm.GetOutputFileName())

		self.outputFileButton = wx.Button(self.panel, id = 1, label = "Browse", pos = (400, self.tempRow), size = (75, 20))
		self.outputFileButton.Bind(wx.EVT_BUTTON, self.AskOutputFileName)
		self.outputFileButton.SetToolTip(wx.ToolTip("Browse for an output image"))

	def AskInputFileName(self, event):
		self.inputFileNameEntry.SetValue(self.AskFileName("Choose a file"))

	def AskOutputFileName(self, event):
		self.outputFileNameEntry.SetValue(self.AskFileName("Choose a file"))

	def AskFileName(self,title):			# File Dialog widget
		dlg = wx.FileDialog(self, title, self.dirname, "", "*.*", wx.OPEN)
		if dlg.ShowModal() == wx.ID_OK:
			self.filename = dlg.GetFilename()
			self.dirname = dlg.GetDirectory()
			return (os.path.join(self.dirname, self.filename))
	
	def SeriesIndices(self):				# Indices for the image stack
		self.NextRow()
		self.startIndexLabel = wx.StaticText(self.panel, -1, "Start Index", pos = (20, self.tempRow))
		self.startIndexEntry = wx.TextCtrl(self.panel, -1, pos = (200, self.tempRow))

		self.NextRow()
		self.endIndexLabel = wx.StaticText(self.panel, -1, "End Index", pos = (20, self.tempRow))
		self.endIndexEntry = wx.TextCtrl(self.panel, -1, pos = (200, self.tempRow))

		self.NextRow()
		self.incrementIndexLabel = wx.StaticText(self.panel, -1, "Increment Index", pos = (20, self.tempRow))
		self.incrementIndexEntry = wx.TextCtrl(self.panel, -1, pos = (200, self.tempRow))

		self.DisableFileIndexWidgets()
	
	def DisableFileIndexWidgets(self):          # File indices only for input/output stack images.
		self.startIndexEntry.Enable(False)
		self.endIndexEntry.Enable(False)
		self.incrementIndexEntry.Enable(False)

		self.startIndexLabel.Enable(False)
		self.endIndexLabel.Enable(False)
		self.incrementIndexLabel.Enable(False)

	def EnableFileIndexWidgets(self):
		self.startIndexEntry.Enable(True)
		self.endIndexEntry.Enable(True)
		self.incrementIndexEntry.Enable(True)

		self.startIndexLabel.Enable(True)
		self.endIndexLabel.Enable(True)
		self.incrementIndexLabel.Enable(True)

	def StartFiltering(self, event):			# When Filter is pressed, get the algorithm details.
                                                  
		self.algorithm.SetInputFileName( self.inputFileNameEntry.GetValue() )
		self.algorithm.SetOutputFileName( self.outputFileNameEntry.GetValue() )

		self.algorithm.SetInputImageStack(self.inputImageStack)
		self.algorithm.SetOutputImageStack(self.outputImageStack)

		if self.inputImageStack or self.outputImageStack:
			self.algorithm.SetStartIndex( eval( self.startIndexEntry.GetValue() ) )
			self.algorithm.SetEndIndex( eval( self.endIndexEntry.GetValue() ) )
			self.algorithm.SetIncrementIndex( eval( self.incrementIndexEntry.GetValue() ) )
		for num in range(self.algorithm.numberParameters):
			if not self.parameterEntries[num].GetValue():
				print "Please enter all the parameters"
				return
			else:
				 self.algorithm.AssignParameter(
						 self.keys[num],
						 self.parameterEntries[num].GetValue() )
		self.algorithm.CheckParameterValidity()		# Not implemented currently			
																				
		self.filterButton.Enable(False)
		xmlObject = FilterObjectXML.FilterObjectXML()
		XMLFileName = "parameters.xml"
				# XXX To be generalized
		xmlObject.UpdateDocument(self.algorithm, XMLFileName)
		ExecuteAlgorithmModule.defaultAlgorithmExecution(XMLFileName)

		self.filterButton.Enable(True)

