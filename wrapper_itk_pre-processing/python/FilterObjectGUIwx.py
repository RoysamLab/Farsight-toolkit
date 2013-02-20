# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# FilterObjectGUIwx - 
# This python module contains definition of the class that defines the GUI required for the 
# parameters of different algorithms, using the wxPython module.

# Author 	: Adarsh K. Ramasubramonian
# Date 		: 31 May, 2009
# 
# Last modified on 3 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import wx
import HyperLinkModule

class FilterObjectApp(wx.App):
	def __init__(self, state, algorithm):
		wx.App.__init__(self,state)
		app = FilterObjectGUIwx(None, -1, 'Algorithm Details', algorithm)
		app.CreateWidgets()
		app.Show(True)


class FilterObjectGUIwx(wx.Frame):
	""" GUI for parameters of FilterObject """
	def __init__(self, parent, id, title, algorithm):
		wx.Frame.__init__(self, parent, id, title, size = (600,600))
										# Super-class constructor
		self.parameterEntries 	= []	# Default values
		self.tempRow 			= 0
		self.algorithm 			= algorithm
		self.tempRow 			= 0
		self.panel 				= wx.Panel(self, -1) 
		self.__link 			= None
		self.Show(True)

	def CreateWidgets(self):
		self.DisplayAlgorithmLabel()
		self.DisplayAlgorithmName()
		self.DisplayAlgorithmKey()
		self.ParametersWidgets()
		self.OutputImagePixelTypeWidgets()
		self.PressFilterButtonWidget()

	def DisplayAlgorithmLabel(self):
		""" Label of the algorithm """
		self.NextRow()
#		import pdb
#		pdb.set_trace()
		key = self.algorithm.GetKey()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Name of the algorithm")
		wx.StaticText(self.panel, -1, pos = (200, self.tempRow), label = self.algorithm.GetLabel())

	def DisplayAlgorithmName(self):
		""" Name of the algorithm """
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "ITK algorithm name")
		wx.StaticText(self.panel, -1, pos = (200, self.tempRow), label = self.algorithm.GetName())

	def DisplayAlgorithmKey(self):
		""" Key of the algorithm """
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Key")
		wx.StaticText(self.panel, -1, pos = (200, self.tempRow), label = self.algorithm.GetKey())
	
	def ParametersWidgets(self):
		""" List of parameters and their text entries """
		self.parameters = self.algorithm.GetParameters()
		self.keys = self.parameters.keys()
		for num in range(self.algorithm.numberParameters):
			# Display label and entry
			textLabel = self.keys[num]
			if self.algorithm.IsTupleParameter(self.keys[num]):
				textLabel += "(separated by commas as x,y)"

			self.NextRow()
			wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = textLabel)
			tempEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))	
			self.parameterEntries.append(tempEntry)

	def PressFilterButtonWidget(self):
		""" Filter Button to continue """
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Press 'Filter' to filter")
		self.filterButton = wx.Button(self.panel, id = 1, label = "Filter", pos = (200, self.tempRow), size = (75, 40))
		self.filterButton.Bind(wx.EVT_BUTTON, self.StartFiltering)
		self.filterButton.SetToolTip(wx.ToolTip("Click 'Filter' to filter"))

	def OutputImagePixelTypeWidgets(self):
		""" Choose if output image is to be 8 or 16 bit, or real type """
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Output Image Pixel Type")
		
		radioUC = wx.RadioButton(self.panel, -1, '8-bit', pos = (200, self.tempRow), style = wx.RB_GROUP)
		radioUS = wx.RadioButton(self.panel, -1, '16-bit', pos = (250, self.tempRow))
		radioF = wx.RadioButton(self.panel, -1, 'float', pos = (300, self.tempRow))

		for eachRadio in [radioUC, radioUS, radioF]:
			self.Bind(wx.EVT_RADIOBUTTON, self.SetOutputImagePixelType, eachRadio)

	def SetOutputImagePixelType(self, event):
		radioSelected = event.GetEventObject()
		if radioSelected.GetLabel() == '8-bit':
			self.algorithm.SetOutputImagePixelType("UC")
		elif radioSelected.GetLabel() == '16-bit':
			self.algorithm.SetOutputImagePixelType("US")
		elif radioSelected.GetLabel() == 'float':
			self.algorithm.SetOutputImagePixelType("F")

					
	def StartFiltering(self, event):
		""" Once filter is pressed, record the algorithm details and return. """
		flag = False
		for num in range(self.algorithm.numberParameters):
			if not self.parameterEntries[num].GetValue():
				print "Please enter all the parameters"
				flag = True
			else:
				 self.algorithm.AssignParameter(
						 self.keys[num],
						 self.parameterEntries[num].GetValue() )
		if not flag:
			self.Destroy()

	def NextRow(self):
		self.tempRow += 30
