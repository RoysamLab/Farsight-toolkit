# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
# EnterAlgorithm.py
#
# This module contains all the GUI required for inserting a particular algorithm in the list of algorithms available,
# with a user interface written in wxPython.
#
# Author 	: Adarsh K. Ramasubramonian
# Date 		: 25 June, 2009
#
# Last modified on: 30 June 2009.
# Last modified on: 24 July 2009.
# 		- when the number of parameters of an existing algorithm is changed, the parameter entries are not displayed properly.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import wx
import FilterAlgorithm
import pdb
import sys
import AlgorithmXML

class EnterAlgorithmGUIApp(wx.App):
	def __init__(self, state):
		""" Constructor for EnterAlgorithmGUIApp """
		wx.App.__init__(self, state)
		app = EnterAlgorithmGUIwx(None, -1, "Algorithm", (600,600))
		app.CreateWidgets()
		app.Show(True)

class EnterAlgorithmGUIwx(wx.Frame):
	""" GUI for entering the information about an algorithm """
	def __init__(self, parent, id, title, size):
		wx.Frame.__init__(self, parent, id, title, size = (600, 600))
#pdb.set_trace()
		self.algorithmNumberParametersOldValue = 0
		self.algorithmNumberParameters = 0
		self.hasStructuringElement = 0
		self.tempRow = 0
		self.panel = wx.Panel(self, -1)
		self.algorithmParameterEntry = []

	def CreateWidgets(self):
		self.GetAlgorithmName()
		self.DisplaySubmitButton()
		self.GetAlgorithmKey()
		self.DisplayQuitButton()							  
		self.GetAlgorithmLabel()															  
		self.DisplayDisplayExistingAlgorithmButton()								  
		self.GetHasStructuringElement()
		self.GetHelpURL()
		self.GetAdvancedHelpURL()
		self.GetOutputType()				 
		self.GetRescalerWidget()
		self.GetNumberOfParameters()
		self.DisplayEnterButton()
	
	def GetAlgorithmName(self):
		# Algorithm name - label and box
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "ITK Algorithm Name")
		self.algorithmNameEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
		
	def GetAlgorithmKey(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Algorithm Key")
		self.algorithmKeyEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))

	def GetAlgorithmLabel(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Algorithm Name")
		self.algorithmLabelEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))																						  
	def GetHasStructuringElement(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Has structuring element?")
		self.hasStructuringElementCheckBox = wx.CheckBox(self.panel, -1, "", pos = (200, self.tempRow))
	
	def GetHelpURL(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Wiki URL for help")
		self.helpURLEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
	
	def GetAdvancedHelpURL(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "URL for Advanced help")
		self.advancedHelpURLEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))

	def GetOutputType(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "If no separate output type, select.")
		self.notOutputTypeCheckBox = wx.CheckBox(self.panel, 1, pos = (200, self.tempRow))

	def GetRescalerWidget(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Don't use Rescaler (rescale intensities)?")
		self.notUseRescalerCheckBox = wx.CheckBox(self.panel, 1, pos = (200, self.tempRow))

	def GetNumberOfParameters(self):
		self.NextRow()
		wx.StaticText(self.panel, -1, pos = (20, self.tempRow), label = "Number of Parameters")
		self.numberOfParametersEntry = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
		self.NextRow()																				
		wx.StaticText(self.panel, 1, pos = (20, self.tempRow), label = "Name of parameter")
		wx.StaticText(self.panel, 1, pos = (120, self.tempRow), label = "Tuple?")
		wx.StaticText(self.panel, 1, pos = (200, self.tempRow), label = "Length")
		wx.StaticText(self.panel, 1, pos = (300, self.tempRow), label = "Type")

		self.lastEntryRow = self.tempRow

	def DisplayEnterButton(self):
		self.numberOfParamtersButton = wx.Button(self.panel, id = 1, label = "Enter", pos = (400, self.tempRow), size = (75,20))
		self.numberOfParamtersButton.Bind(wx.EVT_BUTTON, self.DisplayEntriesForParameters)
		self.numberOfParamtersButton.SetToolTip(wx.ToolTip("Press enter after typing number of parameters"))

	def DisplaySubmitButton(self):
		self.submitButton = wx.Button(self.panel, id = 1, label = "Submit", pos = (400, self.tempRow))
		self.submitButton.Bind(wx.EVT_BUTTON, self.AlgorithmSubmit)
		self.submitButton.SetToolTip(wx.ToolTip("Press submit to store the details of the algorithm"))
	
	def DisplayQuitButton(self):
		self.quitButton = wx.Button(self.panel, id = 1, label = "Quit", pos = (400, self.tempRow))
		self.quitButton.Bind(wx.EVT_BUTTON, self.Quit)
		self.quitButton.SetToolTip(wx.ToolTip("Press Quit to exit"))

	def DisplayDisplayExistingAlgorithmButton(self):
		self.displayButton = wx.Button(self.panel, id = 1, label = "Display existing", pos = (400,self.tempRow))
		self.displayButton.Bind(wx.EVT_BUTTON, self.DisplayExistingAlgorithm)
		self.displayButton.SetToolTip(wx.ToolTip("Display details if algorithm exists already"))

	def Quit(self, event):
		self.Destroy()

	# For parameter, by default, we will allow only one entry. However, if the checkbox is checked, then we will allow parameters to be set for this parameter also
	def DisplayEntriesForParameters(self, event):
#		pdb.set_trace()
		self.tempRow = self.lastEntryRow		# May not be necessary
		oldValue = self.algorithmNumberParameters
		newValue = int(self.numberOfParametersEntry.GetValue())
		if newValue < 0:
			print "Please enter positive value for number of parameters"
			return

#		print "Number of parameters are ", self.algorithmNumberParameters
		differenceValue = newValue - oldValue
		if not differenceValue:
			""" No change. """
		elif differenceValue > 0:
			""" Add entries """
			for i in range(differenceValue):
				self.algorithmParameterEntry.append(self.AddNewParameterRow())
		else:
			""" Delete entries """
			for i in range(-differenceValue):
				self.DeleteParameterRow(self.algorithmParameterEntry[newValue])
				del self.algorithmParameterEntry[newValue]

		self.algorithmNumberParameters = newValue
		self.lastEntryRow = self.tempRow		# May not be necessary

	def AddNewParameterRow(self):
		entry = {}
		self.NextRow()
		entry["nameEntry"] = wx.TextCtrl(self.panel, 1, pos = (20, self.tempRow))
		entry["isTupleParameter"] = wx.CheckBox(self.panel, 1,  pos = (150, self.tempRow))
		entry["isTupleParameter"].SetValue(wx.CHK_UNCHECKED)
		entry["length"] = wx.TextCtrl(self.panel, 1, pos = (200, self.tempRow))
		entry["type"] = wx.TextCtrl(self.panel, 1, pos = (300, self.tempRow))
		entry["length"].Enable(False)
		entry["type"].Enable(False)
		self.Bind(wx.EVT_CHECKBOX, self.EnableEntries, entry["isTupleParameter"])
		return entry
		
	def DeleteParameterRow(self, entry):
#pdb.set_trace()
		entry["nameEntry"].Destroy()
		entry["isTupleParameter"].Destroy()
		entry["length"].Destroy()
		entry["type"].Destroy()
		self.PreviousRow()

	def EnableEntries(self, event):
#pdb.set_trace
		flag = 0
		for entry in self.algorithmParameterEntry:
			if event.GetEventObject() == entry["isTupleParameter"]:
				flag = 1;
				if entry["isTupleParameter"].GetValue():
					entry["length"].Enable(True)
					entry["type"].Enable(True)
			  	else:
					entry["length"].Enable(False)
					entry["type"].Enable(False)
		if not flag:
			print "Checkbox object not found. Please check !!!"
			sys.exit()

	def NextRow(self):
		self.tempRow += 30
	def PreviousRow(self):
		self.tempRow -= 30

	def AlgorithmSubmit(self, event):		
		print "Submitting algorithm details"
		# Summarizing algorithm
		thisAlgorithm = FilterAlgorithm.FilterAlgorithm()
		thisAlgorithm.SetName(self.algorithmNameEntry.GetValue())
		thisAlgorithm.SetKey(self.algorithmKeyEntry.GetValue())
		thisAlgorithm.SetLabel(self.algorithmLabelEntry.GetValue())
		for entryNumber in range(self.algorithmNumberParameters):
			if self.algorithmParameterEntry[entryNumber]["isTupleParameter"].GetValue():
				thisAlgorithm.AddTupleParameter(
						self.algorithmParameterEntry[entryNumber]["nameEntry"].GetValue(),
						self.algorithmParameterEntry[entryNumber]["length"].GetValue(),
						self.algorithmParameterEntry[entryNumber]["type"].GetValue()
						)
																				
			else:																						
				thisAlgorithm.AddParameter(self.algorithmParameterEntry[entryNumber]["nameEntry"].GetValue())
		if self.notOutputTypeCheckBox.GetValue():
			thisAlgorithm.SetOutputType(None)

		thisAlgorithm.SetHasStructuringElement(
				self.hasStructuringElementCheckBox.GetValue())
		thisAlgorithm.SetUseRescaler(
				not self.notUseRescalerCheckBox.GetValue())

		thisAlgorithm.SetHelpURL(self.helpURLEntry.GetValue())
		thisAlgorithm.SetAdvancedHelpURL(self.advancedHelpURLEntry.GetValue())
		
		# Now we need to update the file SmoothingAlgorithms.xml
#pdb.set_trace()
		xmlObject = AlgorithmXML.AlgorithmXML()
		# Check if algorithm already exists
		if xmlObject.CheckIfAlgorithmExists(thisAlgorithm):
			if self.AskUserIfReplace():
				xmlObject.AddAlgorithm(thisAlgorithm,True)
			else:
				print "Algorithm not added"
				return
		else:
			xmlObject.AddAlgorithm(thisAlgorithm, False)

	def AskUserIfReplace(self):
		dialogBox = wx.MessageDialog(None, 'Algorithm already exists. Replace?', '', wx.YES_NO| wx.YES_DEFAULT| wx.ICON_QUESTION)
		answer = dialogBox.ShowModal()
		if answer == wx.ID_YES:
			return True
		else: 
			return False

	def DisplayExistingAlgorithm(self, event):
		key = self.algorithmKeyEntry.GetValue()
		xmlObject = AlgorithmXML.AlgorithmXML()
		algorithm = xmlObject.GetAlgorithm(key)
#pdb.set_trace()
		if algorithm:	
			self.algorithmNameEntry.SetValue(algorithm.GetName()) 
			self.algorithmLabelEntry.SetValue(algorithm.GetLabel())
	
			if algorithm.GetHasStructuringElement():
				self.hasStructuringElementCheckBox.SetValue(wx.CHK_CHECKED)
			else:
				self.hasStructuringElementCheckBox.SetValue(wx.CHK_UNCHECKED)


			if not algorithm.GetOutputType():
				self.notOutputTypeCheckBox.SetValue(wx.CHK_CHECKED)
			else:
				self.notOutputTypeCheckBox.SetValue(wx.CHK_UNCHECKED)

			if not algorithm.GetUseRescaler():
				self.notUseRescalerCheckBox.SetValue(wx.CHK_CHECKED)
			else:
				self.notUseRescalerCheckBox.SetValue(wx.CHK_UNCHECKED)
	
			# Delete existing parameters			
			for entry in self.algorithmParameterEntry:
				self.DeleteParameterRow(entry)
			numberOfItems = len(self.algorithmParameterEntry)
			for i in range(numberOfItems):
				self.algorithmParameterEntry.pop()
			parameters = algorithm.GetParameters()
			self.algorithmNumberParameters = algorithm.numberParameters
			self.numberOfParametersEntry.SetValue(str(self.algorithmNumberParameters))
	
			keys = parameters.keys()
			for key in keys:
				entry = self.AddNewParameterRow()
				self.algorithmParameterEntry.append(entry)
				entry["nameEntry"].SetValue(key)
				if algorithm.IsTupleParameter(key):
					entry["isTupleParameter"].SetValue(wx.CHK_CHECKED)
					entry["length"].SetValue(str(2))		# may be obsolete
					entry["type"].SetValue(parameters[key]["type"])
					entry["length"].Enable(True)													   
					entry["type"].Enable(True)													   
	
	
			self.helpURLEntry.SetValue(algorithm.GetHelpURL())
			self.advancedHelpURLEntry.SetValue(algorithm.GetAdvancedHelpURL())
			self.lastEntryRow = self.tempRow
		

# main

root = EnterAlgorithmGUIApp(0)
root.MainLoop()



		

