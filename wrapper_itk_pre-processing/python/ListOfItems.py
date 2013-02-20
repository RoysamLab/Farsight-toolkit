# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ListOfItems.py
#				- this class creates a general purpose listbox. The user
# 				  has to choose one, and press the Select Key to continue. #

# Author 	: Adarsh K. Ramasubramonian
# Date 		: 31 May, 2009
# 				  
# Last modified on 3 June, 2009.			
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import wx				# Import wxPython library for GUI.	
import AlgorithmXML

[wxID_FRAMELISTBOX] = [wx.NewId() for _init_ctrls in range(1)]

class ListOfItems(wx.Frame):
	def __init__(self, parent, id, title, listOfAlgorithms, selection):
		wx.Frame.__init__(self, parent, id, title, size = (500,300))

		self.listOfAlgorithms 	= listOfAlgorithms
		self.__selection 		= selection
		self.panel 				= wx.Panel(self, -1)

		self.CreateWidgets()		
		self.Show(True)

	def CreateWidgets(self):
		self.DisplayMessage()			# Ask user to select
		self.DisplayListofItems()		# List of Items
		self.CreateSelectButton()		# Select Button
	

	def CreateSelectButton(self):
		self.button = wx.Button(self.panel, id = -1, 
				label = "Select", pos = (210,110), size = (75,30))
		self.button.Bind(wx.EVT_BUTTON, self.ButtonClick)
		self.button.SetToolTip(wx.ToolTip("Click to select parameters"))

	def DisplayMessage(self ):
		self.label = wx.StaticText(self.panel, -1, pos = (10,20), 
				label = "Select an algorithm and press Select to continue")

	def DisplayListofItems(self):
		self.listBox = wx.ListBox(choices=[], id=wxID_FRAMELISTBOX, name = "listBox", 
				parent = self.panel, pos = wx.Point(8,48), 
				size = wx.Size(184,156), style = 0)
		self.listBox.SetBackgroundColour(wx.Colour(255,255,128))

		for name in self.listOfAlgorithms:
			self.listBox.Append(name)
		self.listBox.Bind(wx.EVT_LISTBOX, self.OnListBoxSelect, id = wxID_FRAMELISTBOX)

	def OnListBoxSelect(self, event):
		# NOT IMPLEMENTED
		return

	# Algorithm chosen, proceed to get the parameters.
	def ButtonClick(self, event):					
		selectedAlgorithm = self.listBox.GetStringSelection()
		if selectedAlgorithm:
			self.__selection.append(selectedAlgorithm)
			print "You chose the algorithm ", selectedAlgorithm	
			self.Destroy()
		else:
			print "Please choose an algorithm"

	def GetSelection(self):
		return self.__selection

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 		

