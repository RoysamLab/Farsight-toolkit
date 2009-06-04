# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ListOfItems.py
#				- this class was initially meant to create a general purpose listbox, 
#			 	  but later turned into a specific case for algorithm selection. This 
#	              accepts a list of Algorithms, an diplays them in a list box. The user
# 				  has to choose one, and press the Select Key to continue. Help on 
# 				  the pre-processing algorithms is also given in terms of links to 
# 				  Wiki pages.
#
# Author 	: Adarsh K. Ramasubramonian
# Date 		: 31 May, 2009
# 				  
# Last modified on 3 June, 2009.			
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import wx				# Import wxPython library for GUI.	
import ListOfURLs		# For URLs to the algorithm's wiki page.

[wxID_FRAMELISTBOX] = [wx.NewId() for _init_ctrls in range(1)]

class ListOfItems(wx.Frame):
	def __init__(self, parent, id, title, listOfAlgorithms, selection):
		wx.Frame.__init__(self, parent, id, title, size = (500,300))

		self.listOfAlgorithms 	= listOfAlgorithms
		self.__selection 		= selection
		self.__link 			= None
		self.panel 				= wx.Panel(self, -1)
		self.mainURL 			= "http://www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python"

		self.CreateWidgets()		

		self.Show(True)

	def CreateWidgets(self):
		self.DisplayMessage()			# Ask user to select
		self.DisplayListofItems()		# List of Items
		self.DisplayURLForHelp()		# URL for help - to mainpage
		self.CreateSelectButton()		# Select Button
	
	def DisplayURLForHelp(self):
		wx.lib.agw.hyperlink.HyperLinkCtrl(
				self.panel, -1, 
				label = "Help on Pre-processing algorithms", 
				URL = ListOfURLs.URLOfAlgorithmHelp().GetURLOfAlgorithm("mainpage"), pos = (250,20))	

	def CreateSelectButton(self):
		self.button = wx.Button(self.panel, id = -1, 
				label = "Select", pos = (210,110), size = (75,30))
		self.button.Bind(wx.EVT_BUTTON, self.buttonClick)
		self.button.SetToolTip(wx.ToolTip("Click to select parameters"))

	def DisplayMessage(self ):
		self.label = wx.StaticText(self.panel, -1, pos = (10,20), 
				label = "Select an algorithm and press OK to continue")

	def DisplayListofItems(self):
		self.listBox = wx.ListBox(choices=[], id=wxID_FRAMELISTBOX, name = "listBox", 
				parent = self.panel, pos = wx.Point(8,48), 
				size = wx.Size(184,156), style = 0)
		self.listBox.SetBackgroundColour(wx.Colour(255,255,128))

		for name in self.listOfAlgorithms:
			self.listBox.Append(name)
		self.listBox.Bind(wx.EVT_LISTBOX, self.OnListBoxListbox, id = wxID_FRAMELISTBOX)

	# When an algorithm is chosen, we will display a URL to the algorithm's wiki page.
	def OnListBoxListbox(self, event):
		selName = self.listBox.GetStringSelection()
		if self.__link:
			self.__link.Destroy()
		self.__link = wx.lib.agw.hyperlink.HyperLinkCtrl(
				self.panel, -1, label = "Help on " + selName + " Filter", 
				URL = ListOfURLs.URLOfAlgorithmHelp().GetURLOfAlgorithm(selName), pos = (250,40))

		# XXX Include a brief description of the algorithm

	# Algorithm chosen, proceed to get the parameters.
	def buttonClick(self, event):					
		selectedAlgorithm = self.listBox.GetStringSelection()
		if selectedAlgorithm:
			self.__selection.append(selectedAlgorithm)
			print "You chose the algorithm ", selectedAlgorithm	
			self.Destroy()
		else:
			print "Please choose an algorithm"

	def GetSelection(self):
		return self.__selection

# - - - - - - - - - - - -  End of ListOfItems - - - - - - - - - - - - - - - 

# Create a frame and display list of algorithms.
def DisplayListOfAlgorithms(listOfAlgorithms):
	root = wx.App(False)
	selection = []
	frame = ListOfItems(None, -1, "Choose the algorithm", listOfAlgorithms, selection)
	root.MainLoop()
	return selection[0]	

				   
