# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# ListOfAlgorithmsWithURL.py - 
#		- this class is used for providing a list of algorithms and their associated URL.s The class is derived from the base class ListOfItems.
#
# Author 	: Adarsh K. Ramasubramonian
# Date 		: 30 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import wx
import ListOfItems 
import AlgorithmXML
import HyperLinkModule

class ListOfAlgorithmsWithURL(ListOfItems.ListOfItems):
	def __init__(self, parent, id, title, listOfAlgorithms, selection):
		self.mainURL 			= "http://www.farsight-toolkit.org/wiki/ITK_Pre-Processing_Algorithm_Wrappers_in_Python"
		self.__helplink 		= None
		self.__advancedhelplink = None
		super(ListOfAlgorithmsWithURL, self).__init__(parent, id, title, listOfAlgorithms, selection)

	def CreateWidgets(self):
		super(ListOfAlgorithmsWithURL, self).CreateWidgets()
		self.DisplayURLForHelp()		# URL for help - to mainpage

	def DisplayURLForHelp(self):
		HyperLinkModule.hyperlink.HyperLinkCtrl(
				self.panel, -1, 
				label = "Help on Pre-processing algorithms", 
				URL = self.mainURL,
                pos = (250,20))	

	# When an algorithm is chosen, we will display a URL to the algorithm's wiki page.
	def OnListBoxSelect(self, event):
		selName = self.listBox.GetStringSelection()
		xmlObject = AlgorithmXML.AlgorithmXML()                                 
		url1, url2 = xmlObject.GetURLOfAlgorithm(selName.lower())

		if self.__helplink:
			self.__helplink.Destroy()
		self.__helplink = HyperLinkModule.hyperlink.HyperLinkCtrl(
				self.panel, -1, label = "Help on " + selName + " Filter", 
                URL = url1, pos = (250, 40))

		if self.__advancedhelplink:
			self.__advancedhelplink.Destroy()
		self.__advancedhelplink = HyperLinkModule.hyperlink.HyperLinkCtrl(
				self.panel, -1, label = "Advanced help on " + selName + " Filter", 
                URL = url2, pos = (250, 70))
# - - - - - - - - - - - -  End of ListOfAlgorithmsWithURL - - - - - - - - - - - - - - - 

# Create a frame and display list of algorithms.
def DisplayList(listOfAlgorithms):
	root = wx.App(False)
	selection = []
	frame = ListOfAlgorithmsWithURL(None, -1, "Choose the algorithm", listOfAlgorithms, selection)
	root.MainLoop()
	return selection[0]	

				   
