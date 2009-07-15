# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This file contains the class that will write any algorithm into an XML file that contains the list of algorithms. This will be called from the EnterAlgorithm.py function. 

# Author    : Adarsh K. Ramasubramonian
# Date      : 25 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

import xml.dom.minidom as dom
import os.path
import FilterAlgorithm
from FilterObjectXML import *
import sys

class AlgorithmXML(ObjectXML):
    def __init__(self):
        super(AlgorithmXML, self).__init__()
        if os.path.exists("SmoothingAlgorithms.xml"):
            self.ReadDocument("SmoothingAlgorithms.xml")
            self.AllAlgorithms = self.GetFirstChild()
                                               
        else:
            self.AllAlgorithms = self.GetNewElement("Algorithm","algorithms")
            self.AppendChildToDoc(self.AllAlgorithms)
         

    def CheckIfAlgorithmExists(self, algorithm):
        children = self.AllAlgorithms.childNodes
        returnChild = None
        replaceChild = False                          
        for child in children:
            if algorithm.GetKey() == child.localName:
                while not returnChild:
                    character = raw_input("Algorithm exists already. Replace? (Y/N)")  
                    if character.capitalize() == "Y":
                        returnChild = child
                        replaceChild = True
                    elif character.capitalize() == "N":
                        returnChild = child
                        replaceChild = False
                    else:
                        print "Press Y/N."
        return returnChild, replaceChild                        


    def AddAlgorithm(self, algorithm):
        newAlgorithm = self.GetNewElement("Algorithm",algorithm.GetKey())
        child, replaceChild = self.CheckIfAlgorithmExists(algorithm)                   
        if not child:
            self.AllAlgorithms.appendChild(newAlgorithm)
        elif replaceChild:
            self.AllAlgorithms.replaceChild(newAlgorithm, child)
        else:
            return

        self.NewAttributeAndSet(newAlgorithm, "name", algorithm.GetName())
        self.NewAttributeAndSet(newAlgorithm, "key", algorithm.GetKey())
        self.NewAttributeAndSet(newAlgorithm, "label", algorithm.GetLabel())
        self.NewAttributeAndSet(newAlgorithm, "hasStructuringElement", str(algorithm.GetHasStructuringElement()))
        self.NewAttributeAndSet(newAlgorithm, "outputType", str(algorithm.GetOutputType()))
        self.NewAttributeAndSet(newAlgorithm, "useCaster", str(algorithm.GetUseCaster()))
		
        parameters = self.GetNewElement("Algorithm","parameters")
        newAlgorithm.appendChild(parameters)

        algorithmParameters = algorithm.GetParameters()
        keys = algorithmParameters.keys()
        for key in keys:
            if algorithm.IsTupleParameter(key):
                tempObject = self.GetNewElement("Algorithm",key)
                parameters.appendChild(tempObject)
                self.NewAttributeAndSet(tempObject, "type", algorithmParameters[key]["type"])
            else:
                parameters.setAttributeNode(self.GetNewAttribute(key))

        self.NewAttributeAndSet(newAlgorithm, "helpURL", algorithm.GetHelpURL())          
        self.NewAttributeAndSet(newAlgorithm, "advancedHelpURL", algorithm.GetAdvancedHelpURL())
        # XXX May need to replace child here.
        self.AllAlgorithms.appendChild(newAlgorithm)
        self.WriteDocument("SmoothingAlgorithms.xml")

    def GetAlgorithm(self, key):
		algorithm = FilterAlgorithm.FilterAlgorithm()

		children = self.AllAlgorithms.childNodes
		for child in children:
			if key == child.nodeName:
				""" Found the algorithm. Get the parameters """
				self.GetAllStringAttributes(["name","key","label","helpURL","advancedHelpURL"], child, algorithm)
				self.GetAllOtherAttributes(["hasStructuringElement","outputType","useCaster"], child, algorithm)

				parametersList = child.getElementsByTagName("parameters")
				for parameters in parametersList:
					attributes = parameters.attributes
					keys = attributes.keys()
					for key in keys:	
						attr = attributes[key]
						algorithm.AddParameter(attr.localName)
                        
					for tupleparameters in parameters.childNodes:
						algorithm.AddTupleParameter(tupleparameters.localName, tupleparameters.getAttribute("length"), tupleparameters.getAttribute("type"))
				return algorithm						
		return None
    def GetURLOfAlgorithm(self, key):
        """ We need to return the URLs for the algorithm whose key is passed as a parameter. """
        children = self.AllAlgorithms.childNodes
        for child in children:
            if key == child.nodeName:
                return child.getAttribute("helpURL"), child.getAttribute("advancedHelpURL")
                
        # key not present - something is wrong.
        print "Key ", key," not present in the XML file. Please check."
        sys.exit()

    def CapitalizeFirstLettter(self, str):
        newString = str[0].upper() + str[1:] 
        return newString
	
    def GetListOfAlgorithms(self):
	  listOfAlgorithms = []
#pdb.set_trace()
	  children = self.AllAlgorithms.childNodes
	  for child in children:
	      listOfAlgorithms.append(child.getAttribute("label"))
	  return listOfAlgorithms																

	# For each attribute in the list, invoke the corresponding Set() function.
    def GetAllStringAttributes(self, listOfAttributes, element, thisalgorithm):
        for name in listOfAttributes:
            titleName = self.CapitalizeFirstLettter(name)
            eval( "thisalgorithm.Set" + titleName + "( element.getAttribute(\"" + name + "\"))" )
		
    def GetAllOtherAttributes(self, listOfAttributes, element, thisalgorithm):
        for name in listOfAttributes:
            titleName = self.CapitalizeFirstLettter(name)
            eval( "thisalgorithm.Set" + titleName + "( eval( element.getAttribute(\"" + name + "\") ) )" )
