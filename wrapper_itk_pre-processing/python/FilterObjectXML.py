# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# FilterObjectXML.py

# This module contains a class that is involved in creating an XML file from the 
# FilterObject/Algorithm object that is passed to it as a paramter. The document 
# created is then written to an XML file.

# Author 	: Adarsh K. Ramasubramonian 
# Date 		: 22 May, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
import xml.dom.minidom as dom
import pdb


class ObjectXML(object):
    def __init__(self):
        self.__doc = dom.Document()

    def WriteDocument(self, filename):		# Write to a file
        fileObject = open(filename, "w")
        self.__doc.writexml(fileObject)
        fileObject.close()

    def ReadDocument(self, filename):		# Read from a file
        fileObject = open(filename, "r")
        self.__doc = dom.parse(fileObject)
        fileObject.close()

    def PrintToScreen(self):				# Print to screen
        print(self.__doc.toprettyxml(indent='\t',newl='\n'))

    def NewAttributeAndSet(self, nodeObject, attributeName, value):		# Create an attribute and add it to the element
        nodeObject.setAttributeNode(self.__doc.createAttribute(attributeName))
        nodeObject.setAttribute(attributeName, value)

    def GetFirstChild(self):
        return self.__doc.firstChild

    def GetNewElement(self, namespace, name):
        return self.__doc.createElementNS(namespace, name)

    def GetNewAttribute(self, key):
        return self.__doc.createAttribute(key)

    def AppendChildToDoc(self, child):
        self.__doc.appendChild(child)

class FilterObjectXML(ObjectXML):
	# We define a set of attributes of the FilterAlgorithm class. Note that we have 
	# split them into two sets because some of the them are strings and the other are 
	# int/float/bool type objects. When parsing the document, we would have to assign 
	# the algorithm with the help of "eval" functions, and this separation helps.

	listOfStringAttributes = ["inputFileName", 
							  "outputFileName",
							  "outputImagePixelType",
							  "dimension",
							  ]

	listOfOtherAttributes = ["hasStructuringElement",
							 "inputImageStack",
							 "outputImageStack",
							 "startIndex",
							 "endIndex",
							 "incrementIndex",
							 "outputType",
							 "useRescaler"
							 ]


	listOfAttributes = listOfStringAttributes + listOfOtherAttributes

	def __init__(self):								# Constructor
        	super(FilterObjectXML, self).__init__()

	def UpdateDocument(self, algorithm, filename):	# Write the algorithm details to the XML file, the argument
													# algorithm is a FilterAlgorithm object.  

		# Add name
		xmlalgorithm = self.GetNewElement("Algorithm","xmlalgorithm")
		self.AppendChildToDoc(xmlalgorithm)	

		# Assign attributes
		self.AssignAllAttributes(FilterObjectXML.listOfAttributes, xmlalgorithm, algorithm)

		# Create an element filter - and add its attributes
		filter = self.GetNewElement("Algorithm","filter")
		xmlalgorithm.appendChild(filter)
		self.AssignAllAttributes(["key","name","label"], filter, algorithm)

#		self.PrintToScreen()

		# Create an element parameters - and add its attributes
		parameters = filter.appendChild(self.GetNewElement("Algorithm","parameters"))
		algorithmParameters = algorithm.GetParameters()
		keys = algorithmParameters.keys()
		for key in keys:
			# The tuple parameters (that need to an object to be defined) have an additional 
			# attribute "type". For example, a radius object will be an object of class itkSize2()
			# for a 2-D image.
			if algorithm.IsTupleParameter(key):
				tempObject = self.GetNewElement("Algorithm",key)
				parameters.appendChild(tempObject)
				self.NewAttributeAndSet( tempObject, "type", algorithmParameters[key]["type"] )
				self.NewAttributeAndSet( tempObject, key, str(algorithmParameters[key]["value"]) )
													# The value of attributes should be strings.
			else:
				self.NewAttributeAndSet( parameters, key, str(algorithmParameters[key]["value"]) )


#		self.PrintToScreen()
		self.WriteDocument(filename)


	# Read the document and return the algorithm object.
	def ParseDocument(self, algorithm, filename):

		self.ReadDocument(filename)
#		self.PrintToScreen()
		xmlalgorithm = self.GetFirstChild()

		# Get all the attributes of the algorithm.
		self.GetAllStringAttributes(FilterObjectXML.listOfStringAttributes, xmlalgorithm, algorithm)
		self.GetAllOtherAttributes(FilterObjectXML.listOfOtherAttributes, xmlalgorithm, algorithm)

		# Get the filter element
		filterList = xmlalgorithm.getElementsByTagName("filter")
		for filter in filterList:
			self.GetAllStringAttributes(["key", "name"], filter, algorithm)

		# Get the parameter element
		parameterlist = filter.getElementsByTagName("parameters")
		for parameters in parameterlist:
			attributesOfParameters = parameters.attributes
			keys = attributesOfParameters.keys()	

			for key in keys:
				attr = attributesOfParameters[key]
				algorithm.AddParameter(attr.localName)
				algorithm.AssignParameter(attr.localName , attr.nodeValue)
	
			for child in parameters.childNodes:
				algorithm.AddTupleParameter(child.localName, len(child.getAttribute(child.localName)), child.getAttribute("type"))
				algorithm.AssignParameter(child.localName, child.getAttribute(child.localName))
			
		
	# For each attribute in the list, invoke the corresponding Get() function and update the attribute.
	def AssignAllAttributes(self, listOfAttributes, xmlalgorithm, thisalgorithm):
		for name in listOfAttributes:
			titleName = self.CapitalizeFirstLettter(name)
			self.NewAttributeAndSet(
					xmlalgorithm,
					name,
					str( eval( "thisalgorithm.Get" + titleName + "()" ) ) )

	def CapitalizeFirstLettter(self, str):
		newString = str[0].upper() + str[1:] 
		return newString

	# For each attribute in the list, invoke the corresponding Set() function - these are not related to ITK.
	def GetAllStringAttributes(self, listOfAttributes, element, thisalgorithm):
		for name in listOfAttributes:
			titleName = self.CapitalizeFirstLettter(name)
			eval( "thisalgorithm.Set" + titleName + "( element.getAttribute(\"" + name + "\"))" )
	
	def GetAllOtherAttributes(self, listOfAttributes, element, thisalgorithm):
		for name in listOfAttributes:
			titleName = self.CapitalizeFirstLettter(name)
			eval( "thisalgorithm.Set" + titleName + "( eval( element.getAttribute(\"" + name + "\") ) )" )
	
		
