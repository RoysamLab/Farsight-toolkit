# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This class contains the definition of the class filter object that will contain 
# all the member functions and variables required for a filter object. 

# Author 	: Adarsh K. Ramasubramonian
# Date		: 14 May, 2009

# Last modified on 3 June, 2009
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


class FilterObject(object):
	""" FilterObject - for storing details about algorithms """
	# Members
	def __init__(self):
		""" FilterObject Constructor """
		self.__key = None
		self.__name = None
		self.__label = None						  
		self.__parameters = {}
		self.__tupleParameters = []

	def __str__(self):	
		returnValue = "\nAlgorithm            : " + str(self.GetName())
		returnValue += "\nKey                  : " + str(self.GetKey())
		returnValue += "\nLabel                : " + str(self.GetLabel())

		returnValue += "\nNumber of parameters : " + str(self.numberParameters)
		keys = self.GetParameters().keys()
		for num in range(self.numberParameters):
			returnValue += "\n" + keys[num] + "         : " + str(self.GetParameters()[keys[num]]["value"])
		return returnValue

	# Public functions for getting the member attributes
	def GetNumberParameters(self):
		return len(self.__parameters)

	def GetKey(self):
		return self.__key

	def GetName(self):
		return self.__name

	def GetLabel(self):
		return self.__label

	def GetParameters(self):
		return self.__parameters

	def GetValue(self, keyOfParameter):
		return self.__parameters[keyOfParameter]

	# Public function for setting the member attributes
	def SetKey(self, nameOfKey):
		self.__key = nameOfKey
	
	def SetName(self, nameOfName):
		self.__name = nameOfName

	def SetLabel(self, label):
		self.__label = label

	def SetParameters(self, listOfParameters):	# Actually a dictionary
		self.__parameters = listOfParameters

	def AddParameter(self, keyOfParameter):
		self.__parameters[keyOfParameter] = {}

	def AddParameters(self, keyOfParameter):
		for key in keyOfParameter:
			self.__parameters[str(key)] = {}

	def AddTupleParameter(self, keyOfParameter, tupleSize, elementType):
		self.__parameters[keyOfParameter] = {}
		self.__tupleParameters.append(keyOfParameter)
		self.__parameters[keyOfParameter]["type"] = elementType

	def AssignParameter(self, keyOfParameter, value, dataType = "int"):
		if keyOfParameter not in self.__parameters.keys():
			print "Parameter ", keyOfParameter, " not found. Check again!"
			return
		self.__parameters[keyOfParameter]["value"] = value
	
	numberParameters = property(GetNumberParameters)

	def IsTupleParameter(self, keyOfParameter):
		if keyOfParameter in self.__tupleParameters:
			return True
		else:
			return False


	

	
