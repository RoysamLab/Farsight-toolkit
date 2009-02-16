#!/usr/bin/env python

import os
import pickle
import re
import sys
import xml.parsers.expat

parameters = []
values = {}

def StartElementHandler(name, attributes):
  if name != "parameter":
    return
  parameterName = attributes["name"]
  parameters.append(attributes)

def EndElementHandler(name):
  pass

def CharacterDataHandler(data):
  pass

#this function parses information about parameters from a specified XML file
def ParseParametersFromXML():
  filename = "Parameters.xml"
  try:
    f = file(filename, "r")
  except IOError:
    sys.stderr.write("Error opening Parameters.xml for reading\n")
    sys.exit(1)
  parser = xml.parsers.expat.ParserCreate()
  parser.StartElementHandler = StartElementHandler
  parser.EndElementHandler = EndElementHandler
  parser.CharacterDataHandler = CharacterDataHandler
  parser.ParseFile(f)

#this function asks the user for parameter values based on information parsed
#from the XML file
def AskForValues():
  if os.path.exists("ParameterValues.pck"):
    print "Parameter values have been previously defined."
    print "Would you like to enter new parameter values?"
    answer = sys.stdin.readline()
    if answer.lower().startswith("n"):
      return
  for parameter in parameters:
    printInstructions = True
    while not parameter["name"] in values.keys():
      goodValueSupplied = True
      if printInstructions:
        print "Please enter a value for %s, or ? for more information" % parameter["name"]
        if parameter["type"] == "number":
          print "Typical range: %s to %s, default value: %s" \
            % (parameter["min"], parameter["max"], parameter["default"])
      answer = sys.stdin.readline().strip()
      if answer.lower().startswith("?"):
        print parameter["description"]
        printInstructions = False
        continue
      elif answer == "":
        answer = parameter["default"]

      #check that the user supplied a valid response for this type of parameter
      if parameter["type"] == "file":
        if not os.path.exists(answer.strip()):
          print "%s does not exist.  Please try again.\n" % answer
          goodValueSupplied = False
          printInstructions = False
      elif parameter["type"] == "number":
        try:
          float(answer)
        except ValueError:
          print "%s requires a numerical response.  Please try again.\n" % parameter["name"]
          goodValueSupplied = False
          printInstructions = False
          
      if goodValueSupplied:
        values[parameter["name"]] = answer
        printInstructions = True
  f = file("ParameterValues.pck", "w")
  pickle.dump(values, f)

def SpecifyParameterValues():
  ParseParametersFromXML()
  AskForValues()
