#!/usr/bin/env python

import os
import pickle
import re
import sys
import xml.parsers.expat

parameters = []
values = {}

#XML handler for reading parameter description files
def StartParameterDescriptionHandler(name, attributes):
  if name != "parameter":
    return
  parameterName = attributes["name"]
  parameters.append(attributes)

#XML handler for reading parameter value files
def StartParameterValueHandler(name, attributes):
  if name != "parameter":
    return
  values[attributes["name"]] = attributes["value"]

def EndElementHandler(name):
  pass

def CharacterDataHandler(data):
  pass

#this function parses information about parameters from a specified XML file
#if descriptions are true it reads from the parameter description file
#otherwise it parses a parameter values file.
def ParseXMLParameterFile(moduleName, descriptions=True):
  parser = xml.parsers.expat.ParserCreate()
  if descriptions:
    filename = "%sParameterDescriptions.xml" % moduleName
    parser.StartElementHandler = StartParameterDescriptionHandler
  else:
    filename= "%sParameterValues.xml" % moduleName
    parser.StartElementHandler = StartParameterValueHandler
  try:
    f = file(filename, "r")
  except IOError:
    sys.stderr.write("Error opening %s for reading\n" % filename)
    sys.exit(1)
  parser.EndElementHandler = EndElementHandler
  parser.CharacterDataHandler = CharacterDataHandler
  parser.ParseFile(f)
  if not descriptions:
    return values

#this function asks the user for parameter values based on information parsed
#from the XML file
def AskForValues(moduleName):
  if os.path.exists("%sParameterValues.xml" % moduleName):
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
  f = file("%sParameterValues.xml" % moduleName, "w")
  f.write('<?xml version="1.0"?>\n')
  f.write('<parameters>\n')
  for name,value in values.iteritems():
    f.write('<parameter name="%s" value="%s"></parameter>\n' % (name,value))
  f.write('</parameters>\n')

def SpecifyParameterValues(moduleName):
  ParseXMLParameterFile(moduleName, True)
  AskForValues(moduleName)
