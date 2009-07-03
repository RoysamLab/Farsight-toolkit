#!/usr/bin/env python

import os
import os.path
import subprocess
import sys

#since we're going to be making use of relative paths we need to be sure
#that we're running from the directory that contains this script
scriptPath = os.path.dirname(sys.argv[0])
if scriptPath != "":
  os.chdir(scriptPath)

sys.path.append("Common")
from ParseParameters import ParseXMLParameterFile, SpecifyParameterValues

SpecifyParameterValues("RPITrace3D")
inputImage = ParseXMLParameterFile("RPITrace3D", False)["input_image"]
outputXML = inputImage[0:inputImage.rfind(".")] + "TracedPoints.xml"
print "Here's the outputXML file: " + outputXML

#change to the binary directory before we start calling executables
os.chdir("../bin")
subprocess.call(["./RPITrace3D",\
                 "../python/XML/RPITrace3DParameterValues.xml"])

#display the results from the first step
#some validation that everything went well would probably be smart...
subprocess.call(["./vtkTraceViewer",\
                 outputXML]) 
