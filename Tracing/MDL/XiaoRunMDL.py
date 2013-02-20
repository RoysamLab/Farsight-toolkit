#!/usr/bin/env python

import os.path
import re
import subprocess
import sys
import time

parameterValues = {}

#since we're going to be making use of relative paths we need to be sure
#that we're running from the directory that contains this script
scriptPath = os.path.dirname(sys.argv[0])
if scriptPath != "":
  os.chdir(scriptPath)

sys.path.append("Common")
from ParseParameters import ParseXMLParameterFile, SpecifyParameterValues

SpecifyParameterValues("MDL")
f = file("MDLParameterValues.xml")

parameterValues = ParseXMLParameterFile("MDL", False)

fullPathToInput = parameterValues["Input image"]
dataDir = os.path.dirname(fullPathToInput)
filename = os.path.basename(fullPathToInput)

match = re.compile(r".*?\.(\d+)x(\d+)x(\d+)\.").search(filename)
sizeX = match.group(1)
sizeY = match.group(2)
sizeZ = match.group(3)
intensityThreshold = parameterValues["Intensity threshold"]
connectedComponentsSize = parameterValues["Connected components size"]

#change to the binary directory before we start calling executables
##os.chdir("${CMAKE_BINARY_DIR}")
os.chdir("D:/MDL102/bin/release")

cmd = "./volumeProcess.exe %s %s %s %s 1 %s %s" % \
  (fullPathToInput, sizeX, sizeY, sizeZ, dataDir+"/volume_Processed.raw",\
   intensityThreshold)
print "Step #1, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./volumeProcess",\
                  fullPathToInput, sizeX, sizeY, sizeZ, \
                  dataDir+"/volume_Processed.raw", intensityThreshold])
timeElapsed = time.time() - timeBefore
print "Step #1 completed in %d seconds\n\n" % timeElapsed

cmd = "./ConnCompntwFldfill %s %s %s %s %s %s" % \
  (dataDir+"/volume_Processed.raw", sizeX, sizeY, sizeZ,\
   dataDir+"/components_Connected.raw", connectedComponentsSize)
print "Step #2, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./ConnCompntwFldfill",\
                  dataDir+"/volume_Processed.raw", sizeX, sizeY, sizeZ,\
                  dataDir+"/components_Connected.raw",\
                  connectedComponentsSize])
timeElapsed = time.time() - timeBefore
print "Step #2 completed in %d seconds\n\n" % timeElapsed

cmd = "./AnisoDiffuse %s %s %s %s %s" %\
  (dataDir+"/components_Connected.raw", sizeX, sizeY, sizeZ,\
   dataDir+"/Aniso_Diffused.raw")
print "Step #3, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./AnisoDiffuse",\
                  dataDir+"/components_Connected.raw", sizeX, sizeY, sizeZ,\
                  dataDir+"/Aniso_Diffused.raw"])
timeElapsed = time.time() - timeBefore
print "Step #3 completed in %d seconds\n\n" % timeElapsed

cmd = "./GradientVecField %s %s %s %s %s" % \
  (dataDir+"/Aniso_Diffused.raw", sizeX, sizeY, sizeZ, dataDir+"/out.vec")
print "Step #4, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./GradientVecField",\
                  dataDir+"/Aniso_Diffused.raw", sizeX, sizeY, sizeZ,\
                  dataDir+"/out.vec"])
timeElapsed = time.time() - timeBefore
print "Step #4 completed in %d seconds\n\n" % timeElapsed
cmd = "./Integratedskel %s %s %s %s 0.05 %s %s" % \
  (dataDir+"/out.vec", sizeX, sizeY, sizeZ, dataDir+"/out.seed",\
   dataDir+"/out.skel")
print "Step #5, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./Integratedskel",\
                  dataDir+"/out.vec", sizeX, sizeY, sizeZ, "0.05",\
                  dataDir+"/out.seed", dataDir+"/out.skel"])
timeElapsed = time.time() - timeBefore
print "Step #5 completed in %d seconds\n\n" % timeElapsed
cmd = "./MinSpanTree %s %s %s %s %s %s 12 4 70 0.70 %s %s %s" % \
  (dataDir+"/", "out.skel", "components_Connected.raw", sizeX, sizeY, sizeZ, \
   dataDir+"/Backbone.vtk", dataDir+"/out.txt", dataDir+"/Aniso_Diffused.raw")
print "Step #6, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./MinSpanTree",\
                  dataDir+"/", "out.skel", "components_Connected.raw", sizeX,\
                  sizeY, sizeZ, "12", "4", "70", "0.7", dataDir+"/Backbone.vtk",\
                  dataDir+"/out.txt", dataDir+"/Aniso_Diffused.raw"])
timeElapsed = time.time() - timeBefore
print "Step #6 completed in %d seconds\n\n" % timeElapsed
cmd = "./SpineMinSpanTree %s %s %s %s %s %s 12 4 70 0.70 %s %s %s" % \
  (dataDir+"/", "out.skel", "components_Connected.raw", sizeX, sizeY, sizeZ, \
   dataDir+"/Spine.vtk", dataDir+"/out.txt", dataDir+"/Aniso_Diffused.raw")
print "Step #7, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./SpineMinSpanTree",\
                  dataDir+"/", "out.skel", "components_Connected.raw", sizeX,\
                  sizeY, sizeZ, "12", "4", "70", "0.7", dataDir+"/Spine.vtk",\
                  dataDir+"/out.txt", dataDir+"/Aniso_Diffused.raw"])
timeElapsed = time.time() - timeBefore
print "Step #7 completed in %d seconds\n\n" % timeElapsed
cmd = "./BSplineFitting.exe %s %s %s %s %s %s %s %s" % \
  (dataDir+"/", "Aniso_Diffused.raw", "Backbone.vtk", sizeX, sizeY, sizeZ, \
   dataDir+"/SmoothBacbone.vtk", dataDir+"/ExtraSpine.vtk")
print "Step #8, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./BSplineFitting",\
                  dataDir+"/", "Aniso_Diffused.raw", "Backbone.vtk", sizeX,sizeY,sizeZ,\
                  dataDir+"/SmoothBacbone.vtk", dataDir+"/ExtraSpine.vtk"])
timeElapsed = time.time() - timeBefore
print "Step #8 completed in %d seconds\n\n" % timeElapsed
