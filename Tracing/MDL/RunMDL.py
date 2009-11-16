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
fullPathToOutput = parameterValues["Output image"]
dataDir = os.path.dirname(fullPathToInput)
if dataDir == "":
  dataDir = "."
filename = os.path.basename(fullPathToInput)

rawInput = False
if fullPathToInput.endswith(".raw"):
  rawInput = True
  match = re.compile(r".*?\.(\d+)x(\d+)x(\d+)\.").search(filename)
  sizeX = match.group(1)
  sizeY = match.group(2)
  sizeZ = match.group(3)

connectedComponentsSize = parameterValues["Connected components size"]

#change to the binary directory before we start calling executables
os.chdir("C:/Farsight/Bin/Tracing/MDL/release")

timeBefore = time.time()
if rawInput:
  cmd = "./volumeProcess %s %s %s %s %s" % \
    (fullPathToInput, sizeX, sizeY, sizeZ, dataDir+"/volume_Processed.raw")
  print "Step #1, executing the following command:\n%s" % cmd
  subprocess.call(["./volumeProcess",\
                    fullPathToInput, sizeX, sizeY, sizeZ, \
                    dataDir+"/volume_Processed.raw"])
  timeElapsed = time.time() - timeBefore
  print "Step #1 completed in %d seconds\n\n" % timeElapsed
else:
  cmd = "./volumeProcess %s %s " % \
    (fullPathToInput, dataDir+"/volume_Processed.raw")
  print "Step #1, executing the following command:\n%s" % cmd
  timeBefore = time.time()
  output = subprocess.Popen(["./volumeProcess", fullPathToInput,\
    dataDir+"/volume_Processed.raw"], stdout=subprocess.PIPE).communicate()[0]
  print output
  match = re.compile(r"input image size: \[(\d+), (\d+), (\d+)\]").search(output)
  sizeX = match.group(1)
  sizeY = match.group(2)
  sizeZ = match.group(3)
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

outBaseName = os.path.basename(fullPathToOutput);
backbonePath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".vtk")] + "-backbone.vtk"
spinesPath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".vtk")] + "-spines.vtk"
cmd = "./MinSpanTree %s %s %s %s %s %s 10 10 50 0.70 %s %s %s" % \
  (dataDir+"/", "out.skel", "components_Connected.raw", sizeX, sizeY, sizeZ, \
   dataDir+"/Aniso_Diffused.raw", backbonePath, spinesPath)
print "Step #6, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./MinSpanTree",\
                  dataDir+"/", "out.skel", "components_Connected.raw", sizeX,\
                  sizeY, sizeZ, "10", "10", "50", "0.7",
                  dataDir+"/Aniso_Diffused.raw", backbonePath, spinesPath])
timeElapsed = time.time() - timeBefore
print "Step #6 completed in %d seconds\n\n" % timeElapsed

cmd = "./BSplineFitting %s %s %s %s %s %s %s/smoothBB.vtk %s/smoothSpine.vtk" %\
  (dataDir+"/", filename, os.path.basename(backbonePath), sizeX, sizeY, sizeZ, dataDir, dataDir) 
print "Step #7, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./BSplineFitting",\
                dataDir+"/", filename, os.path.basename(backbonePath), sizeX,\
                sizeY, sizeZ, dataDir+"/smoothBB.vtk", dataDir+"/smoothSpine.vtk"])
timeElapsed = time.time() - timeBefore
print "Step #7 completed in %d seconds\n\n" % timeElapsed

cmd = "./AddSpinesToSmoothedBackbone %s/smoothBB.vtk %s %s" % \
  (dataDir, spinesPath, fullPathToOutput)
print "Step #8, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./AddSpinesToSmoothedBackbone",\
                  dataDir+"/smoothBB.vtk", spinesPath, fullPathToOutput])
timeElapsed = time.time() - timeBefore
print "Step #8 completed in %d seconds\n\n" % timeElapsed


outBaseName = os.path.basename(fullPathToOutput);
RefineSkelPath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".Skel")] + "RefineSkel.Skel"

cmd = "./RefiningSkeleton %s %s %s %s %s 0" % \
  (dataDir+"/", "smoothBB.vtk", "smoothSpine.vtk", "-Spine.vtk", \
   RefineSkelPath)
print "Step #9, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./RefiningSkeleton",\
                  dataDir+"/", "smoothBB.vtk", "smoothSpine.vtk", "-spines.vtk",\
                  RefineSkelPath,"0"])
timeElapsed = time.time() - timeBefore
print "Step #9 completed in %d seconds\n\n" % timeElapsed

outBaseName = os.path.basename(fullPathToOutput);
backbonePath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".vtk")] + "Refinebackbone.vtk"
spinesPath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".vtk")] + "Refinespines.vtk"
cmd = "./MinSpanTree %s %s %s %s %s %s 10 10 50 0.70 %s %s %s" % \
  (dataDir+"/", "RefineSkel.Skel", "components_Connected.raw", sizeX, sizeY, sizeZ, \
   dataDir+"/Aniso_Diffused.raw", backbonePath, spinesPath)
print "Step #10, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./MinSpanTree",\
                  dataDir+"/", "RefineSkel.Skel", "components_Connected.raw", sizeX,\
                  sizeY, sizeZ, "10", "10", "50", "0.7",
                  dataDir+"/Aniso_Diffused.raw", backbonePath, spinesPath])
timeElapsed = time.time() - timeBefore
print "Step #10 completed in %d seconds\n\n" % timeElapsed

outBaseName = os.path.basename(fullPathToOutput);
backbonePath = dataDir + "/" + outBaseName[0:outBaseName.rfind(".vtk")] + "Refinebackbone.vtk"

cmd = "./BSplineFitting %s %s %s %s %s %s %s/smoothBB-2.vtk %s/smoothSpine-2.vtk" %\
  (dataDir+"/", filename, os.path.basename(backbonePath), sizeX, sizeY, sizeZ, dataDir, dataDir) 
print "Step #11, executing the following command:\n%s" % cmd
timeBefore = time.time()
subprocess.call(["./BSplineFitting",\
                dataDir+"/", filename, os.path.basename(backbonePath), sizeX,\
                sizeY, sizeZ, dataDir+"/smoothBB-2.vtk", dataDir+"/smoothSpine-2.vtk"])
timeElapsed = time.time() - timeBefore
print "Step #11 completed in %d seconds\n\n" % timeElapsed

                 
