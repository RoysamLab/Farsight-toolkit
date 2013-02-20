#!/usr/bin/env python

import os
import sys

if len(sys.argv) < 3:
  print sys.argv[0] + " <MDLParameters.xml> <input image 1> [input image 2...]"
  sys.exit(1)

paramFileName = sys.argv[1]

for i in range(2, len(sys.argv)):
  cmd = "./MDL2Wizard %s %s" % (paramFileName, sys.argv[i])
  print "Running this command: " + cmd
  os.system(cmd)
