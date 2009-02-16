#
#  Convenience wrapper of the SmoothingRecursiveGaussianImageFilter
#

from InsightToolkit import *

from sys import argv
import os.path

if len(sys.argv) < 5:
  print "filter <input file> <output file> <output min> <output max> <sigma>"
  sys.exit(1)

inputFile = sys.argv[1]

if not os.path.exists(inputFile):
  print "Error: %s does not exist" % inputFile
  sys.exit(1)

outputFile = sys.argv[2]
outputMin = sys.argv[3]
outputMax = sys.argv[4]
sigma = sys.argv[5]

reader = itkImageFileReaderF2_New()
writer = itkImageFileWriterUS2_New()

outputCast = itkRescaleIntensityImageFilterF2US2_New()

filter  = itkSmoothingRecursiveGaussianImageFilterF2F2_New()

filter.SetInput(      reader.GetOutput()   )
outputCast.SetInput(  filter.GetOutput()      )
writer.SetInput(      outputCast.GetOutput()  )

reader.SetFileName( inputFile )
writer.SetFileName( outputFile )

outputCast.SetOutputMinimum( outputMin )
outputCast.SetOutputMaximum( outputMax )

filter.SetSigma( sigma )

writer.Update()


